# MIT License
# Copyright (c) Albert Pamonag Engineering Consultancy

"""
FRP Flexural Strengthening — ACI 440.2R-17 Chapter 10
=====================================================

Computes flexural capacity of RC beams strengthened with
externally bonded FRP systems using strain compatibility.

Reference: ACI 440.2R-17, ACI 318M-14
"""

import math


# ──────────────────────────────────────────────
# ACI 440.2R-17 Table 9.4 — Environmental Reduction Factor
# ──────────────────────────────────────────────
ENV_REDUCTION = {
    "interior":   {"carbon": 0.95, "glass": 0.75, "aramid": 0.85},
    "exterior":   {"carbon": 0.85, "glass": 0.65, "aramid": 0.75},
    "aggressive": {"carbon": 0.85, "glass": 0.50, "aramid": 0.70},
}

# ACI 440.2R-17 Table 10.2.9 — Sustained plus cyclic service load limits
CREEP_RUPTURE_LIMIT = {"glass": 0.20, "aramid": 0.30, "carbon": 0.55}


def get_env_reduction(exposure, fibertype):
    """Environmental reduction factor CE per ACI 440.2R-17 Table 9.4."""
    row = ENV_REDUCTION.get(exposure)
    if row is None:
        raise ValueError(f"Unknown exposure: {exposure}")
    ce = row.get(fibertype)
    if ce is None:
        raise ValueError(f"Unknown fibertype: {fibertype}")
    return ce


# ──────────────────────────────────────────────
# Private helpers
# ──────────────────────────────────────────────
def _beta1_aci(fc):
    """Whitney stress block factor per ACI 318M-14 Section 22.2.2.4.3."""
    if fc <= 28:
        return 0.85
    return max(0.85 - ((fc - 28) / 7) * 0.05, 0.65)


def _neutral_axis_solver(b, d, df, As, Af, fc, fy, Ef, Ec, Es, eps_bi, eps_fd):
    """Iterative solver for compression block depth c (20 iterations).

    Uses Hognestad parabolic stress block parameters (alpha1, beta1)
    derived from concrete strain at each iteration.
    ACI 440.2R-17 Section 10.2.10.
    """
    iterations = []
    c = 0.20 * d  # initial guess

    for i in range(1, 21):
        # FRP effective strain (Eq 10.2.3)
        eps_fe = min(0.003 * ((df - c) / c) - eps_bi, eps_fd)
        # Concrete strain at top fiber
        eps_c = (eps_fe + eps_bi) * (c / (df - c))
        # Steel strain
        eps_s = (eps_fe + eps_bi) * ((d - c) / (df - c))
        # Steel stress (capped at yield)
        fs = min(Es * eps_s, fy)
        # Hognestad peak strain
        eps_c_prime = (1.7 * fc) / Ec
        # Stress block parameters from Hognestad model
        b1 = ((4 * eps_c_prime) - eps_c) / ((6 * eps_c_prime) - (2 * eps_c))
        a1 = ((3 * eps_c_prime * eps_c) - eps_c ** 2) / (3 * b1 * eps_c_prime ** 2)
        # FRP effective stress
        ffe = Ef * eps_fe
        # Updated neutral axis from force equilibrium
        c = ((As * fs) + (Af * ffe)) / (a1 * fc * b1 * b)

        iterations.append({
            "i": i, "c": c, "eps_fe": eps_fe, "eps_c": eps_c,
            "eps_s": eps_s, "fs": fs, "eps_c_prime": eps_c_prime,
            "b1": b1, "a1": a1, "ffe": ffe,
        })

    last = iterations[-1]
    return {
        "c": c,
        "beta1": last["b1"],
        "alpha1": last["a1"],
        "fs": last["fs"],
        "ffe": last["ffe"],
        "eps_fe": last["eps_fe"],
        "eps_c": last["eps_c"],
        "eps_s": last["eps_s"],
        "iterations": iterations,
    }


def _compute_service_stress_k(rho_s, rho_f, Es, Ef, Ec, df, d):
    """Neutral axis ratio k for FRP-strengthened section under service loads."""
    t1 = (rho_s * (Es / Ec) + rho_f * (Ef / Ec)) ** 2
    t2 = 2 * (rho_s * (Es / Ec) + rho_f * (Ef / Ec) * (df / d))
    t3 = rho_s * (Es / Ec) + rho_f * (Ef / Ec)
    return math.sqrt(t1 + t2) - t3


def _compute_fss(Ms_Nmm, eps_bi, Af, Ef, df, k, d, Es, As):
    """Steel stress under service loads for FRP-strengthened section."""
    num = (Ms_Nmm + eps_bi * Af * Ef * (df - (k * d) / 3)) * (d - k * d) * Es
    den = (As * Es * (d - (k * d) / 3) * (d - k * d)) + \
          (Af * Ef * (df - (k * d) / 3) * (df - k * d))
    return num / den


def _compute_ffs(fss, Ef, Es, df, k, d, eps_bi):
    """FRP stress under service loads (ACI 440.2R-17 Section 10.2.9)."""
    return fss * (Ef / Es) * ((df - k * d) / (d - k * d)) - eps_bi * Ef


# ──────────────────────────────────────────────
# Main public function
# ──────────────────────────────────────────────
def frp_flexural_strengthening(
    h, b, d, df,
    As, fy, Es,
    fc,
    n_ply, thk_ply, Ef,
    CE, ffu_star, eps_fu_star,
    fibertype,
    moment_dead, moment_live, moment_capacity,
):
    """
    FRP flexural strengthening design per ACI 440.2R-17 Chapter 10.

    Parameters
    ----------
    h : float           Beam total height (mm).
    b : float           Beam width (mm).
    d : float           Effective depth to steel centroid (mm).
    df : float          Effective depth to FRP centroid (mm).
    As : float          Total steel reinforcement area (mm^2).
    fy : float          Steel yield strength (MPa).
    Es : float          Steel modulus of elasticity (MPa).
    fc : float          Concrete compressive strength (MPa).
    n_ply : int         Number of FRP plies.
    thk_ply : float     Thickness per ply (mm).
    Ef : float          FRP modulus of elasticity (MPa).
    CE : float          Environmental reduction factor (Table 9.4).
    ffu_star : float    Manufacturer-reported FRP tensile strength (MPa).
    eps_fu_star : float Manufacturer-reported FRP rupture strain (mm/mm).
    fibertype : str     "carbon", "glass", or "aramid".
    moment_dead : float Dead load moment (kN-m).
    moment_live : float Live load moment (kN-m).
    moment_capacity : float  Required design moment for utilization (kN-m).

    Returns
    -------
    dict
        All intermediate values, iteration history, and pass/fail status.
    """
    # 1.0 Design FRP properties (Section 9.4)
    ffu = CE * ffu_star                     # MPa, design tensile strength
    eps_fu = CE * eps_fu_star               # design rupture strain

    # 2.0 Section properties
    Ec = 4700 * math.sqrt(fc)               # MPa, concrete modulus (ACI 318)
    Af = n_ply * thk_ply * b                # mm^2, FRP area
    rho_s = As / (b * d)                    # steel reinforcement ratio
    rho_f = Af / (b * h)                    # FRP reinforcement ratio
    n_ratio = Es / Ec                       # modular ratio

    # 3.0 Existing beam capacity without FRP
    beta1 = _beta1_aci(fc)                  # Whitney block factor
    c_wo = (As * fy) / (0.85 * fc * b)     # mm, NA depth without FRP
    a_depth = beta1 * c_wo                  # mm, equivalent block depth
    Mn_existing = 0.9 * As * fy * (d - a_depth / 2) / 1e6  # kN-m

    # 4.0 Cracked section properties
    k = math.sqrt(2 * rho_s * n_ratio + (rho_s * n_ratio) ** 2) - rho_s * n_ratio
    j = 1 - k / 3
    Icr = (b * k ** 2 * j * d ** 3) / 2    # mm^4, cracked moment of inertia

    # 5.0 Strain at soffit at time of FRP installation
    MDL_Nmm = moment_dead * 1e6             # N-mm
    sigma_bi = MDL_Nmm * (df - k * d) / Icr  # MPa, stress at soffit
    eps_bi = sigma_bi / Ec                  # initial substrate strain

    # 6.0 Debonding strain (Eq 10.1.1)
    eps_fd = min(
        0.41 * math.sqrt(fc / (n_ply * Ef * thk_ply)),
        0.90 * eps_fu,
    )

    # 7.0 Load combinations
    Ms = moment_dead + moment_live          # kN-m, service moment
    Mu = 1.2 * moment_dead + 1.6 * moment_live  # kN-m, factored moment
    Mlimit = 1.1 * moment_dead + 0.75 * moment_live  # kN-m, strengthening limit

    # 8.0 Iterative neutral axis solution (Section 10.2.10)
    solver = _neutral_axis_solver(b, d, df, As, Af, fc, fy, Ef, Ec, Es, eps_bi, eps_fd)
    c_final = solver["c"]
    b1_final = solver["beta1"]
    fs_final = solver["fs"]
    ffe_final = solver["ffe"]
    eps_fe_final = solver["eps_fe"]

    # 9.0 Flexural strength (Eq 10.2.10d)
    Mns = As * fs_final * (d - (b1_final * c_final) / 2)       # N-mm, steel contribution
    Mnf = Af * ffe_final * (df - (b1_final * c_final) / 2)     # N-mm, FRP contribution
    Mns_kNm = Mns / 1e6
    Mnf_kNm = Mnf / 1e6
    phi_Mn = 0.90 * (Mns_kNm + 0.85 * Mnf_kNm)                # kN-m, design capacity

    # Utilization
    flexure_ratio = moment_capacity / phi_Mn if phi_Mn > 0 else 0.0

    # Strength increase over existing beam
    if Mn_existing > 0:
        strength_increase = ((phi_Mn - Mn_existing) / Mn_existing) * 100
    else:
        strength_increase = 0.0

    # 10.0 Service stress checks (Section 10.2.9)
    k_service = _compute_service_stress_k(rho_s, rho_f, Es, Ef, Ec, h, d)
    Ms_Nmm = Ms * 1e6
    fss = _compute_fss(Ms_Nmm, eps_bi, Af, Ef, h, k_service, d, Es, As)
    fss_limit = 0.8 * fy                   # MPa, steel service limit
    fss_ratio = fss / fss_limit

    ffs = _compute_ffs(fss, Ef, Es, h, k_service, d, eps_bi)
    creep_limit = CREEP_RUPTURE_LIMIT[fibertype]
    ffs_limit = creep_limit * ffu           # MPa, FRP creep-rupture limit
    ffs_ratio = ffs / ffs_limit if ffs_limit > 0 else 0.0

    # 11.0 Pass/fail
    flexure_pass = flexure_ratio <= 1.0
    fss_pass = fss_ratio <= 1.0
    ffs_pass = ffs_ratio <= 1.0
    overall_pass = flexure_pass and fss_pass and ffs_pass

    return {
        # Material design properties
        "CE": CE,
        "ffu": ffu,
        "eps_fu": eps_fu,
        "Ec": Ec,
        "beta1": beta1,
        # Section properties
        "d": d,
        "df": df,
        "As": As,
        "Af": Af,
        "rho_s": rho_s,
        "rho_f": rho_f,
        "n_ratio": n_ratio,
        # Existing beam
        "c_wo_frp": c_wo,
        "a_depth": a_depth,
        "Mn_existing": Mn_existing,
        # Cracked section
        "k": k,
        "j": j,
        "Icr": Icr,
        "sigma_bi": sigma_bi,
        "eps_bi": eps_bi,
        # FRP strain limits
        "eps_fd": eps_fd,
        # Load combinations
        "MDL": moment_dead,
        "Ms": Ms,
        "Mu": Mu,
        "Mlimit": Mlimit,
        # Neutral axis solver
        "c_final": c_final,
        "b1_final": b1_final,
        "fs_final": fs_final,
        "ffe_final": ffe_final,
        "eps_fe_final": eps_fe_final,
        "iterations": solver["iterations"],
        # Flexural strength
        "Mns_kNm": Mns_kNm,
        "Mnf_kNm": Mnf_kNm,
        "phi_Mn": phi_Mn,
        "flexure_ratio": flexure_ratio,
        "strength_increase": strength_increase,
        # Service stresses
        "k_service": k_service,
        "fss": fss,
        "fss_ratio": fss_ratio,
        "fss_limit": fss_limit,
        "ffs": ffs,
        "ffs_ratio": ffs_ratio,
        "creep_limit": creep_limit,
        "ffs_limit": ffs_limit,
        # Status
        "flexure_pass": flexure_pass,
        "fss_pass": fss_pass,
        "ffs_pass": ffs_pass,
        "overall_pass": overall_pass,
    }
