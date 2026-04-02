# MIT License
# Copyright (c) Albert Pamonag Engineering Consultancy

"""
FRP Flexural Strengthening Calculator
=======================================

Computes flexural strengthening of RC beams with externally
bonded FRP reinforcement per ACI 440R-17.

Reference: ACI 440R-17, NSCP 2015
"""

import math


_ERF_TABLE = {
    "interior": {"carbon": 0.95, "glass": 0.75, "aramid": 0.85},
    "exterior": {"carbon": 0.85, "glass": 0.65, "aramid": 0.75},
    "aggressive": {"carbon": 0.85, "glass": 0.50, "aramid": 0.70},
}


def _environmental_reduction(exposure, fiber):
    """Environmental reduction factor CE per ACI 440 Table 9.4."""
    if exposure not in _ERF_TABLE:
        raise ValueError(f"Unknown exposure: {exposure}")
    if fiber not in _ERF_TABLE[exposure]:
        raise ValueError(f"Unknown fiber: {fiber}")
    return _ERF_TABLE[exposure][fiber]


def _beta_one(fc):
    """Whitney stress block factor."""
    if fc <= 28:
        return 0.85
    elif fc < 55:
        return max(0.85 - 0.05 * (fc - 28) / 7, 0.65)
    else:
        return 0.65


def frp_flexural_strengthening(
    b, d, h, fc, fy, as_steel, mu_required,
    tf, ffu_star, efu_star, ef_frp, n_plies,
    exposure="interior", fiber="carbon",
    mdl=0, k_initial=0.334, icr=None,
    es=200000.0, psi_f=0.85,
):
    """
    Compute FRP flexural strengthening of an RC beam.

    Parameters
    ----------
    b : float
        Beam width (mm).
    d : float
        Effective depth to tension steel (mm).
    h : float
        Total beam height (mm).
    fc : float
        Concrete compressive strength (MPa).
    fy : float
        Steel yield strength (MPa).
    as_steel : float
        Area of existing tension reinforcement (mm^2).
    mu_required : float
        Required factored moment capacity (kN-m).
    tf : float
        FRP thickness per ply (mm).
    ffu_star : float
        Manufacturer's ultimate tensile strength of FRP (MPa).
    efu_star : float
        Manufacturer's rupture strain of FRP (mm/mm).
    ef_frp : float
        Modulus of elasticity of FRP laminates (MPa).
    n_plies : int
        Number of FRP plies.
    exposure : str
        'interior', 'exterior', or 'aggressive'.
    fiber : str
        'carbon', 'glass', or 'aramid'.
    mdl : float
        Dead load moment for initial strain (kN-m).
    k_initial : float
        Initial neutral axis depth ratio for existing strain calc.
    icr : float or None
        Cracked moment of inertia (mm^4). Computed if None.
    es : float
        Steel modulus of elasticity (MPa).
    psi_f : float
        FRP strength reduction factor (default 0.85).

    Returns
    -------
    dict
        Keys: ffu, efu, efd, af, c, mns, mnf, phi_mn, mu_required,
        failure_mode, service_stress_steel, service_stress_frp, status
    """
    # Design material properties
    ce = _environmental_reduction(exposure, fiber)
    ffu = ce * ffu_star
    efu = ce * efu_star

    ec = 4700 * math.sqrt(fc)
    beta1 = _beta_one(fc)

    # FRP area
    af = n_plies * tf * b

    # Modular ratio
    n = max(6, math.ceil(es * 1000 / ec))

    # Cracked moment of inertia (approximate if not provided)
    if icr is None:
        a_coeff = b / 2.0
        b_coeff = n * as_steel
        c_coeff = -n * as_steel * d
        disc = b_coeff ** 2 - 4 * a_coeff * c_coeff
        kd = (-b_coeff + math.sqrt(disc)) / (2 * a_coeff)
        icr = b * kd ** 3 / 3.0 + n * as_steel * (d - kd) ** 2

    # Existing strain on soffit
    kd_initial = k_initial * d
    if mdl > 0:
        ebi = (mdl * 1e6 * (h - kd_initial)) / (icr * ec)
    else:
        ebi = 0

    # Design strain of FRP
    efd = 0.41 * math.sqrt(fc / (2 * ef_frp * tf))
    debonding = efd <= 0.9 * efu
    failure_mode = "debonding" if debonding else "FRP rupture"

    # Iterative solver for neutral axis depth c
    c_val = 0.20 * d
    for _ in range(30):
        epsilon_fe = 0.003 * ((h - c_val) / c_val) - ebi
        epsilon_fe = min(epsilon_fe, efd)

        epsilon_c = (epsilon_fe + ebi) * c_val / (h - c_val)
        epsilon_s = (epsilon_fe + ebi) * (d - c_val) / (h - c_val)

        fs = min(es * 1000 * epsilon_s, fy)
        f_fe = ef_frp * epsilon_fe

        e_c_prime = 1.7 * fc / ec
        if epsilon_c > 0 and e_c_prime > 0:
            beta_one_calc = (4 * e_c_prime - epsilon_c) / (6 * e_c_prime - 2 * epsilon_c)
            alpha_one = ((3 * e_c_prime * epsilon_c - epsilon_c ** 2)
                         / (3 * beta_one_calc * e_c_prime ** 2))
        else:
            beta_one_calc = beta1
            alpha_one = 0.85

        denom = alpha_one * fc * beta_one_calc * b
        c_val = (as_steel * fs + af * f_fe) / denom if denom > 0 else c_val

    # Final steel and FRP contributions
    mns = as_steel * fs * (d - (beta_one_calc * c_val) / 2.0)
    mnf = af * f_fe * (h - (beta_one_calc * c_val) / 2.0)

    # Strength reduction factor
    epsilon_t = (epsilon_fe + ebi) * (d - c_val) / (h - c_val)
    if epsilon_t >= 0.005:
        phi = 0.90
    elif epsilon_t <= fy / (es * 1000):
        phi = 0.65
    else:
        phi = 0.65 + 0.25 * (epsilon_t - fy / (es * 1000)) / (0.005 - fy / (es * 1000))

    phi_mn = phi * (mns + psi_f * mnf) / 1e6

    return {
        "ce": ce,
        "ffu": round(ffu, 2),
        "efu": round(efu, 6),
        "efd": round(efd, 6),
        "ebi": round(ebi, 6),
        "af": round(af, 2),
        "c": round(c_val, 2),
        "fs": round(fs, 2),
        "f_fe": round(f_fe, 2),
        "mns_knm": round(mns / 1e6, 2),
        "mnf_knm": round(mnf / 1e6, 2),
        "phi": round(phi, 4),
        "phi_mn": round(phi_mn, 2),
        "mu_required": mu_required,
        "failure_mode": failure_mode,
        "status": "PASS" if phi_mn >= mu_required else "FAIL",
    }
