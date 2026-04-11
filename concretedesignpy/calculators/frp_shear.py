# MIT License
# Copyright (c) Albert Pamonag Engineering Consultancy

"""
FRP Shear Strengthening — ACI 440.2R-17 Chapter 11
===================================================

Computes shear contribution of externally bonded FRP wraps
for beam shear strengthening.

Reference: ACI 440.2R-17 Section 11.4
"""

import math


def frp_shear_strengthening(
    n, tf, wf, sf, Ef,
    eps_fu_star, CE, alpha,
    fc, dfv, wrap,
    shear_demand=0.0, shear_conc=0.0, shear_steel=0.0,
):
    """
    FRP shear strengthening design per ACI 440.2R-17 Chapter 11.

    Parameters
    ----------
    n : int             Number of FRP plies.
    tf : float          Thickness per ply (mm).
    wf : float          Width of FRP wrap strip (mm).
    sf : float          Center-to-center spacing of FRP strips (mm).
    Ef : float          FRP modulus of elasticity (MPa).
    eps_fu_star : float Manufacturer-reported FRP rupture strain (mm/mm).
    CE : float          Environmental reduction factor (Table 9.4).
    alpha : float       Angle of FRP fiber orientation (degrees).
    fc : float          Concrete compressive strength (MPa).
    dfv : float         Effective depth for shear (mm).
    wrap : str          Wrap type: "U" (U-wrap / 3-sided) or "S" (two-sided).
    shear_demand : float  Factored shear demand (kN).
    shear_conc : float  Concrete shear contribution (kN).
    shear_steel : float Steel shear contribution (kN).

    Returns
    -------
    dict
        All intermediate values and pass/fail status.
    """
    # 1.0 Design FRP rupture strain
    efu = CE * eps_fu_star                  # design rupture strain

    # 2.0 Effective bond length (Eq 11.4.1.2c)
    le = 23300 / ((n * tf * Ef) ** 0.58)    # mm

    # 3.0 Bond reduction factors (Section 11.4.1.2)
    k1 = (fc / 27) ** (2 / 3)

    if wrap == "U":
        k2 = (dfv - 2 * le) / dfv          # U-wrap (3-sided)
    else:
        k2 = (dfv - le) / dfv              # Two-sided bonding

    # 4.0 Bond-dependent coefficient kv (Eq 11.4.1.2b, capped at 0.75)
    kv_raw = (k1 * k2 * le) / (11900 * efu)
    kv = min(kv_raw, 0.75)

    # 5.0 Effective FRP strain (max 0.004 per Section 11.4.1.2)
    eps_fe_raw = kv * efu
    eps_fe = min(eps_fe_raw, 0.004)

    # 6.0 Effective FRP tensile stress
    ffe = eps_fe * Ef                       # MPa

    # 7.0 Area of FRP shear reinforcement (Eq 11.4.1.1)
    afv = 2 * n * tf * wf                   # mm^2

    # 8.0 FRP shear contribution Vf (Eq 11.4.1.1)
    alpha_rad = math.radians(alpha)
    Vf_N = (afv * ffe * (math.sin(alpha_rad) + math.cos(alpha_rad)) * dfv) / sf
    Vf = Vf_N / 1000                        # kN

    # 9.0 Design FRP shear strength
    psi = 0.95 if wrap == "U" else 0.85     # additional reduction (Section 11.4)
    phi = 0.75                               # strength reduction factor
    phi_Vnf = phi * psi * Vf                 # kN

    # 10.0 Total shear capacity and utilization
    total_capacity = phi_Vnf + shear_conc + shear_steel  # kN
    ratio = shear_demand / total_capacity if total_capacity > 0 else 0.0
    spacing = sf - wf                        # mm, clear spacing

    overall_pass = ratio <= 1.0

    return {
        # Bond properties
        "le": le,
        "k1": k1,
        "k2": k2,
        # Design FRP strain
        "efu": efu,
        "kv_raw": kv_raw,
        "kv": kv,
        "eps_fe_raw": eps_fe_raw,
        "eps_fe": eps_fe,
        # FRP contribution
        "ffe": ffe,
        "afv": afv,
        "Vf": Vf,
        # Design strength
        "psi": psi,
        "phi": phi,
        "phi_Vnf": phi_Vnf,
        # Total capacity
        "total_capacity": total_capacity,
        "ratio": ratio,
        "spacing": spacing,
        # Pass/fail
        "overall_pass": overall_pass,
        # Echoed demands
        "shear_demand": shear_demand,
        "shear_conc": shear_conc,
        "shear_steel": shear_steel,
    }
