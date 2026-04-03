# MIT License
# Copyright (c) Albert Pamonag Engineering Consultancy

"""
Mander's Confined Concrete Model
==================================

Computes confined concrete strength and stress-strain parameters.

Reference: Mander, Priestley, and Park (1988)
"""

import math


def confined_strength_ratio(fpco, fpl):
    """
    Compute confined concrete strength ratio f'cc / f'co.

    Parameters
    ----------
    fpco : float
        Unconfined concrete compressive strength (MPa).
    fpl : float
        Effective lateral confining stress (MPa).

    Returns
    -------
    float
        Ratio f'cc / f'co.
    """
    if fpco <= 0:
        raise ValueError("fpco must be positive.")
    ratio = fpl / fpco
    return -1.254 + 2.254 * math.sqrt(1 + 7.94 * ratio) - 2 * ratio


def confined_stress_strain(
    fc, fy_transverse, b, h, cover, db_main, db_tie,
    n_bars_x, n_bars_y, s_tie, n_legs_x=2, n_legs_y=2,
):
    """
    Compute confined concrete strength and strain using Mander's model.

    Parameters
    ----------
    fc : float
        Unconfined compressive strength (MPa).
    fy_transverse : float
        Yield strength of transverse reinforcement (MPa).
    b : float
        Section width (mm).
    h : float
        Section depth (mm).
    cover : float
        Clear cover to tie center (mm).
    db_main : float
        Diameter of main bars (mm).
    db_tie : float
        Diameter of ties (mm).
    n_bars_x : int
        Number of bars along x-direction.
    n_bars_y : int
        Number of bars along y-direction.
    s_tie : float
        Spacing of transverse reinforcement (mm).
    n_legs_x : int
        Number of tie legs in x-direction.
    n_legs_y : int
        Number of tie legs in y-direction.

    Returns
    -------
    dict
        Keys: fcc, ecc, bc, dc, ke, fl_x, fl_y, fl_eff, ratio
    """
    # Core dimensions (center-to-center of ties)
    bc = b - 2 * cover + db_tie
    dc = h - 2 * cover + db_tie

    # Area of transverse steel per leg
    as_tie = math.pi * db_tie ** 2 / 4.0

    # Transverse steel ratios
    rho_x = n_legs_x * as_tie / (dc * s_tie)
    rho_y = n_legs_y * as_tie / (bc * s_tie)

    # Clear spacing between bars along each face
    n_total = 2 * (n_bars_x + n_bars_y) - 4
    as_main = math.pi * db_main ** 2 / 4.0
    rho_cc = n_total * as_main / (bc * dc)

    # Clear spacings between bars
    wi_x = []
    if n_bars_x > 1:
        spacing_x = (bc - db_main) / (n_bars_x - 1)
        for _ in range(n_bars_x - 1):
            wi_x.append(spacing_x - db_main)

    wi_y = []
    if n_bars_y > 1:
        spacing_y = (dc - db_main) / (n_bars_y - 1)
        for _ in range(n_bars_y - 1):
            wi_y.append(spacing_y - db_main)

    sum_wi2_x = sum(w ** 2 for w in wi_x)
    sum_wi2_y = sum(w ** 2 for w in wi_y)
    sum_wi2 = sum_wi2_x + sum_wi2_y

    # Core area and effective concrete core
    ac = bc * dc
    acc = ac * (1 - rho_cc)

    # Confinement effectiveness
    s_prime = s_tie - db_tie
    ae_factor1 = 1 - sum_wi2 / (6 * bc * dc)
    ae_factor2 = (1 - s_prime / (2 * bc))
    ae_factor3 = (1 - s_prime / (2 * dc))
    ae = bc * dc * ae_factor1 * ae_factor2 * ae_factor3
    ke = ae / acc

    # Total transverse steel area per direction
    asy = n_legs_x * as_tie
    asz = n_legs_y * as_tie

    # Effective lateral confining stresses
    fl_x = ke * rho_x * fy_transverse
    fl_y = ke * rho_y * fy_transverse

    # Effective lateral stress (biaxial)
    fl_eff = (fl_x + fl_y) / 2.0

    # Confined strength
    ratio = confined_strength_ratio(fc, fl_eff)
    fcc = ratio * fc

    # Confined strain
    eco = 0.002
    ecc = eco * (1 + 5 * (ratio - 1))

    # Secant modulus at peak
    e_sec = fcc / ecc
    ec = 5000 * math.sqrt(fc)

    # Mander stress-strain curve data
    r = ec / (ec - e_sec)
    ecu = 0.004 + 1.4 * ((fl_x + fl_y) / 2) * 0.01 / fcc if fcc > 0 else 0.004
    ecu = max(ecu, 3 * ecc)
    n_pts = 200
    strain_pts = [i * ecu / n_pts for i in range(n_pts + 1)]
    stress_pts = []
    for eps in strain_pts:
        x = eps / ecc if ecc > 0 else 0
        denom = r - 1 + x ** r
        stress = fcc * x * r / denom if denom != 0 else 0
        stress_pts.append(round(stress, 4))
    strain_pts = [round(e, 8) for e in strain_pts]

    return {
        "fcc": round(fcc, 2),
        "ecc": round(ecc, 6),
        "bc": round(bc, 2),
        "dc": round(dc, 2),
        "ac": round(ac, 2),
        "acc": round(acc, 2),
        "ae": round(ae, 2),
        "ke": round(ke, 4),
        "s_prime": round(s_prime, 2),
        "fl_x": round(fl_x, 4),
        "fl_y": round(fl_y, 4),
        "fl_eff": round(fl_eff, 4),
        "ratio": round(ratio, 4),
        "rho_x": round(rho_x, 6),
        "rho_y": round(rho_y, 6),
        "rho_cc": round(rho_cc, 6),
        "n_total": n_total,
        "as_main": round(as_main, 2),
        "as_tie": round(as_tie, 2),
        "asy": round(asy, 2),
        "asz": round(asz, 2),
        "wi_x": [round(w, 2) for w in wi_x],
        "wi_y": [round(w, 2) for w in wi_y],
        "sum_wi2_x": round(sum_wi2_x, 2),
        "sum_wi2_y": round(sum_wi2_y, 2),
        "eco": eco,
        "e_sec": round(e_sec, 2),
        "ec": round(ec, 2),
        "ecu": round(ecu, 6),
        "r": round(r, 4),
        "curve_strain": strain_pts,
        "curve_stress": stress_pts,
    }
