# MIT License
# Copyright (c) Albert Pamonag Engineering Consultancy

"""
Column P-M Interaction Diagram Generator
==========================================

Generates axial load vs moment interaction data
for rectangular RC columns.

Reference: NSCP 2015 / ACI 318-19
"""

import math


def _beta_one(fc):
    """Whitney stress block factor."""
    if fc <= 28:
        return 0.85
    elif fc < 55:
        return max(0.85 - 0.05 * (fc - 28) / 7, 0.65)
    else:
        return 0.65


def generate_interaction_diagram(
    fc, fy, b, h, n_bars, d_bar, n_bars_side=0,
    cover=40, c_start=None, c_end=80, c_step=5,
):
    """
    Generate P-M interaction diagram data for a rectangular column.

    Parameters
    ----------
    fc : float
        Concrete compressive strength (MPa).
    fy : float
        Steel yield strength (MPa).
    b : float
        Column width (mm).
    h : float
        Column depth (mm).
    n_bars : int
        Total number of bars.
    d_bar : float
        Main bar diameter (mm).
    n_bars_side : int
        Number of bars on each side (tension/compression face).
        Defaults to n_bars // 4.
    cover : float
        Concrete cover to tie center (mm).
    c_start : float or None
        Starting neutral axis depth (mm). Default is h + buffer.
    c_end : float
        Ending neutral axis depth (mm).
    c_step : float
        Step size for neutral axis (mm).

    Returns
    -------
    dict
        Keys: points (list of dicts), pure_compression, balanced_point
    """
    es_mod = 200000.0
    ecu = 0.003

    # Tie diameter
    tie = 10 if d_bar <= 32 else 12

    # Effective depths
    d_eff = h - cover - tie - d_bar / 2.0
    d_prime = cover + tie + d_bar / 2.0

    beta1 = _beta_one(fc)

    if n_bars_side == 0:
        n_bars_side = max(n_bars // 4, 2)

    # Steel areas on each face
    ast = math.pi * d_bar ** 2 / 4.0 * n_bars_side
    asc = ast

    ey = fy / es_mod

    if c_start is None:
        c_start = int(h * 1.5)

    points = []
    balanced_point = None

    c_values = list(range(int(c_start), int(c_end), -int(c_step)))

    for c in c_values:
        # Strains
        es_tension = ecu * (d_eff - c) / c
        es_comp = ecu * (c - d_prime) / c
        a = beta1 * c

        # Stresses (capped at fy)
        fs = min(abs(es_tension) * es_mod, fy) * (1 if es_tension >= 0 else -1)
        fs_prime = min(abs(es_comp) * es_mod, fy) * (1 if es_comp >= 0 else -1)

        # Strength reduction factor
        if abs(es_tension) >= 0.005:
            phi = 0.90
            classify = "tension-controlled"
        elif abs(es_tension) <= ey:
            phi = 0.65
            classify = "compression-controlled"
        else:
            phi = 0.65 + 0.25 * (abs(es_tension) - ey) / (0.005 - ey)
            classify = "transition"

        # Forces
        t_force = ast * min(fy, abs(fs)) / 1000.0
        cc = 0.85 * fc * a * b / 1000.0
        cs = asc * (min(fy, abs(fs_prime)) - 0.85 * fc) / 1000.0

        pn = cc + cs - t_force
        pu = phi * pn

        # Moment about centroid of tension steel
        if abs(pn) > 0.1:
            ecc = (cc * (d_eff - a / 2.0) + cs * (d_eff - d_prime)) / pn - (h / 2.0 - d_prime)
        else:
            ecc = 0

        if ecc <= 1.5 * h:
            mn = pn * ecc / 1000.0
            mu = phi * mn

            point = {
                "c": c,
                "pn": round(pn, 2),
                "pu": round(pu, 2),
                "mn": round(mn, 2),
                "mu": round(mu, 2),
                "phi": round(phi, 4),
                "ecc": round(ecc, 2),
                "classification": classify,
                "es_tension": round(es_tension, 6),
            }
            points.append(point)

            if balanced_point is None and abs(es_tension) >= ey:
                balanced_point = point

    # Pure compression
    ag = b * h
    ast_total = math.pi * d_bar ** 2 / 4.0 * n_bars
    po = 0.85 * fc * (ag - ast_total) + fy * ast_total
    po_kn = po / 1000.0

    return {
        "points": points,
        "pure_compression_kn": round(po_kn, 2),
        "balanced_point": balanced_point,
        "beta1": round(beta1, 4),
        "d_eff": round(d_eff, 2),
        "d_prime": round(d_prime, 2),
    }
