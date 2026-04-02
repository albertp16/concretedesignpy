# MIT License
# Copyright (c) Albert Pamonag Engineering Consultancy

"""
Moment-Curvature Analysis
===========================

Computes the 6-point moment-curvature relationship for
rectangular RC sections.

Points:
    1. Just before cracking
    2. Just after cracking
    3. Elastic limit (0.45 f'c)
    4. Steel yields
    5. Concrete at peak strain
    6. Concrete at ultimate strain

Reference: NSCP 2015 / ACI 318-19
"""

import math


def moment_curvature_analysis(b, h, d, fc, fy, as_tension, es=200000.0):
    """
    Compute 6-point moment-curvature relationship.

    Parameters
    ----------
    b : float
        Section width (mm).
    h : float
        Section total height (mm).
    d : float
        Effective depth (mm).
    fc : float
        Concrete compressive strength (MPa).
    fy : float
        Steel yield strength (MPa).
    as_tension : float
        Area of tension reinforcement (mm^2).
    es : float
        Steel modulus of elasticity (MPa).

    Returns
    -------
    dict
        Keys: points (list of 6 dicts with phi and moment),
              section_properties
    """
    if b <= 0 or h <= 0 or d <= 0 or fc <= 0 or fy <= 0:
        raise ValueError("All section/material parameters must be positive.")

    ec = 4700 * math.sqrt(fc)
    n = es / ec
    fr = 0.62 * math.sqrt(fc)
    ig = b * h ** 3 / 12.0
    yt = h / 2.0
    ey = fy / es

    # Transformed section properties (uncracked)
    # Neutral axis of uncracked transformed section (from top)
    at = b * h + (n - 1) * as_tension
    y_bar = (b * h * h / 2.0 + (n - 1) * as_tension * d) / at
    it = (b * h ** 3 / 12.0 + b * h * (h / 2.0 - y_bar) ** 2
          + (n - 1) * as_tension * (d - y_bar) ** 2)

    # Cracked section neutral axis (from top)
    a_coeff = b / 2.0
    b_coeff = n * as_tension
    c_coeff = -n * as_tension * d
    disc = b_coeff ** 2 - 4 * a_coeff * c_coeff
    kd = (-b_coeff + math.sqrt(disc)) / (2 * a_coeff)
    icr = b * kd ** 3 / 3.0 + n * as_tension * (d - kd) ** 2

    points = []

    # Point 1: Just before cracking
    mcr = fr * it / (h - y_bar)
    phi1 = mcr / (ec * it) if it > 0 else 0
    points.append({
        "point": 1,
        "event": "Before cracking",
        "phi": phi1,
        "moment_knm": mcr / 1e6,
    })

    # Point 2: Just after cracking
    phi2 = mcr / (ec * icr) if icr > 0 else 0
    points.append({
        "point": 2,
        "event": "After cracking",
        "phi": phi2,
        "moment_knm": mcr / 1e6,
    })

    # Point 3: Elastic limit (0.45 f'c)
    fc_limit = 0.45 * fc
    ec_limit = fc_limit / ec
    c3 = kd
    m3 = 0.5 * fc_limit * b * c3 * (d - c3 / 3.0)
    phi3 = ec_limit / c3 if c3 > 0 else 0
    points.append({
        "point": 3,
        "event": "Elastic limit (0.45fc)",
        "phi": phi3,
        "moment_knm": m3 / 1e6,
    })

    # Point 4: Steel yields
    ec_yield = ey * kd / (d - kd) if (d - kd) > 0 else 0
    c4 = kd
    fs4 = fy
    cc4 = 0.5 * ec_yield * ec * b * c4
    m4 = cc4 * (d - c4 / 3.0)
    phi4 = ey / (d - c4) if (d - c4) > 0 else 0
    points.append({
        "point": 4,
        "event": "Steel yields",
        "phi": phi4,
        "moment_knm": m4 / 1e6,
    })

    # Point 5: Concrete at peak strain (epsilon_co = 0.002)
    eco = 0.002
    # Find c such that equilibrium holds: 0.85*fc*beta1*c*b = As*fs
    beta1 = 0.85 if fc <= 28 else max(0.85 - 0.05 * (fc - 28) / 7, 0.65)
    c5 = as_tension * fy / (0.85 * fc * beta1 * b)
    a5 = beta1 * c5
    m5 = as_tension * fy * (d - a5 / 2.0)
    phi5 = eco / c5 if c5 > 0 else 0
    points.append({
        "point": 5,
        "event": "Concrete peak strain",
        "phi": phi5,
        "moment_knm": m5 / 1e6,
    })

    # Point 6: Concrete at ultimate strain (epsilon_cu = 0.003)
    ecu = 0.003
    c6 = c5  # Same equilibrium for singly reinforced
    a6 = beta1 * c6
    m6 = as_tension * fy * (d - a6 / 2.0)
    phi6 = ecu / c6 if c6 > 0 else 0
    points.append({
        "point": 6,
        "event": "Concrete ultimate",
        "phi": phi6,
        "moment_knm": m6 / 1e6,
    })

    # Round all values
    for p in points:
        p["phi"] = round(p["phi"], 8)
        p["moment_knm"] = round(p["moment_knm"], 2)

    return {
        "points": points,
        "section_properties": {
            "ig": round(ig, 2),
            "it": round(it, 2),
            "icr": round(icr, 2),
            "y_bar": round(y_bar, 2),
            "kd": round(kd, 2),
            "ec": round(ec, 2),
            "fr": round(fr, 4),
            "n": round(n, 4),
            "beta1": round(beta1, 4),
        },
    }
