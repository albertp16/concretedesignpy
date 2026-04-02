# MIT License
# Copyright (c) Albert Pamonag Engineering Consultancy

"""
Beam Deflection Calculator
============================

Computes short-term and long-term deflection of RC beams
per NSCP 2015 Section 424.2 / ACI 318-19.
"""

import math


def _max_deflection_limit(clearspan, member_type):
    """
    Maximum allowable deflection per NSCP 2015 Table 424.2.2.

    Parameters
    ----------
    clearspan : float
        Clear span (mm).
    member_type : str
        One of: 'flat_roof', 'floor', 'roof_with_partitions',
        'floor_with_partitions'.

    Returns
    -------
    float
        Allowable deflection (mm).
    """
    limits = {
        "flat_roof": 180,
        "floor": 360,
        "roof_with_partitions": 480,
        "floor_with_partitions": 240,
    }
    if member_type not in limits:
        raise ValueError(f"Unknown member_type: {member_type}")
    return clearspan / limits[member_type]


def _time_dependent_factor(sustained_duration):
    """
    Time-dependent factor per NSCP 2015 Table 424.2.4.1.3.

    Parameters
    ----------
    sustained_duration : str
        One of: '3_months', '6_months', '12_months', '5_years_or_more'.

    Returns
    -------
    float
    """
    factors = {
        "3_months": 1.0,
        "6_months": 1.2,
        "12_months": 1.4,
        "5_years_or_more": 2.0,
    }
    if sustained_duration not in factors:
        raise ValueError(f"Unknown duration: {sustained_duration}")
    return factors[sustained_duration]


def deflection_computation(
    b, h, d, fc, fy, clearspan, as_tension, n_bars_comp,
    db_comp, es=200000.0, ma=None, point_load=0, uniform_load=0,
    beam_type="simply_supported", member_type="floor",
    sustained_duration="5_years_or_more",
):
    """
    Compute beam deflection (short-term and long-term).

    Parameters
    ----------
    b : float
        Beam width (mm).
    h : float
        Total height (mm).
    d : float
        Effective depth (mm).
    fc : float
        Concrete compressive strength (MPa).
    fy : float
        Steel yield strength (MPa).
    clearspan : float
        Clear span (mm).
    as_tension : float
        Area of tension reinforcement (mm^2).
    n_bars_comp : int
        Number of compression bars.
    db_comp : float
        Diameter of compression bars (mm).
    es : float
        Steel modulus of elasticity (MPa).
    ma : float or None
        Maximum applied moment (N-mm). Computed from loads if None.
    point_load : float
        Concentrated load at midspan (N).
    uniform_load : float
        Uniform load (N/mm).
    beam_type : str
        One of: 'simply_supported', 'cantilever', 'fixed_both',
        'propped_cantilever'.
    member_type : str
        Deflection limit category.
    sustained_duration : str
        Duration for long-term factor.

    Returns
    -------
    dict
        Keys: ig, fr, mcr, n, icr, ie, short_deflection, long_deflection,
        total_deflection, max_allowable, status
    """
    if b <= 0 or h <= 0 or d <= 0 or fc <= 0:
        raise ValueError("b, h, d, fc must all be positive.")

    # Gross moment of inertia
    ig = b * h ** 3 / 12.0

    # Modulus of rupture
    fr = 0.62 * math.sqrt(fc)

    # Cracking moment
    yt = h / 2.0
    mcr = fr * ig / yt

    # Concrete modulus
    ec = 4700.0 * math.sqrt(fc)

    # Modular ratio
    n = max(6, math.ceil(es / ec))

    # Cracked moment of inertia (simplified for singly reinforced)
    # Using transformed section: b*kd^3/3 + n*As*(d-kd)^2
    # Solve for kd from: b*kd^2/2 = n*As*(d - kd)
    a_coeff = b / 2.0
    b_coeff = n * as_tension
    c_coeff = -n * as_tension * d
    discriminant = b_coeff ** 2 - 4 * a_coeff * c_coeff
    kd = (-b_coeff + math.sqrt(discriminant)) / (2 * a_coeff)
    icr = b * kd ** 3 / 3.0 + n * as_tension * (d - kd) ** 2

    # Maximum moment from loads if not provided
    if ma is None:
        l = clearspan
        if beam_type == "simply_supported":
            ma = point_load * l / 4.0 + uniform_load * l ** 2 / 8.0
        elif beam_type == "cantilever":
            ma = point_load * l + uniform_load * l ** 2 / 2.0
        elif beam_type == "fixed_both":
            ma = point_load * l / 8.0 + uniform_load * l ** 2 / 12.0
        elif beam_type == "propped_cantilever":
            ma = point_load * l / 4.0 + uniform_load * l ** 2 / 8.0
        else:
            raise ValueError(f"Unknown beam_type: {beam_type}")

    # Effective moment of inertia (Branson's equation)
    if ma <= 0:
        ie = ig
    else:
        ratio = (mcr / ma) ** 3
        ie = ratio * ig + (1 - ratio) * icr
        ie = min(ie, ig)

    # Short-term deflection
    l = clearspan
    if beam_type == "simply_supported":
        delta_p = point_load * l ** 3 / (48.0 * ec * ie) if point_load else 0
        delta_w = 5.0 * uniform_load * l ** 4 / (384.0 * ec * ie) if uniform_load else 0
    elif beam_type == "cantilever":
        delta_p = point_load * l ** 3 / (3.0 * ec * ie) if point_load else 0
        delta_w = uniform_load * l ** 4 / (8.0 * ec * ie) if uniform_load else 0
    elif beam_type == "fixed_both":
        delta_p = point_load * l ** 3 / (192.0 * ec * ie) if point_load else 0
        delta_w = uniform_load * l ** 4 / (384.0 * ec * ie) if uniform_load else 0
    else:
        delta_p = point_load * l ** 3 / (48.0 * ec * ie) if point_load else 0
        delta_w = 5.0 * uniform_load * l ** 4 / (384.0 * ec * ie) if uniform_load else 0

    short_defl = delta_p + delta_w

    # Long-term deflection
    time_factor = _time_dependent_factor(sustained_duration)
    as_comp = n_bars_comp * math.pi * (db_comp / 2.0) ** 2
    rho_prime = as_comp / (b * d) if (b * d) > 0 else 0
    lambda_delta = time_factor / (1 + 50 * rho_prime)
    long_defl = lambda_delta * short_defl

    total_defl = short_defl + long_defl
    max_allow = _max_deflection_limit(clearspan, member_type)

    return {
        "ig": round(ig, 2),
        "fr": round(fr, 4),
        "mcr": round(mcr, 2),
        "mcr_knm": round(mcr / 1e6, 2),
        "n": n,
        "kd": round(kd, 2),
        "icr": round(icr, 2),
        "ie": round(ie, 2),
        "ec": round(ec, 2),
        "short_deflection": round(short_defl, 4),
        "long_deflection": round(long_defl, 4),
        "total_deflection": round(total_defl, 4),
        "max_allowable": round(max_allow, 4),
        "lambda_delta": round(lambda_delta, 4),
        "time_factor": time_factor,
        "status": "OK" if total_defl <= max_allow else "EXCEEDS LIMIT",
    }
