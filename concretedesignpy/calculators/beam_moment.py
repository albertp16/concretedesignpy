# MIT License
# Copyright (c) Albert Pamonag Engineering Consultancy

"""
Beam Flexural Capacity Calculator
==================================

Computes nominal and ultimate moment capacity of RC beams
using iterative neutral axis solver with strain compatibility.

Reference: NSCP 2015 / ACI 318-19
"""

import math


def _calculate_beta_one(fc):
    """Whitney stress block factor per NSCP 2015 Section 422.2.2.4.3."""
    if fc <= 17:
        return 0.85
    elif fc <= 28:
        return 0.85
    elif fc < 55:
        return max(0.85 - 0.05 * (fc - 28) / 7, 0.65)
    else:
        return 0.65


def _strength_reduction_factor(epsilon_t, epsilon_ty):
    """Strength reduction factor (phi) based on tensile strain."""
    if epsilon_t >= 0.005:
        return 0.90, "tension-controlled"
    elif epsilon_t <= epsilon_ty:
        return 0.65, "compression-controlled"
    else:
        phi = 0.65 + 0.25 * (epsilon_t - epsilon_ty) / (0.005 - epsilon_ty)
        return phi, "transition"


def calculate_beam_moment(rebar_list, fc, fy, b, h, es=200000.0):
    """
    Calculate moment capacity of an RC beam section.

    Parameters
    ----------
    rebar_list : list of dict
        Each dict has keys:
            - 'd': depth from top (mm)
            - 'diam': bar diameter (mm)
            - 'num': number of bars
    fc : float
        Concrete compressive strength (MPa).
    fy : float
        Steel yield strength (MPa).
    b : float
        Beam width (mm).
    h : float
        Beam total height (mm).
    es : float, optional
        Steel modulus of elasticity (MPa). Default 200000.

    Returns
    -------
    dict
        Keys: neutral_axis, mn, mu, phi, classification, rebar_forces
    """
    if not rebar_list:
        raise ValueError("rebar_list must not be empty.")
    if fc <= 0 or fy <= 0 or b <= 0 or h <= 0:
        raise ValueError("fc, fy, b, h must all be positive.")

    beta1 = _calculate_beta_one(fc)
    epsilon_ty = fy / es
    ecu = 0.003

    # Process rebar data
    rebars = []
    for bar in rebar_list:
        area = (math.pi / 4) * bar["diam"] ** 2 * bar["num"]
        rebars.append({"d": bar["d"], "area": area})

    # Phase 1: coarse search for neutral axis
    step_coarse = 1.0
    c = h
    prev_ratio = None

    while c > 0:
        a = beta1 * c
        fc_concrete = 0.85 * fc * a * b

        fs_total = 0.0
        for rb in rebars:
            strain = ecu * (rb["d"] - c) / c
            stress = max(-fy, min(fy, strain * es))
            fs_total += stress * rb["area"]

        ratio = fc_concrete / abs(fs_total) if abs(fs_total) > 1e-9 else 999.0

        if prev_ratio is not None and ratio <= 1.0 and prev_ratio > 1.0:
            break
        prev_ratio = ratio
        c -= step_coarse

    # Phase 2: fine refinement
    c_start = c + step_coarse
    step_fine = 0.01
    c = c_start

    for _ in range(int(step_coarse / step_fine) + 200):
        a = beta1 * c
        fc_concrete = 0.85 * fc * a * b

        fs_total = 0.0
        rebar_forces = []
        for rb in rebars:
            strain = ecu * (rb["d"] - c) / c
            stress = max(-fy, min(fy, strain * es))
            force = stress * rb["area"]
            fs_total += force
            rebar_forces.append({
                "d": rb["d"],
                "area": rb["area"],
                "strain": strain,
                "stress": stress,
                "force": force,
            })

        if abs(fc_concrete - fs_total) < 0.5:
            break
        if fc_concrete < fs_total:
            break
        c -= step_fine

    # Compute moment about top fiber
    a = beta1 * c
    mn = 0.0
    for rf in rebar_forces:
        mn += rf["force"] * (rf["d"] - a / 2)
    mn_knm = mn / 1e6

    # Tensile strain at extreme tension bar
    d_max = max(rb["d"] for rb in rebars)
    epsilon_t = ecu * (d_max - c) / c
    phi, classification = _strength_reduction_factor(epsilon_t, epsilon_ty)
    mu_knm = phi * mn_knm

    return {
        "neutral_axis": round(c, 2),
        "a": round(a, 2),
        "beta1": round(beta1, 4),
        "mn": round(mn_knm, 2),
        "mu": round(mu_knm, 2),
        "phi": round(phi, 4),
        "classification": classification,
        "epsilon_t": round(epsilon_t, 6),
        "fc_concrete": round(fc_concrete / 1000, 2),
        "fs_total": round(fs_total / 1000, 2),
        "rebar_forces": rebar_forces,
    }
