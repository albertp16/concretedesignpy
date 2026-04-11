# MIT License
# Copyright (c) Albert Pamonag Engineering Consultancy

"""
Moment-Curvature Analysis
===========================

Computes the moment-curvature relationship for rectangular RC sections.

Modes:
    1. Quick (6-Point): Closed-form key points using Whitney block
    2. Advanced (Incremental): Fiber approach with selectable concrete model

Concrete Models:
    - Hognestad (1951): Parabolic ascending + linear descending (unconfined)
    - Mander et al. (1988): Power-law curve for confined concrete

Features:
    - Compression steel support (doubly-reinforced sections)
    - Ductility ratio computation (mu = phi_u / phi_y)
    - Event detection (cracking, yield, peak, ultimate)

Reference:
    NSCP 2015 / ACI 318-19
    Anwar & Najam (2017) Ch. 6, Eqs. 6.1-6.7
    Hognestad (1951) — Univ. of Illinois Bulletin No. 399
    Mander, Priestley & Park (1988) — J. Structural Eng., ASCE
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


# ---------------------------------------------------------------------------
# Advanced incremental M-phi with Hognestad concrete & axial load
# ---------------------------------------------------------------------------

def _hognestad_stress(ec_strain, fc, eco, ecu):
    """
    Hognestad concrete stress-strain model (compression positive).

    Ascending: parabolic up to eco
    Descending: linear from eco to ecu (stress drops to 0.85*fc)
    """
    if ec_strain <= 0:
        return 0.0
    if ec_strain <= eco:
        ratio = ec_strain / eco
        return fc * (2 * ratio - ratio ** 2)
    if ec_strain <= ecu:
        return fc * (1 - 0.15 * (ec_strain - eco) / (ecu - eco))
    return 0.0


def _mander_stress(ec_strain, fcc, ecc, ecu, r):
    """
    Mander confined concrete stress-strain model.

    Reference: Mander, Priestley & Park (1988), Eq. 6.15a
    """
    if ec_strain <= 0 or ec_strain > ecu:
        return 0.0
    x = ec_strain / ecc if ecc > 0 else 0
    denom = r - 1 + x ** r
    if denom == 0:
        return 0.0
    return fcc * x * r / denom


def _steel_stress(strain, fy, es):
    """Elastic-perfectly-plastic steel model."""
    fs = es * strain
    if fs > fy:
        return fy
    if fs < -fy:
        return -fy
    return fs


def moment_curvature_advanced(
    b, h, d, fc, fy, as_tension, es=200000.0,
    axial_load=0.0, n_increments=60, n_fibers=50,
    d_prime=0.0, as_compression=0.0,
    concrete_model="hognestad", mander_params=None,
):
    """
    Compute moment-curvature using incremental fiber approach.

    Supports Hognestad (unconfined) or Mander (confined) concrete,
    compression steel, ductility computation, and event detection.

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
    axial_load : float
        Axial load in kN (positive = compression).
    n_increments : int
        Number of curvature increments.
    n_fibers : int
        Number of concrete fiber strips.
    d_prime : float
        Distance from compression face to compression steel (mm).
    as_compression : float
        Area of compression reinforcement (mm^2).
    concrete_model : str
        'hognestad' or 'mander'.
    mander_params : dict or None
        Result from confined_stress_strain() when using Mander model.
        Required keys: fcc, ecc, ecu, r

    Returns
    -------
    dict
        Keys: points, section_properties, hognestad_params or
        mander_curve_params, ductility, events
    """
    if b <= 0 or h <= 0 or d <= 0 or fc <= 0 or fy <= 0:
        raise ValueError("All section/material parameters must be positive.")

    ec_mod = 4700 * math.sqrt(fc)  # concrete elastic modulus
    fr = 0.62 * math.sqrt(fc)
    ey = fy / es

    # Select concrete model parameters
    use_mander = (concrete_model == "mander" and mander_params is not None)
    if use_mander:
        fcc = mander_params['fcc']
        ecc = mander_params['ecc']
        ecu = mander_params['ecu']
        r_mander = mander_params['r']
        eco = ecc  # peak strain for Mander
    else:
        eco = 2 * fc / ec_mod           # strain at peak stress
        ecu = 0.003                     # ultimate concrete strain

    # Axial load in N (input is kN)
    p_axial = axial_load * 1000.0

    # Fiber geometry: strips from bottom (y=0) to top (y=h)
    fiber_h = h / n_fibers
    fiber_y = [(i + 0.5) * fiber_h for i in range(n_fibers)]

    # Steel layers
    steel_y_tens = h - d  # tension steel y from bottom
    steel_y_comp = h - d_prime if (as_compression > 0 and d_prime > 0) else None

    centroid = h / 2.0

    points = []

    # Event detection state
    events = []
    cracking_detected = False
    yield_detected = False
    peak_moment = 0.0
    peak_moment_point = None

    for step in range(1, n_increments + 1):
        ec_top = step * ecu / n_increments

        # Binary search for neutral axis depth c (from top)
        c_lo, c_hi = 1.0, h * 2.0
        c_found = None
        final_total_moment = 0.0
        final_strain_steel = 0.0

        for _ in range(100):
            c = (c_lo + c_hi) / 2.0

            # Concrete forces (fiber approach)
            total_force = 0.0
            total_moment = 0.0
            for fy_i in fiber_y:
                dist_from_top = h - fy_i
                strain_i = ec_top * (c - dist_from_top) / c if c > 0 else 0
                if use_mander:
                    stress_i = _mander_stress(strain_i, fcc, ecc, ecu, r_mander)
                else:
                    stress_i = _hognestad_stress(strain_i, fc, eco, ecu)
                force_i = stress_i * b * fiber_h
                total_force += force_i
                total_moment += force_i * (fy_i - centroid)

            # Tension steel
            strain_steel = ec_top * (c - d) / c if c > 0 else 0
            stress_steel = _steel_stress(-strain_steel, fy, es)
            force_steel = -stress_steel * as_tension
            total_force += force_steel
            total_moment += force_steel * (steel_y_tens - centroid)

            # Compression steel
            if steel_y_comp is not None:
                strain_comp = ec_top * (c - d_prime) / c if c > 0 else 0
                stress_comp = _steel_stress(strain_comp, fy, es)
                force_comp = stress_comp * as_compression
                total_force += force_comp
                total_moment += force_comp * (steel_y_comp - centroid)

            # Axial load
            net_force = total_force + p_axial

            if abs(net_force) < 0.5:
                c_found = c
                final_total_moment = total_moment
                final_strain_steel = strain_steel
                break

            if net_force > 0:
                c_hi = c
            else:
                c_lo = c

        if c_found is None:
            continue

        phi = ec_top / c_found if c_found > 0 else 0

        points.append({
            "ec_top": round(ec_top, 8),
            "phi": round(phi, 10),
            "moment_knm": round(final_total_moment / 1e6, 2),
            "c": round(c_found, 2),
        })

        # ── Event detection ──
        moment_val = abs(final_total_moment / 1e6)

        # Cracking: tension fiber strain exceeds fr/Ec
        if not cracking_detected:
            max_tension_strain = ec_top * (c_found - h) / c_found if c_found > 0 else 0
            if abs(max_tension_strain) * ec_mod >= fr:
                cracking_detected = True
                events.append({
                    "event": "Cracking",
                    "phi": round(phi, 10),
                    "moment_knm": round(final_total_moment / 1e6, 2),
                    "step": step,
                })

        # Steel yield: tension steel strain >= ey
        if not yield_detected:
            if abs(final_strain_steel) >= ey:
                yield_detected = True
                events.append({
                    "event": "Steel Yield",
                    "phi": round(phi, 10),
                    "moment_knm": round(final_total_moment / 1e6, 2),
                    "step": step,
                })

        # Track peak moment
        if moment_val > peak_moment:
            peak_moment = moment_val
            peak_moment_point = {
                "event": "Peak Moment",
                "phi": round(phi, 10),
                "moment_knm": round(final_total_moment / 1e6, 2),
                "step": step,
            }

    # Add peak moment event
    if peak_moment_point:
        events.append(peak_moment_point)

    # Add ultimate event (last converged point)
    if points:
        last = points[-1]
        events.append({
            "event": "Ultimate",
            "phi": last["phi"],
            "moment_knm": last["moment_knm"],
            "step": len(points),
        })

    # ── Ductility computation ──
    ductility = None
    yield_evt = next((e for e in events if e["event"] == "Steel Yield"), None)
    ult_evt = next((e for e in events if e["event"] == "Ultimate"), None)
    if yield_evt and ult_evt and yield_evt["phi"] > 0:
        phi_y = yield_evt["phi"]
        phi_u = ult_evt["phi"]
        mu = phi_u / phi_y
        ei_yield = (yield_evt["moment_knm"] * 1e6) / phi_y if phi_y > 0 else 0
        ei_ultimate = (ult_evt["moment_knm"] * 1e6) / phi_u if phi_u > 0 else 0
        ductility = {
            "phi_yield": phi_y,
            "moment_yield_knm": yield_evt["moment_knm"],
            "phi_ultimate": phi_u,
            "moment_ultimate_knm": ult_evt["moment_knm"],
            "mu": round(mu, 2),
            "ei_yield": round(ei_yield, 0),
            "ei_ultimate": round(ei_ultimate, 0),
        }

    # ── Stress-strain curve data for plotting ──
    if use_mander:
        n_curve = 200
        curve_strains = [i * ecu / n_curve for i in range(n_curve + 1)]
        curve_stresses = [round(_mander_stress(e, fcc, ecc, ecu, r_mander), 4)
                          for e in curve_strains]
        curve_strains = [round(e, 8) for e in curve_strains]
        model_params = {
            "model": "mander",
            "fcc": mander_params['fcc'],
            "ecc": mander_params['ecc'],
            "ecu": mander_params['ecu'],
            "r": mander_params['r'],
            "fc": fc,
            "eco_unconfined": round(2 * fc / ec_mod, 6),
            "strains": curve_strains,
            "stresses": curve_stresses,
        }
    else:
        hog_strains = [i * ecu / 100 for i in range(101)]
        hog_stresses = [round(_hognestad_stress(e, fc, eco, ecu), 4)
                        for e in hog_strains]
        hog_strains = [round(e, 6) for e in hog_strains]
        model_params = {
            "model": "hognestad",
            "eco": round(eco, 6),
            "ecu": ecu,
            "fc": fc,
            "strains": hog_strains,
            "stresses": hog_stresses,
        }

    return {
        "points": points,
        "section_properties": {
            "ec_mod": round(ec_mod, 2),
            "n": round(es / ec_mod, 4),
            "fr": round(fr, 4),
        },
        "concrete_model_params": model_params,
        "events": events,
        "ductility": ductility,
        "compression_steel": {
            "d_prime": d_prime,
            "as_compression": as_compression,
        } if as_compression > 0 else None,
        # Backward compat: keep hognestad_params for quick+adv overlay
        "hognestad_params": model_params if not use_mander else {
            "eco": round(2 * fc / ec_mod, 6),
            "ecu": 0.003,
            "fc": fc,
            "strains": [round(i * 0.003 / 100, 6) for i in range(101)],
            "stresses": [round(_hognestad_stress(i * 0.003 / 100, fc,
                         2 * fc / ec_mod, 0.003), 4) for i in range(101)],
        },
    }
