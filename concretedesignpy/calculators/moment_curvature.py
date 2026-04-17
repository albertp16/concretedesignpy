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


# ---------------------------------------------------------------------------
# Polygon geometry helpers (for Special / custom-polygon sections)
# ---------------------------------------------------------------------------

def polygon_area(vertices):
    """Signed-area via shoelace; returns absolute area."""
    n = len(vertices)
    if n < 3:
        return 0.0
    a = 0.0
    for i in range(n):
        x1, y1 = vertices[i]
        x2, y2 = vertices[(i + 1) % n]
        a += x1 * y2 - x2 * y1
    return abs(a) / 2.0


def polygon_centroid(vertices):
    """Centroid (x_bar, y_bar) of a simple polygon."""
    n = len(vertices)
    a = 0.0
    cx = 0.0
    cy = 0.0
    for i in range(n):
        x1, y1 = vertices[i]
        x2, y2 = vertices[(i + 1) % n]
        cross = x1 * y2 - x2 * y1
        a += cross
        cx += (x1 + x2) * cross
        cy += (y1 + y2) * cross
    a /= 2.0
    if abs(a) < 1e-12:
        return 0.0, 0.0
    return cx / (6.0 * a), cy / (6.0 * a)


def polygon_inertia_x(vertices, y_axis=None):
    """
    Second moment of area about a horizontal axis at y=y_axis.
    If y_axis is None, uses the centroidal axis.
    """
    n = len(vertices)
    # Integral of y^2 dA via Green's theorem
    ix_origin = 0.0
    for i in range(n):
        x1, y1 = vertices[i]
        x2, y2 = vertices[(i + 1) % n]
        cross = x1 * y2 - x2 * y1
        ix_origin += cross * (y1 * y1 + y1 * y2 + y2 * y2)
    ix_origin = abs(ix_origin) / 12.0
    area = polygon_area(vertices)
    _, y_bar = polygon_centroid(vertices)
    ix_centroid = ix_origin - area * y_bar * y_bar
    if y_axis is None:
        return ix_centroid
    return ix_centroid + area * (y_bar - y_axis) ** 2


def polygon_width_at(vertices, y):
    """
    Total horizontal extent of the polygon at height y
    (sum of interval widths for non-convex shapes).
    """
    intersections = []
    n = len(vertices)
    eps = 1e-9
    for i in range(n):
        x1, y1 = vertices[i]
        x2, y2 = vertices[(i + 1) % n]
        if abs(y2 - y1) < eps:
            continue
        y_lo, y_hi = (y1, y2) if y1 < y2 else (y2, y1)
        if y < y_lo - eps or y > y_hi + eps:
            continue
        t = (y - y1) / (y2 - y1)
        intersections.append(x1 + t * (x2 - x1))
    intersections.sort()
    width = 0.0
    for k in range(0, len(intersections) - 1, 2):
        width += intersections[k + 1] - intersections[k]
    return width


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
    section_vertices=None,
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

    # Fiber geometry: strips from bottom (y=0) to top (y=h).
    # For polygon sections, width varies with y via polygon_width_at().
    use_polygon = section_vertices is not None and len(section_vertices) >= 3
    fiber_h = h / n_fibers
    fiber_y = [(i + 0.5) * fiber_h for i in range(n_fibers)]
    if use_polygon:
        fiber_b = [polygon_width_at(section_vertices, y) for y in fiber_y]
    else:
        fiber_b = [b] * n_fibers

    # Steel layers
    steel_y_tens = h - d  # tension steel y from bottom
    steel_y_comp = h - d_prime if (as_compression > 0 and d_prime > 0) else None

    if use_polygon:
        _, centroid = polygon_centroid(section_vertices)
    else:
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
            for idx_f, fy_i in enumerate(fiber_y):
                dist_from_top = h - fy_i
                strain_i = ec_top * (c - dist_from_top) / c if c > 0 else 0
                if use_mander:
                    stress_i = _mander_stress(strain_i, fcc, ecc, ecu, r_mander)
                else:
                    stress_i = _hognestad_stress(strain_i, fc, eco, ecu)
                force_i = stress_i * fiber_b[idx_f] * fiber_h
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

    # ── Fiber data at ultimate for plotting ──
    fiber_plot = None
    if points:
        last_pt = points[-1]
        ec_ult = last_pt["ec_top"]
        c_ult = last_pt["c"]
        fiber_strains = []
        fiber_stresses = []
        for fy_i in fiber_y:
            dist_from_top = h - fy_i
            strain_i = ec_ult * (c_ult - dist_from_top) / c_ult if c_ult > 0 else 0
            if use_mander:
                stress_i = _mander_stress(strain_i, fcc, ecc, ecu, r_mander)
            else:
                stress_i = _hognestad_stress(strain_i, fc, eco, ecu)
            fiber_strains.append(round(strain_i, 8))
            fiber_stresses.append(round(stress_i, 4))
        # Rebar strains/stresses at ultimate
        rebar_data = []
        # Tension steel
        strain_t = ec_ult * (c_ult - d) / c_ult if c_ult > 0 else 0
        stress_t = _steel_stress(-strain_t, fy, es)
        rebar_data.append({
            "y": round(h - d, 1), "strain": round(-strain_t, 8),
            "stress": round(stress_t, 2), "label": "Tension",
        })
        # Compression steel
        if as_compression > 0 and d_prime > 0:
            strain_c = ec_ult * (c_ult - d_prime) / c_ult if c_ult > 0 else 0
            stress_c = _steel_stress(strain_c, fy, es)
            rebar_data.append({
                "y": round(h - d_prime, 1), "strain": round(strain_c, 8),
                "stress": round(stress_c, 2), "label": "Compression",
            })

        fiber_plot = {
            "y": [round(fy_i, 1) for fy_i in fiber_y],
            "widths": [round(w, 2) for w in fiber_b],
            "strains": fiber_strains,
            "stresses": fiber_stresses,
            "c": c_ult,
            "ec_top": ec_ult,
            "n_fibers": n_fibers,
            "fiber_h": round(fiber_h, 2),
            "rebars": rebar_data,
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

    section_props = {
        "ec_mod": round(ec_mod, 2),
        "n": round(es / ec_mod, 4),
        "fr": round(fr, 4),
    }
    if use_polygon:
        area_poly = polygon_area(section_vertices)
        ig_poly = polygon_inertia_x(section_vertices)
        section_props.update({
            "area": round(area_poly, 2),
            "ig": round(ig_poly, 2),
            "y_bar": round(centroid, 2),
            "section_type": "polygon",
            "vertices": [[round(x, 3), round(y, 3)] for x, y in section_vertices],
        })

    return {
        "points": points,
        "section_properties": section_props,
        "concrete_model_params": model_params,
        "events": events,
        "ductility": ductility,
        "fiber_plot": fiber_plot,
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


# ---------------------------------------------------------------------------
# Moment-Curvature via OpenSeesPy (zero-length section element)
# ---------------------------------------------------------------------------

def moment_curvature_opensees(
    b, h, d, fc, fy, as_tension, es=200000.0,
    axial_load=0.0, n_increments=100,
    d_prime=0.0, as_compression=0.0,
    cover=40.0, n_bars_tension=4, n_bars_comp=0,
    concrete_model="hognestad", mander_params=None,
):
    """
    Moment-curvature analysis using OpenSeesPy.

    Uses a zero-length section element with fiber discretisation,
    Concrete01 (Kent-Park) for concrete and Steel01 for reinforcement.

    Parameters
    ----------
    b, h, d, fc, fy, as_tension, es :
        Same as moment_curvature_advanced().
    axial_load : float
        Axial load in kN (positive = compression).
    n_increments : int
        Number of curvature steps.
    d_prime : float
        Distance from compression face to compression steel (mm).
    as_compression : float
        Compression steel area (mm^2).
    cover : float
        Clear cover to ties/stirrups (mm).
    n_bars_tension, n_bars_comp : int
        Number of tension/compression bars.
    concrete_model : str
        'hognestad' or 'mander'.
    mander_params : dict or None
        Result from confined_stress_strain() with keys fcc, ecc, ecu, r.

    Returns
    -------
    dict
        Compatible with moment_curvature_advanced() output format.
    """
    try:
        import openseespy.opensees as ops
    except (ImportError, RuntimeError, OSError) as e:
        raise ImportError(
            "openseespy is not available: " + str(e)
            + " — Install a compatible version with: pip install openseespy"
        )

    if b <= 0 or h <= 0 or d <= 0 or fc <= 0 or fy <= 0:
        raise ValueError("All section/material parameters must be positive.")

    ec_mod = 4700 * math.sqrt(fc)
    fr = 0.62 * math.sqrt(fc)
    ey = fy / es
    eco = 2 * fc / ec_mod

    use_mander = concrete_model == "mander" and mander_params is not None

    # ── Build OpenSees model ──
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 3)

    # -- Material definitions (Concrete01: fpc, epsc0, fpcu, epsU) --
    # All values negative per OpenSees sign convention for compression.
    if use_mander:
        fcc = mander_params['fcc']
        ecc = mander_params['ecc']
        ecu_conf = mander_params['ecu']
        # Core – confined
        ops.uniaxialMaterial('Concrete01', 1,
                             -fcc, -ecc, -0.2 * fcc, -ecu_conf)
        # Cover – unconfined
        ops.uniaxialMaterial('Concrete01', 2,
                             -fc, -eco, 0.0, -0.006)
        ecu_analysis = ecu_conf
    else:
        ecu_hog = 0.003
        # Core (lightly confined)
        ops.uniaxialMaterial('Concrete01', 1,
                             -fc, -eco, -0.2 * fc, -0.005)
        # Cover (unconfined)
        ops.uniaxialMaterial('Concrete01', 2,
                             -fc, -eco, 0.0, -0.006)
        ecu_analysis = ecu_hog

    # Steel01: Fy, E0, b (strain-hardening ratio)
    ops.uniaxialMaterial('Steel01', 3, fy, es, 0.01)

    # -- Fiber section --
    y1 = h / 2.0   # half-depth
    z1 = b / 2.0   # half-width
    n_core = max(20, n_increments // 2)

    ops.section('Fiber', 1)

    # Core concrete (inside cover)
    ops.patch('rect', 1, n_core, 1,
              cover - y1, cover - z1, y1 - cover, z1 - cover)

    # Cover patches (top, bottom, left, right)
    ops.patch('rect', 2, 2, 1, y1 - cover, -z1, y1, z1)
    ops.patch('rect', 2, 2, 1, -y1, -z1, -(y1 - cover), z1)
    ops.patch('rect', 2, n_core, 1,
              -(y1 - cover), -z1, y1 - cover, -(z1 - cover))
    ops.patch('rect', 2, n_core, 1,
              -(y1 - cover), z1 - cover, y1 - cover, z1)

    # Tension reinforcement layer
    y_tens = y1 - d   # negative: below centroid
    if n_bars_tension > 0 and as_tension > 0:
        a_bar = as_tension / n_bars_tension
        if n_bars_tension == 1:
            ops.fiber(y_tens, 0.0, a_bar, 3)
        else:
            ops.layer('straight', 3, n_bars_tension, a_bar,
                      y_tens, -(z1 - cover), y_tens, z1 - cover)

    # Compression reinforcement layer
    if n_bars_comp > 0 and as_compression > 0 and d_prime > 0:
        y_comp = y1 - d_prime
        a_bar_c = as_compression / n_bars_comp
        if n_bars_comp == 1:
            ops.fiber(y_comp, 0.0, a_bar_c, 3)
        else:
            ops.layer('straight', 3, n_bars_comp, a_bar_c,
                      y_comp, -(z1 - cover), y_comp, z1 - cover)

    # -- Nodes & element --
    ops.node(1, 0.0, 0.0)
    ops.node(2, 0.0, 0.0)
    ops.fix(1, 1, 1, 1)
    ops.fix(2, 0, 1, 0)
    ops.element('zeroLengthSection', 1, 1, 2, 1)

    # -- Phase 1: apply axial load --
    p_axial = axial_load * 1000.0  # kN → N
    ops.timeSeries('Constant', 1)
    ops.pattern('Plain', 1, 1)
    ops.load(2, p_axial, 0.0, 0.0)

    ops.integrator('LoadControl', 0)
    ops.system('BandGeneral')
    ops.test('NormUnbalance', 1.0e-9, 10)
    ops.numberer('Plain')
    ops.constraints('Plain')
    ops.algorithm('Newton')
    ops.analysis('Static')
    ops.analyze(1)

    ops.loadConst('-time', 0.0)

    # -- Phase 2: displacement-controlled curvature --
    c_approx = d * 0.4
    max_k = ecu_analysis / c_approx
    dk = max_k / n_increments

    ops.timeSeries('Linear', 2)
    ops.pattern('Plain', 2, 2)
    ops.load(2, 0.0, 0.0, 1.0)
    ops.integrator('DisplacementControl', 2, 3, dk)
    ops.analysis('Static')

    phi_list = []
    moment_list = []

    for i in range(n_increments):
        ok = ops.analyze(1)
        if ok != 0:
            ops.algorithm('ModifiedNewton')
            ok = ops.analyze(1)
        if ok != 0:
            ops.algorithm('Newton', '-initial')
            ok = ops.analyze(1)
        if ok != 0:
            break  # use data collected so far
        ops.algorithm('Newton')

        phi = ops.nodeDisp(2, 3)
        forces = ops.eleForce(1)
        moment = -forces[2]  # N·mm, sign convention

        phi_list.append(phi)
        moment_list.append(moment / 1e6)  # → kN·m

    ops.wipe()

    if not phi_list:
        raise RuntimeError("OpenSeesPy analysis failed to converge.")

    # ── Build points list ──
    points = []
    for idx in range(len(phi_list)):
        points.append({
            "ec_top": 0.0,
            "phi": round(phi_list[idx], 10),
            "moment_knm": round(moment_list[idx], 4),
            "c": 0.0,
        })

    # ── Event detection ──
    events = []

    # Cracking moment
    ig = b * h ** 3 / 12
    mcr = fr * ig / (h / 2) / 1e6  # kN·m
    for idx, m in enumerate(moment_list):
        if abs(m) >= mcr:
            events.append({
                "event": "Cracking",
                "phi": round(phi_list[idx], 10),
                "moment_knm": round(m, 4),
                "step": idx + 1,
            })
            break

    # Steel yield: approximate phi_y = ey / (d - c_approx)
    # Detect by slope change: when tangent stiffness drops below
    # 20% of the initial stiffness
    if len(phi_list) > 2:
        ei_init = abs(moment_list[1] - moment_list[0]) / \
                  abs(phi_list[1] - phi_list[0]) if abs(phi_list[1] - phi_list[0]) > 0 else 0
        for idx in range(2, len(phi_list)):
            dphi = phi_list[idx] - phi_list[idx - 1]
            dm = moment_list[idx] - moment_list[idx - 1]
            if dphi <= 0:
                continue
            ei_local = abs(dm / dphi)
            if ei_init > 0 and ei_local < 0.20 * ei_init:
                events.append({
                    "event": "Steel Yield",
                    "phi": round(phi_list[idx], 10),
                    "moment_knm": round(moment_list[idx], 4),
                    "step": idx + 1,
                })
                break

    # Peak moment
    peak_idx = max(range(len(moment_list)), key=lambda i: abs(moment_list[i]))
    events.append({
        "event": "Peak Moment",
        "phi": round(phi_list[peak_idx], 10),
        "moment_knm": round(moment_list[peak_idx], 4),
        "step": peak_idx + 1,
    })

    # Ultimate (last converged point)
    events.append({
        "event": "Ultimate",
        "phi": round(phi_list[-1], 10),
        "moment_knm": round(moment_list[-1], 4),
        "step": len(phi_list),
    })

    # ── Ductility ──
    ductility = None
    yield_evt = next((e for e in events if e["event"] == "Steel Yield"), None)
    ult_evt = events[-1]  # Ultimate
    if yield_evt and yield_evt["phi"] > 0:
        phi_y = yield_evt["phi"]
        phi_u = ult_evt["phi"]
        mu = phi_u / phi_y
        ei_y = (yield_evt["moment_knm"] * 1e6) / phi_y
        ei_u = (ult_evt["moment_knm"] * 1e6) / phi_u if phi_u > 0 else 0
        ductility = {
            "phi_yield": phi_y,
            "moment_yield_knm": yield_evt["moment_knm"],
            "phi_ultimate": phi_u,
            "moment_ultimate_knm": ult_evt["moment_knm"],
            "mu": round(mu, 2),
            "ei_yield": round(ei_y, 0),
            "ei_ultimate": round(ei_u, 0),
        }

    # ── Concrete model curve data (for stress-strain chart) ──
    if use_mander:
        n_curve = 200
        curve_strains = [i * ecu_conf / n_curve for i in range(n_curve + 1)]
        curve_stresses = [
            round(_mander_stress(e, fcc, ecc, ecu_conf,
                                 mander_params['r']), 4)
            for e in curve_strains
        ]
        model_params = {
            "model": "mander",
            "fcc": fcc, "ecc": ecc, "ecu": ecu_conf,
            "r": mander_params['r'], "fc": fc,
            "eco_unconfined": round(eco, 6),
            "strains": [round(e, 8) for e in curve_strains],
            "stresses": curve_stresses,
        }
    else:
        hog_strains = [i * 0.003 / 100 for i in range(101)]
        hog_stresses = [
            round(_hognestad_stress(e, fc, eco, 0.003), 4)
            for e in hog_strains
        ]
        model_params = {
            "model": "hognestad",
            "eco": round(eco, 6), "ecu": 0.003, "fc": fc,
            "strains": [round(e, 6) for e in hog_strains],
            "stresses": hog_stresses,
        }

    return {
        "points": points,
        "events": events,
        "ductility": ductility,
        "engine": "opensees",
        "concrete_model_params": model_params,
        "section_properties": {
            "ec_mod": round(ec_mod, 2),
            "n": round(es / ec_mod, 4),
            "fr": round(fr, 4),
        },
        "fiber_plot": None,
        "hognestad_params": model_params if not use_mander else {
            "eco": round(eco, 6), "ecu": 0.003, "fc": fc,
            "strains": [round(i * 0.003 / 100, 6) for i in range(101)],
            "stresses": [round(_hognestad_stress(i * 0.003 / 100, fc,
                         eco, 0.003), 4) for i in range(101)],
        },
    }


# ---------------------------------------------------------------------------
# Moment-Rotation (M-theta) from M-phi
# ---------------------------------------------------------------------------

def plastic_hinge_length(method, h, fy=0, db=0, shear_span=0):
    """
    Compute plastic hinge length Lp using common empirical formulas.

    Parameters
    ----------
    method : str
        'manual' (returns 0), 'priestley', or 'park_paulay'.
    h : float
        Section depth (mm).
    fy : float
        Steel yield strength (MPa).
    db : float
        Tension bar diameter (mm).
    shear_span : float
        Shear span / member length (mm). Required for Priestley.

    Returns
    -------
    float
        Plastic hinge length Lp (mm).
    """
    if method == "priestley":
        # Priestley (1992): Lp = 0.08*L + 0.022*fy*db
        if shear_span <= 0 or fy <= 0 or db <= 0:
            raise ValueError("Priestley method requires shear_span, fy, db > 0.")
        return 0.08 * shear_span + 0.022 * fy * db
    elif method == "park_paulay":
        # Park & Paulay (1975): Lp = 0.5*h
        return 0.5 * h
    else:
        return 0.0


def moment_rotation_from_mphi(mphi_result, lp):
    """
    Convert moment-curvature results to moment-rotation.

    Applies theta = phi * Lp for each point and event.

    Parameters
    ----------
    mphi_result : dict
        Result from moment_curvature_advanced().
    lp : float
        Plastic hinge length (mm).

    Returns
    -------
    dict
        points, events, ductility_rotation, lp
    """
    if lp <= 0:
        raise ValueError("Plastic hinge length Lp must be positive.")

    # Convert M-phi points to M-theta
    points = []
    for p in mphi_result.get("points", []):
        phi = p.get("phi", 0)
        theta = phi * lp
        points.append({
            "theta_rad": round(theta, 10),
            "theta_mrad": round(theta * 1000, 6),
            "moment_knm": p["moment_knm"],
            "phi": p["phi"],
            "c": p.get("c", 0),
        })

    # Convert events
    events = []
    for e in mphi_result.get("events", []):
        theta = e["phi"] * lp
        events.append({
            "event": e["event"],
            "theta_rad": round(theta, 10),
            "theta_mrad": round(theta * 1000, 6),
            "moment_knm": e["moment_knm"],
            "phi": e["phi"],
        })

    # Rotation ductility
    ductility_rotation = None
    mphi_duct = mphi_result.get("ductility")
    if mphi_duct:
        theta_y = mphi_duct["phi_yield"] * lp
        theta_u = mphi_duct["phi_ultimate"] * lp
        mu_theta = theta_u / theta_y if theta_y > 0 else 0
        ductility_rotation = {
            "theta_y_rad": round(theta_y, 10),
            "theta_y_mrad": round(theta_y * 1000, 6),
            "theta_u_rad": round(theta_u, 10),
            "theta_u_mrad": round(theta_u * 1000, 6),
            "moment_yield_knm": mphi_duct["moment_yield_knm"],
            "moment_ultimate_knm": mphi_duct["moment_ultimate_knm"],
            "mu_theta": round(mu_theta, 2),
        }

    return {
        "points": points,
        "events": events,
        "ductility_rotation": ductility_rotation,
        "lp": round(lp, 2),
    }


# ---------------------------------------------------------------------------
# ASCE 41-17 / FEMA 356 Backbone Curve
# ---------------------------------------------------------------------------

# Table 10-7 Condition i — beams controlled by flexure
# Columns: (rho_ratio_max, v_ratio_max, conforming, a, b, c, IO, LS_pri, LS_sec, CP_pri, CP_sec)
_TABLE_10_7 = [
    # Conforming transverse reinforcement
    (0.0, 3.0, True,  0.025, 0.050, 0.20, 0.010, 0.020, 0.025, 0.025, 0.050),
    (0.0, 6.0, True,  0.020, 0.040, 0.20, 0.005, 0.010, 0.020, 0.020, 0.040),
    (0.5, 3.0, True,  0.020, 0.030, 0.20, 0.005, 0.010, 0.020, 0.020, 0.030),
    (0.5, 6.0, True,  0.015, 0.020, 0.20, 0.005, 0.005, 0.015, 0.015, 0.020),
    # Nonconforming transverse reinforcement
    (0.0, 3.0, False, 0.020, 0.030, 0.20, 0.005, 0.010, 0.020, 0.020, 0.030),
    (0.0, 6.0, False, 0.010, 0.015, 0.20, 0.0015, 0.005, 0.010, 0.010, 0.015),
    (0.5, 3.0, False, 0.010, 0.015, 0.20, 0.005, 0.010, 0.010, 0.010, 0.015),
    (0.5, 6.0, False, 0.005, 0.010, 0.20, 0.0015, 0.005, 0.005, 0.005, 0.010),
]


def _interp(x, x0, x1, y0, y1):
    """Linear interpolation clamped to [x0, x1]."""
    if x1 == x0:
        return y0
    t = max(0.0, min(1.0, (x - x0) / (x1 - x0)))
    return y0 + t * (y1 - y0)


def _lookup_table_10_7(rho_ratio, v_ratio, conforming):
    """
    Interpolate ASCE 41-17 Table 10-7 for beams controlled by flexure.

    Returns (a, b, c, IO, LS_pri, LS_sec, CP_pri, CP_sec).
    """
    conf = conforming
    # Filter rows by conforming
    rows = [r for r in _TABLE_10_7 if r[2] == conf]

    # Clamp inputs
    rho_ratio = max(0.0, min(0.5, rho_ratio))
    v_ratio = max(3.0, min(6.0, v_ratio))

    # 2D interpolation: first by v_ratio at rho=0 and rho=0.5, then by rho
    def get_at_rho(rho_max):
        r_lo = [r for r in rows if r[0] == rho_max and r[1] == 3.0]
        r_hi = [r for r in rows if r[0] == rho_max and r[1] == 6.0]
        if not r_lo or not r_hi:
            return None
        lo, hi = r_lo[0], r_hi[0]
        return tuple(_interp(v_ratio, 3.0, 6.0, lo[i], hi[i]) for i in range(3, 11))

    vals_0 = get_at_rho(0.0)
    vals_5 = get_at_rho(0.5)
    if vals_0 is None or vals_5 is None:
        return (0.02, 0.03, 0.2, 0.005, 0.01, 0.02, 0.02, 0.03)

    result = tuple(_interp(rho_ratio, 0.0, 0.5, vals_0[i], vals_5[i])
                   for i in range(8))
    return result


def asce41_backbone_beam(
    my_knm, theta_y_rad,
    fc, fy, b, d,
    as_tension, as_compression=0.0,
    stirrup_spacing=150.0, db_stirrup=10.0,
    v_demand=0.0,
    code="asce41",
):
    """
    Generate ASCE 41-17 (or FEMA 356) backbone curve for RC beams.

    Constructs the generalized force-deformation relation per Fig. 10-1
    with modeling parameters from Table 10-7 (Condition i — flexure).

    Parameters
    ----------
    my_knm : float
        Yield moment from M-phi analysis (kN-m).
    theta_y_rad : float
        Yield rotation (radians).
    fc : float
        Concrete compressive strength (MPa).
    fy : float
        Steel yield strength (MPa).
    b : float
        Section width (mm).
    d : float
        Effective depth (mm).
    as_tension : float
        Tension reinforcement area (mm^2).
    as_compression : float
        Compression reinforcement area (mm^2).
    stirrup_spacing : float
        Transverse reinforcement spacing (mm).
    db_stirrup : float
        Stirrup bar diameter (mm).
    v_demand : float
        Shear demand at section (kN). If 0, computed from My/Ln.
    code : str
        'asce41' or 'fema356'.

    Returns
    -------
    dict
        backbone_points, params, acceptance, calculation_steps
    """
    steps = []

    # Step 1: Beta1
    if fc <= 28:
        beta1 = 0.85
    elif fc < 55:
        beta1 = max(0.85 - 0.05 * (fc - 28) / 7, 0.65)
    else:
        beta1 = 0.65
    steps.append({
        "step": 1, "title": "Whitney Stress Block Factor",
        "eq": "beta1", "value": round(beta1, 4),
        "detail": f"beta1 = {beta1:.4f} (f'c = {fc} MPa)",
    })

    # Step 2: Balanced reinforcement ratio
    rho_bal = 0.85 * beta1 * (fc / fy) * (0.003 / (0.003 + fy / 200000))
    steps.append({
        "step": 2, "title": "Balanced Reinforcement Ratio",
        "eq": "rho_bal = 0.85 * beta1 * (fc/fy) * (ecu / (ecu + ey))",
        "value": round(rho_bal, 6),
        "detail": f"rho_bal = {rho_bal:.6f}",
    })

    # Step 3: Reinforcement ratios
    rho = as_tension / (b * d)
    rho_prime = as_compression / (b * d) if as_compression > 0 else 0
    rho_ratio = (rho - rho_prime) / rho_bal if rho_bal > 0 else 0
    rho_ratio = max(0.0, rho_ratio)
    steps.append({
        "step": 3, "title": "Reinforcement Ratio",
        "eq": "(rho - rho') / rho_bal",
        "value": round(rho_ratio, 4),
        "detail": (f"rho = {rho:.6f}, rho' = {rho_prime:.6f}, "
                   f"(rho - rho')/rho_bal = {rho_ratio:.4f}"),
    })

    # Step 4: Shear ratio
    if v_demand > 0:
        v_ratio = (v_demand * 1000) / (b * d * math.sqrt(fc))
    else:
        v_ratio = 3.0  # default low shear
    v_ratio_clamped = max(3.0, min(6.0, v_ratio))
    steps.append({
        "step": 4, "title": "Shear Stress Ratio",
        "eq": "V / (bw * d * sqrt(f'c))",
        "value": round(v_ratio, 4),
        "detail": f"V_ratio = {v_ratio:.4f} (clamped to {v_ratio_clamped:.4f})",
    })

    # Step 5: Transverse reinforcement classification
    conforming = stirrup_spacing <= d / 2
    conf_label = "Conforming (C)" if conforming else "Nonconforming (NC)"
    steps.append({
        "step": 5, "title": "Transverse Reinforcement",
        "eq": "s <= d/2 → Conforming",
        "value": conf_label,
        "detail": f"s = {stirrup_spacing} mm, d/2 = {d/2:.1f} mm → {conf_label}",
    })

    # Step 6: Table 10-7 lookup
    a, b_val, c, io, ls_pri, ls_sec, cp_pri, cp_sec = _lookup_table_10_7(
        rho_ratio, v_ratio_clamped, conforming,
    )
    table_ref = "ASCE 41-17 Table 10-7" if code == "asce41" else "FEMA 356 Table 6-7"
    steps.append({
        "step": 6, "title": f"Modeling Parameters ({table_ref})",
        "eq": "Interpolated from table",
        "value": f"a={a:.4f}, b={b_val:.4f}, c={c:.2f}",
        "detail": (f"a = {a:.4f} rad, b = {b_val:.4f} rad, c = {c:.2f} "
                   f"(Condition i — flexure-controlled)"),
    })

    # Step 7: Acceptance criteria
    steps.append({
        "step": 7, "title": "Acceptance Criteria",
        "eq": "Plastic rotation limits (rad)",
        "value": f"IO={io:.4f}, LS(Pri)={ls_pri:.4f}, CP(Pri)={cp_pri:.4f}",
        "detail": (f"IO = {io:.4f}, LS(Pri) = {ls_pri:.4f}, "
                   f"LS(Sec) = {ls_sec:.4f}, CP(Pri) = {cp_pri:.4f}, "
                   f"CP(Sec) = {cp_sec:.4f} rad"),
    })

    # Step 8: Construct backbone points
    # ASCE 41-17 Fig. 10-1: A-B-C-D-E
    theta_y = theta_y_rad
    backbone_points = [
        {"label": "A", "theta_rad": 0, "theta_mrad": 0,
         "moment_knm": 0, "q_ratio": 0},
        {"label": "B", "theta_rad": round(theta_y, 8),
         "theta_mrad": round(theta_y * 1000, 4),
         "moment_knm": round(my_knm, 2), "q_ratio": 1.0},
        {"label": "C", "theta_rad": round(theta_y + a, 8),
         "theta_mrad": round((theta_y + a) * 1000, 4),
         "moment_knm": round(my_knm, 2), "q_ratio": 1.0},
        {"label": "D", "theta_rad": round(theta_y + a, 8),
         "theta_mrad": round((theta_y + a) * 1000, 4),
         "moment_knm": round(my_knm * c, 2), "q_ratio": c},
        {"label": "E", "theta_rad": round(theta_y + b_val, 8),
         "theta_mrad": round((theta_y + b_val) * 1000, 4),
         "moment_knm": round(my_knm * c, 2), "q_ratio": c},
    ]

    # Acceptance criteria as absolute rotation (theta_y + plastic_rotation)
    acceptance = {
        "io_rad": round(theta_y + io, 8),
        "io_mrad": round((theta_y + io) * 1000, 4),
        "ls_pri_rad": round(theta_y + ls_pri, 8),
        "ls_pri_mrad": round((theta_y + ls_pri) * 1000, 4),
        "ls_sec_rad": round(theta_y + ls_sec, 8),
        "ls_sec_mrad": round((theta_y + ls_sec) * 1000, 4),
        "cp_pri_rad": round(theta_y + cp_pri, 8),
        "cp_pri_mrad": round((theta_y + cp_pri) * 1000, 4),
        "cp_sec_rad": round(theta_y + cp_sec, 8),
        "cp_sec_mrad": round((theta_y + cp_sec) * 1000, 4),
        "io_plastic": round(io, 4),
        "ls_pri_plastic": round(ls_pri, 4),
        "cp_pri_plastic": round(cp_pri, 4),
    }

    return {
        "backbone_points": backbone_points,
        "params": {
            "a": round(a, 4), "b": round(b_val, 4), "c": round(c, 2),
            "rho": round(rho, 6), "rho_prime": round(rho_prime, 6),
            "rho_bal": round(rho_bal, 6), "rho_ratio": round(rho_ratio, 4),
            "v_ratio": round(v_ratio, 4),
            "conforming": conforming, "conf_label": conf_label,
            "beta1": round(beta1, 4),
            "theta_y_rad": round(theta_y, 8),
            "my_knm": round(my_knm, 2),
            "code": code, "table_ref": table_ref,
        },
        "acceptance": acceptance,
        "calculation_steps": steps,
    }
