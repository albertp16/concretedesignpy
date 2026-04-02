# MIT License
# Copyright (c) Albert Pamonag Engineering Consultancy

"""
Column P-M Interaction Diagram Generator
==========================================

Strain-compatibility analysis for rectangular RC columns.
Aligned with ShortCol methodology (yakpol.net).

Reference: NSCP 2015 / ACI 318-19
"""

import math


def _beta_one(fc):
    """Whitney stress block factor β₁.

    Parameters
    ----------
    fc : float
        Concrete compressive strength (MPa).
    """
    if fc <= 28:
        return 0.85
    elif fc < 55:
        return max(0.85 - 0.05 * (fc - 28) / 7, 0.65)
    else:
        return 0.65


def _phi_factor(es_tension, fy, es_mod, confinement="tied"):
    """Strength reduction factor per ACI 318-19 Table 21.2.2.

    Parameters
    ----------
    es_tension : float
        Net tensile strain in extreme tension steel.
    fy : float
        Steel yield strength (MPa).
    es_mod : float
        Steel modulus of elasticity (MPa).
    confinement : str
        'tied' or 'spiral'.

    Returns
    -------
    tuple
        (phi, classification)
    """
    ey = fy / es_mod
    ecl = 0.002   # compression-controlled strain limit
    etl = 0.005   # tension-controlled strain limit

    if confinement == "spiral":
        phi_c = 0.75
    else:
        phi_c = 0.65

    phi_b = 0.90

    et = abs(es_tension)
    if et >= etl:
        return phi_b, "tension-controlled"
    elif et <= ecl:
        return phi_c, "compression-controlled"
    else:
        phi = phi_c + (phi_b - phi_c) * (et - ecl) / (etl - ecl)
        return phi, "transition"


def _generate_bar_coords(b, h, cover, tie, d_bar, n_bars, n_bars_side=0):
    """Generate bar coordinates for rectangular column.

    Bars are placed around the perimeter. Coordinates are distances
    from the compression face (top) measured downward.

    Parameters
    ----------
    b : float
        Column width (mm).
    h : float
        Column depth in bending direction (mm).
    cover : float
        Clear cover to stirrup (mm).
    tie : float
        Stirrup/tie diameter (mm).
    d_bar : float
        Main bar diameter (mm).
    n_bars : int
        Total number of bars.
    n_bars_side : int
        Bars on each side face (excluding corners). If 0, auto-calculated.

    Returns
    -------
    list of float
        Distance from compression face (top) for each bar.
    """
    d_prime = cover + tie + d_bar / 2.0
    d_eff = h - d_prime

    if n_bars <= 4:
        # Corner bars only: 2 top, 2 bottom
        n_top = n_bars // 2
        n_bot = n_bars - n_top
        return [d_prime] * n_top + [d_eff] * n_bot

    # Number of bars per face
    if n_bars_side == 0:
        # Auto: distribute evenly around perimeter
        # Top face + bottom face get remaining bars
        n_bars_side = max((n_bars - 4) // 4, 0)
        # Remaining for top/bottom
        remaining = n_bars - 4 - 4 * n_bars_side
        n_top_extra = remaining // 2
        n_bot_extra = remaining - n_top_extra
    else:
        n_top_extra = 0
        n_bot_extra = 0
        # Recalculate: total = 2*(n_top_face) + 2*(n_bars_side) where top_face includes corners
        pass

    # Top face bars (at d_prime from top)
    n_top_face = 2 + n_top_extra  # 2 corners + extra
    # Bottom face bars (at d_eff from top)
    n_bot_face = 2 + n_bot_extra

    bars = []
    # Top face
    bars.extend([d_prime] * n_top_face)

    # Side bars - evenly spaced between d_prime and d_eff
    total_side = n_bars - n_top_face - n_bot_face
    n_side_per = total_side // 2 if total_side > 0 else 0
    if n_side_per > 0:
        for j in range(1, n_side_per + 1):
            ds = d_prime + j * (d_eff - d_prime) / (n_side_per + 1)
            bars.append(ds)  # left side
            bars.append(ds)  # right side

    # Bottom face
    bars.extend([d_eff] * n_bot_face)

    return bars


def generate_interaction_diagram(
    fc, fy, b, h, n_bars, d_bar, n_bars_side=0,
    cover=40, c_step=5, confinement="tied",
    bar_coords=None, bar_areas=None,
    n_points=32, ecu=0.003,
):
    """
    Generate P-M interaction diagram using strain compatibility.

    Follows the ShortCol methodology: for each neutral axis position,
    compute strain in every bar, cap stress at ±fy, compute concrete
    compression via Whitney stress block, sum forces and moments about
    the gross section centroid.

    Parameters
    ----------
    fc : float
        Concrete compressive strength (MPa).
    fy : float
        Steel yield strength (MPa).
    b : float
        Column width (mm).
    h : float
        Column depth in bending direction (mm).
    n_bars : int
        Total number of bars.
    d_bar : float
        Main bar diameter (mm).
    n_bars_side : int
        Bars on each side face (excluding corners).
    cover : float
        Clear cover to stirrup/tie (mm).
    c_step : float
        Step size for neutral axis sweep (mm). Used if n_points=0.
    confinement : str
        'tied' or 'spiral'.
    bar_coords : list of float or None
        Custom bar distances from compression face (mm).
        If None, auto-generated from n_bars/n_bars_side.
    bar_areas : list of float or None
        Custom bar areas (mm²). If None, all bars use π/4 * d_bar².
    n_points : int
        Number of evenly-spaced points along the interaction curve.
        Set to 0 to use c_step sweep instead.
    ecu : float
        Ultimate concrete strain (default 0.003).

    Returns
    -------
    dict
        points, pure_compression_kn, pure_tension_kn,
        balanced_point, beta1, d_eff, d_prime, bar_depths,
        confinement
    """
    es_mod = 200000.0

    # Tie diameter for auto bar layout
    tie = 10 if d_bar <= 32 else 12

    # Effective depths
    d_prime = cover + tie + d_bar / 2.0
    d_eff = h - d_prime

    beta1 = _beta_one(fc)

    # Bar layout
    if bar_coords is not None:
        bar_ds = list(bar_coords)
    else:
        bar_ds = _generate_bar_coords(b, h, cover, tie, d_bar, n_bars, n_bars_side)

    n_actual = len(bar_ds)

    if bar_areas is not None:
        bar_as = list(bar_areas)
    else:
        a_bar = math.pi * d_bar ** 2 / 4.0
        bar_as = [a_bar] * n_actual

    # Gross section properties
    ag = b * h
    ast_total = sum(bar_as)

    # Section centroid (from compression face)
    yc = h / 2.0

    # Pure compression (c → ∞): all steel yields in compression
    # Pn,max = 0.85*f'c*(Ag - Ast) + fy*Ast
    po = 0.85 * fc * (ag - ast_total) + fy * ast_total

    # Maximum factored axial (ACI 318-19 Eq. 22.4.2.1)
    if confinement == "spiral":
        k1 = 0.85
    else:
        k1 = 0.80
    pn_max = k1 * po

    # Pure tension: all steel yields in tension
    pt = -fy * ast_total

    # Generate c values for the sweep
    # Sweep from large c (pure compression) down to small c (pure tension)
    c_min = 1.0  # near pure tension
    c_max = h * 1.5  # beyond pure compression

    if n_points > 0:
        c_values = []
        for i in range(n_points + 1):
            ci = c_max - i * (c_max - c_min) / n_points
            c_values.append(ci)
    else:
        c_values = []
        ci = c_max
        while ci >= c_min:
            c_values.append(ci)
            ci -= c_step

    points = []
    balanced_point = None
    max_moment_point = None
    balanced_steel_detail = None

    for c in c_values:
        a = beta1 * c

        # --- Concrete compression force ---
        # Whitney stress block: 0.85*f'c over depth 'a', width 'b'
        a_eff = min(a, h)  # stress block can't exceed section height
        cc = 0.85 * fc * a_eff * b  # N (compression, positive)

        # Centroid of concrete compression from top
        ycc = a_eff / 2.0

        # --- Steel forces ---
        steel_forces = []
        for k in range(n_actual):
            ds_k = bar_ds[k]  # distance from compression face
            as_k = bar_as[k]

            # Strain: linear from ecu at compression face to 0 at NA
            es_k = ecu * (c - ds_k) / c  # positive = compression

            # Stress: cap at ±fy
            fs_k = es_k * es_mod
            if fs_k > fy:
                fs_k = fy
            elif fs_k < -fy:
                fs_k = -fy

            # Subtract concrete stress where bar displaces concrete
            # (only if bar is within the compression block)
            fc_at_bar = 0.85 * fc if ds_k <= a_eff else 0.0

            # Net bar force (positive = compression)
            f_bar = as_k * (fs_k - fc_at_bar)

            steel_forces.append({
                "ds": ds_k,
                "es": es_k,
                "fs": fs_k,
                "force": f_bar,
            })

        # Sum of steel forces
        fs_total = sum(sf["force"] for sf in steel_forces)

        # Total axial force (positive = compression)
        pn = cc + fs_total  # N

        # Moment about gross section centroid
        # Concrete: Cc acts at ycc from top, arm = yc - ycc
        mn_concrete = cc * (yc - ycc)

        # Steel: each bar at ds_k from top, arm = yc - ds_k
        mn_steel = sum(sf["force"] * (yc - sf["ds"]) for sf in steel_forces)

        mn = mn_concrete + mn_steel  # N-mm

        # Extreme tension bar strain (for phi calculation)
        # Find bar farthest from compression face
        max_ds = max(bar_ds)
        es_extreme = ecu * (c - max_ds) / c  # positive = compression
        # Net tensile strain (negative es_extreme means tension)
        net_tensile_strain = -es_extreme if es_extreme < 0 else 0

        # Phi factor
        phi, classify = _phi_factor(net_tensile_strain, fy, es_mod, confinement)

        # Cap Pn at Pn,max for compression
        if pn > 0:
            pn_capped = min(pn, pn_max)
        else:
            pn_capped = pn

        # Factored values
        pu = phi * pn_capped
        mu = phi * mn

        # Convert to kN and kN-m
        pn_kn = pn / 1000.0
        pu_kn = pu / 1000.0
        mn_knm = mn / 1e6
        mu_knm = mu / 1e6

        # Per-bar detail for report
        bar_detail = []
        for k in range(n_actual):
            sf = steel_forces[k]
            bar_detail.append({
                "bar": k + 1,
                "as_mm2": round(bar_as[k], 1),
                "ds_mm": round(sf["ds"], 2),
                "strain": round(sf["es"], 6),
                "fs_mpa": round(sf["fs"], 2),
                "force_kn": round(sf["force"] / 1000.0, 2),
            })

        point = {
            "c": round(c, 2),
            "a": round(min(a, h), 2),
            "pn": round(pn_kn, 2),
            "pu": round(pu_kn, 2),
            "mn": round(mn_knm, 2),
            "mu": round(mu_knm, 2),
            "phi": round(phi, 4),
            "classification": classify,
            "es_tension": round(net_tensile_strain, 6),
            "cc_kn": round(cc / 1000.0, 2),
            "fs_total_kn": round(fs_total / 1000.0, 2),
            "steel_forces": bar_detail,
        }
        points.append(point)

        # Track balanced point (tension steel just yielding)
        ey = fy / es_mod
        if balanced_point is None and net_tensile_strain >= ey:
            balanced_point = point
            balanced_steel_detail = bar_detail

        # Track maximum moment point
        if max_moment_point is None or abs(mu_knm) > abs(max_moment_point["mu"]):
            max_moment_point = point

    # Add pure tension point
    pt_kn = pt / 1000.0
    points.append({
        "c": 0,
        "a": 0,
        "pn": round(pt_kn, 2),
        "pu": round(0.9 * pt_kn, 2),
        "mn": 0,
        "mu": 0,
        "phi": 0.9,
        "classification": "tension-controlled",
        "es_tension": 999,
    })

    # Add pure compression point at top
    po_kn = po / 1000.0
    pnmax_kn = pn_max / 1000.0
    phi_c = 0.75 if confinement == "spiral" else 0.65

    # Gross moment of inertia
    ig = b * h ** 3 / 12.0

    return {
        "points": points,
        "pure_compression_kn": round(po_kn, 2),
        "pure_tension_kn": round(pt_kn, 2),
        "phi_pn_max_kn": round(phi_c * pnmax_kn, 2),
        "balanced_point": balanced_point,
        "balanced_steel_detail": balanced_steel_detail,
        "max_moment_point": max_moment_point,
        "beta1": round(beta1, 4),
        "d_eff": round(d_eff, 2),
        "d_prime": round(d_prime, 2),
        "bar_depths": [round(d, 2) for d in bar_ds],
        "bar_areas": [round(a, 1) for a in bar_as],
        "n_bars_actual": n_actual,
        "ast_total": round(ast_total, 2),
        "ag": round(ag, 0),
        "ig": round(ig, 0),
        "rho": round(ast_total / ag, 6),
        "confinement": confinement,
    }


def check_capacity(interaction_result, pu_demand, mu_demand):
    """
    Check if a given (Pu, Mu) demand falls inside the interaction diagram.

    Given an axial load Pu (kN), find the corresponding φMn capacity
    by linear interpolation of the interaction curve.

    Parameters
    ----------
    interaction_result : dict
        Result from generate_interaction_diagram().
    pu_demand : float
        Factored axial load demand (kN), positive = compression.
    mu_demand : float
        Factored moment demand (kN-m), always positive.

    Returns
    -------
    dict
        phi_mn_capacity, dc_ratio, status
    """
    points = interaction_result["points"]

    # Sort points by pu descending
    sorted_pts = sorted(points, key=lambda p: -p["pu"])

    # Find two points that bracket the demand Pu
    phi_mn = None
    for i in range(len(sorted_pts) - 1):
        p1 = sorted_pts[i]
        p2 = sorted_pts[i + 1]
        if p1["pu"] >= pu_demand >= p2["pu"]:
            # Linear interpolation
            if abs(p1["pu"] - p2["pu"]) < 0.01:
                phi_mn = max(abs(p1["mu"]), abs(p2["mu"]))
            else:
                t = (pu_demand - p2["pu"]) / (p1["pu"] - p2["pu"])
                phi_mn = abs(p2["mu"]) + t * (abs(p1["mu"]) - abs(p2["mu"]))
            break

    if phi_mn is None:
        # Demand outside range
        if pu_demand > sorted_pts[0]["pu"]:
            return {
                "phi_mn_capacity": 0,
                "dc_ratio": 999,
                "status": "exceeds max axial capacity",
            }
        else:
            return {
                "phi_mn_capacity": 0,
                "dc_ratio": 999,
                "status": "exceeds tensile capacity",
            }

    dc_ratio = mu_demand / phi_mn if phi_mn > 0 else 999

    return {
        "phi_mn_capacity": round(phi_mn, 2),
        "dc_ratio": round(dc_ratio, 4),
        "status": "OK" if dc_ratio <= 1.0 else "FAIL",
    }
