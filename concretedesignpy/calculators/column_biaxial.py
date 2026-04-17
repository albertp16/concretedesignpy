# MIT License
# Copyright (c) Albert Pamonag Engineering Consultancy

"""
Column Biaxial P-M-M Interaction Surface
==========================================

Generates the 3D interaction surface (P, Mx, My) for rectangular
RC columns using a fiber-based strain compatibility approach.

For each neutral axis angle (0-180 deg) and depth, the section is
discretized into a 2D fiber grid. Whitney stress block is applied
to fibers on the compression side of the inclined neutral axis.

Reference:
    Anwar & Najam (2017), Ch. 3 & 5
    NSCP 2015 / ACI 318-19
"""

import math
import numpy as np

from concretedesignpy.calculators.column_interaction import (
    _beta_one,
    _phi_factor,
)


def _generate_bar_coords_2d(b, h, cover, tie, d_bar, nx, ny):
    """
    Generate 2D bar coordinates centered on section centroid (0, 0).

    Parameters
    ----------
    b : float
        Section width (mm).
    h : float
        Section depth (mm).
    cover : float
        Clear cover to tie (mm).
    tie : float
        Tie diameter (mm).
    d_bar : float
        Main bar diameter (mm).
    nx : int
        Bars along width (top and bottom faces).
    ny : int
        Bars along depth (left and right faces).

    Returns
    -------
    list of (float, float)
        Bar coordinates (x, y) from centroid.
    """
    offset = cover + tie + d_bar / 2.0
    x_left = -b / 2.0 + offset
    x_right = b / 2.0 - offset
    y_bot = -h / 2.0 + offset
    y_top = h / 2.0 - offset

    bars = []

    # Top face: nx bars
    if nx >= 2:
        for i in range(nx):
            x = x_left + i * (x_right - x_left) / (nx - 1)
            bars.append((x, y_top))
    elif nx == 1:
        bars.append((0.0, y_top))

    # Bottom face: nx bars
    if nx >= 2:
        for i in range(nx):
            x = x_left + i * (x_right - x_left) / (nx - 1)
            bars.append((x, y_bot))
    elif nx == 1:
        bars.append((0.0, y_bot))

    # Left and right sides: (ny - 2) intermediate bars each
    if ny > 2:
        for j in range(1, ny - 1):
            y = y_bot + j * (y_top - y_bot) / (ny - 1)
            bars.append((x_left, y))
            bars.append((x_right, y))

    return bars


def generate_biaxial_diagram(
    fc, fy, b, h, bar_coords_2d, bar_areas,
    cover=40, confinement="tied",
    n_angles=24, n_c_values=32,
    nx_fibers=20, ny_fibers=20,
    ecu=0.003,
):
    """
    Generate P-M-M biaxial interaction surface using fiber approach.

    Parameters
    ----------
    fc : float
        Concrete compressive strength (MPa).
    fy : float
        Steel yield strength (MPa).
    b : float
        Section width (mm).
    h : float
        Section depth (mm).
    bar_coords_2d : list of (float, float)
        Bar coordinates (x, y) from centroid.
    bar_areas : list of float
        Bar areas (mm^2).
    cover : float
        Clear cover (mm).
    confinement : str
        'tied' or 'spiral'.
    n_angles : int
        Number of NA angle divisions (0 to 180 deg).
    n_c_values : int
        Number of NA depth divisions per angle.
    nx_fibers : int
        Fiber grid divisions along width.
    ny_fibers : int
        Fiber grid divisions along depth.
    ecu : float
        Ultimate concrete strain.

    Returns
    -------
    dict
        surface_points, section_info, pure_compression_kn, etc.
    """
    es_mod = 200000.0
    beta1 = _beta_one(fc)
    ag = b * h
    ast_total = sum(bar_areas)
    a_bar = math.pi / 4.0  # not used, bar_areas provided

    # Pure compression / tension
    po = 0.85 * fc * (ag - ast_total) + fy * ast_total
    k1 = 0.85 if confinement == "spiral" else 0.80
    pn_max = k1 * po
    pt = -fy * ast_total

    phi_c = 0.75 if confinement == "spiral" else 0.65

    # Fiber grid centered on (0, 0)
    dx = b / nx_fibers
    dy = h / ny_fibers
    cx = np.linspace(-b / 2 + dx / 2, b / 2 - dx / 2, nx_fibers)
    cy = np.linspace(-h / 2 + dy / 2, h / 2 - dy / 2, ny_fibers)
    fgx, fgy = np.meshgrid(cx, cy)
    fiber_x = fgx.ravel()
    fiber_y = fgy.ravel()
    fiber_area = dx * dy

    # Bar arrays
    bar_x = np.array([c[0] for c in bar_coords_2d])
    bar_y = np.array([c[1] for c in bar_coords_2d])
    bar_a = np.array(bar_areas)
    n_bars = len(bar_areas)

    # Section corners for extreme compression calculation
    corners_x = np.array([-b / 2, b / 2, b / 2, -b / 2])
    corners_y = np.array([-h / 2, -h / 2, h / 2, h / 2])

    # Angle sweep (0 to 180 degrees, inclusive)
    angles = np.linspace(0, math.pi, n_angles + 1)

    surface_points = []

    for theta in angles:
        cos_t = math.cos(theta)
        sin_t = math.sin(theta)

        # Projection of all corners along NA normal
        corner_proj = corners_x * cos_t + corners_y * sin_t
        d_max = float(np.max(corner_proj))
        d_min = float(np.min(corner_proj))
        section_depth = d_max - d_min

        # Fiber projections (distance from extreme compression face)
        fiber_proj = fiber_x * cos_t + fiber_y * sin_t
        fiber_dist = d_max - fiber_proj  # 0 at compression face

        # Bar projections
        bar_proj = bar_x * cos_t + bar_y * sin_t
        bar_dist = d_max - bar_proj

        # Sweep NA depth
        c_min_val = 1.0
        c_max_val = section_depth * 1.5

        for ci in range(n_c_values + 1):
            c = c_max_val - ci * (c_max_val - c_min_val) / n_c_values
            a = beta1 * c

            # ── Concrete fibers ──
            # Whitney block: fibers within depth a get 0.85*f'c
            in_block = fiber_dist <= a
            concrete_forces = np.where(in_block, 0.85 * fc * fiber_area, 0.0)

            cc_total = float(np.sum(concrete_forces))
            mnx_concrete = float(np.sum(concrete_forces * fiber_y))
            mny_concrete = float(np.sum(concrete_forces * fiber_x))

            # ── Steel bars ──
            bar_strains = ecu * (c - bar_dist) / c  # positive = compression
            bar_stresses = np.clip(bar_strains * es_mod, -fy, fy)

            # Subtract displaced concrete where bar is in Whitney block
            fc_at_bar = np.where(bar_dist <= a, 0.85 * fc, 0.0)
            bar_forces = bar_a * (bar_stresses - fc_at_bar)

            fs_total = float(np.sum(bar_forces))
            mnx_steel = float(np.sum(bar_forces * bar_y))
            mny_steel = float(np.sum(bar_forces * bar_x))

            # ── Totals ──
            pn = cc_total + fs_total  # N
            mnx = mnx_concrete + mnx_steel  # N-mm
            mny = mny_concrete + mny_steel  # N-mm

            # Phi factor from extreme tension bar
            max_bar_dist = float(np.max(bar_dist))
            es_extreme = ecu * (c - max_bar_dist) / c
            net_tensile = -es_extreme if es_extreme < 0 else 0.0
            phi, classify = _phi_factor(net_tensile, fy, es_mod, confinement)

            # Cap Pn
            pn_capped = min(pn, pn_max) if pn > 0 else pn

            # Factored
            pu = phi * pn_capped
            mux = phi * mnx
            muy = phi * mny

            surface_points.append({
                "theta_deg": round(math.degrees(theta), 1),
                "c": round(c, 1),
                "pn": round(pn / 1000, 2),
                "pu": round(pu / 1000, 2),
                "mnx": round(mnx / 1e6, 2),
                "mny": round(mny / 1e6, 2),
                "mux": round(mux / 1e6, 2),
                "muy": round(muy / 1e6, 2),
                "phi": round(phi, 4),
                "classification": classify,
            })

    # Pure tension point
    pt_kn = pt / 1000
    surface_points.append({
        "theta_deg": 0, "c": 0,
        "pn": round(pt_kn, 2), "pu": round(0.9 * pt_kn, 2),
        "mnx": 0, "mny": 0, "mux": 0, "muy": 0,
        "phi": 0.9, "classification": "tension-controlled",
    })

    # ── Fiber plot at theta=0, balanced depth ──
    # Find a representative c at theta=0 near balanced condition
    theta0 = 0.0
    cos0, sin0 = math.cos(theta0), math.sin(theta0)
    corner_proj0 = corners_x * cos0 + corners_y * sin0
    d_max0 = float(np.max(corner_proj0))
    fiber_dist0 = d_max0 - (fiber_x * cos0 + fiber_y * sin0)
    bar_dist0 = d_max0 - (bar_x * cos0 + bar_y * sin0)
    # Use c at about 60% of section depth (near balanced)
    section_d0 = d_max0 - float(np.min(corner_proj0))
    c_rep = section_d0 * 0.6
    a_rep = beta1 * c_rep

    # Fiber stresses
    in_block = fiber_dist0 <= a_rep
    fiber_stresses_plot = np.where(in_block, 0.85 * fc, 0.0)
    fiber_strains_plot = ecu * (c_rep - fiber_dist0) / c_rep

    # Bar stresses
    bar_strains_rep = ecu * (c_rep - bar_dist0) / c_rep
    bar_stresses_rep = np.clip(bar_strains_rep * es_mod, -fy, fy)

    fiber_plot = {
        "fiber_x": [round(float(x), 1) for x in fiber_x],
        "fiber_y": [round(float(y), 1) for y in fiber_y],
        "stresses": [round(float(s), 4) for s in fiber_stresses_plot],
        "strains": [round(float(s), 8) for s in fiber_strains_plot],
        "nx": nx_fibers,
        "ny": ny_fibers,
        "dx": round(dx, 2),
        "dy": round(dy, 2),
        "c": round(c_rep, 1),
        "theta_deg": 0,
        "bar_x": [round(float(x), 1) for x in bar_x],
        "bar_y": [round(float(y), 1) for y in bar_y],
        "bar_stresses": [round(float(s), 2) for s in bar_stresses_rep],
    }

    return {
        "surface_points": surface_points,
        "fiber_plot": fiber_plot,
        "section_info": {
            "b": b, "h": h,
            "ag": round(ag, 0),
            "ast_total": round(ast_total, 2),
            "rho": round(ast_total / ag, 6),
            "n_bars": n_bars,
            "beta1": round(beta1, 4),
            "confinement": confinement,
            "bar_coords": [(round(x, 1), round(y, 1))
                           for x, y in bar_coords_2d],
            "bar_areas": [round(a, 1) for a in bar_areas],
        },
        "pure_compression_kn": round(po / 1000, 2),
        "pure_tension_kn": round(pt_kn, 2),
        "phi_pn_max_kn": round(phi_c * pn_max / 1000, 2),
        "n_angles": n_angles,
        "n_c_values": n_c_values,
    }


def extract_contour_at_pu(biaxial_result, pu_level):
    """
    Extract Mx-My contour at a specified Pu level.

    Interpolates the interaction surface at the given factored
    axial load to produce a closed Mx-My envelope.

    Parameters
    ----------
    biaxial_result : dict
        Result from generate_biaxial_diagram().
    pu_level : float
        Factored axial load (kN) at which to slice.

    Returns
    -------
    list of dict
        Contour points {theta_deg, mux, muy}.
    """
    pts = biaxial_result["surface_points"]

    # Group points by theta
    by_theta = {}
    for p in pts:
        theta = p["theta_deg"]
        if theta not in by_theta:
            by_theta[theta] = []
        by_theta[theta].append(p)

    contour = []
    for theta in sorted(by_theta.keys()):
        group = sorted(by_theta[theta], key=lambda p: -p["pu"])

        # Find two points bracketing pu_level
        mux_interp = None
        muy_interp = None
        for i in range(len(group) - 1):
            p1 = group[i]
            p2 = group[i + 1]
            if p1["pu"] >= pu_level >= p2["pu"]:
                if abs(p1["pu"] - p2["pu"]) < 0.01:
                    mux_interp = p1["mux"]
                    muy_interp = p1["muy"]
                else:
                    t = (pu_level - p2["pu"]) / (p1["pu"] - p2["pu"])
                    mux_interp = p2["mux"] + t * (p1["mux"] - p2["mux"])
                    muy_interp = p2["muy"] + t * (p1["muy"] - p2["muy"])
                break

        if mux_interp is not None:
            contour.append({
                "theta_deg": theta,
                "mux": round(mux_interp, 2),
                "muy": round(muy_interp, 2),
            })

    # Mirror to 180-360 for full contour
    mirrored = []
    for p in reversed(contour):
        if p["theta_deg"] in (0, 180):
            continue
        mirrored.append({
            "theta_deg": round(p["theta_deg"] + 180, 1),
            "mux": round(-p["mux"], 2),
            "muy": round(-p["muy"], 2),
        })
    contour.extend(mirrored)

    # Close the contour
    if contour:
        contour.append(contour[0])

    return contour


def check_biaxial_capacity(biaxial_result, pu_demand, mux_demand, muy_demand):
    """
    Check if a (Pu, Mux, Muy) demand falls inside the interaction surface.

    Uses the contour method: extracts the Mx-My envelope at the demand
    Pu level, then checks if the demand moment vector is inside.

    Parameters
    ----------
    biaxial_result : dict
        Result from generate_biaxial_diagram().
    pu_demand : float
        Factored axial load demand (kN).
    mux_demand : float
        Factored moment demand about x-axis (kN-m).
    muy_demand : float
        Factored moment demand about y-axis (kN-m).

    Returns
    -------
    dict
        demand_mn, capacity_mn, dc_ratio, status, demand_angle_deg
    """
    demand_mn = math.sqrt(mux_demand ** 2 + muy_demand ** 2)
    demand_angle = math.atan2(muy_demand, mux_demand)  # radians

    contour = extract_contour_at_pu(biaxial_result, pu_demand)

    if not contour:
        # Demand Pu outside surface range
        return {
            "demand_mn": round(demand_mn, 2),
            "capacity_mn": 0,
            "dc_ratio": 999,
            "status": "exceeds axial capacity",
            "demand_angle_deg": round(math.degrees(demand_angle), 1),
        }

    # Find capacity at the demand angle by interpolating contour
    # Convert contour points to polar (angle, radius)
    best_capacity = 0
    for i in range(len(contour) - 1):
        p1 = contour[i]
        p2 = contour[i + 1]
        a1 = math.atan2(p1["muy"], p1["mux"])
        a2 = math.atan2(p2["muy"], p2["mux"])
        r1 = math.sqrt(p1["mux"] ** 2 + p1["muy"] ** 2)
        r2 = math.sqrt(p2["mux"] ** 2 + p2["muy"] ** 2)

        # Check if demand_angle is between a1 and a2
        # Normalize angles
        da = demand_angle
        if a2 < a1:
            a1, a2, r1, r2 = a2, a1, r2, r1
        if a1 <= da <= a2 and abs(a2 - a1) < math.pi:
            if abs(a2 - a1) < 1e-8:
                best_capacity = max(r1, r2)
            else:
                t = (da - a1) / (a2 - a1)
                best_capacity = r1 + t * (r2 - r1)
            break

    # Fallback: find closest contour point
    if best_capacity == 0 and contour:
        min_angle_diff = math.pi * 2
        for p in contour:
            a_p = math.atan2(p["muy"], p["mux"])
            diff = abs(a_p - demand_angle)
            if diff > math.pi:
                diff = 2 * math.pi - diff
            if diff < min_angle_diff:
                min_angle_diff = diff
                best_capacity = math.sqrt(p["mux"] ** 2 + p["muy"] ** 2)

    dc_ratio = demand_mn / best_capacity if best_capacity > 0 else 999

    return {
        "demand_mn": round(demand_mn, 2),
        "capacity_mn": round(best_capacity, 2),
        "dc_ratio": round(dc_ratio, 4),
        "status": "OK" if dc_ratio <= 1.0 else "FAIL",
        "demand_angle_deg": round(math.degrees(demand_angle), 1),
    }
