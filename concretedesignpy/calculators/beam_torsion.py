# MIT License
# Copyright (c) Albert Pamonag Engineering Consultancy

"""
Beam Torsion Design Checker
=============================

Performs torsion design checks per NSCP 2015 / ACI 318-19.

Checks: threshold, cracking, crushing limit, combined stress,
transverse reinforcement, spacing, and longitudinal steel.
"""

import math


def torsion_design(width, height, cover, db, tf, beff, phi_torsion,
                   fc, fy, tu, vc, ds, smax_shear, s_actual, av, s):
    """
    Perform torsion design checks.

    Parameters
    ----------
    width : float
        Beam width (mm).
    height : float
        Beam height (mm).
    cover : float
        Concrete cover (mm).
    db : float
        Main bar diameter (mm).
    tf : float
        Flange thickness (mm), 0 for rectangular beams.
    beff : float
        Effective flange width (mm), equals width for rectangular.
    phi_torsion : float
        Strength reduction factor for torsion.
    fc : float
        Concrete compressive strength (MPa).
    fy : float
        Steel yield strength (MPa).
    tu : float
        Factored torsion (kN-m).
    vc : float
        Shear carried by concrete (kN).
    ds : float
        Stirrup diameter (mm).
    smax_shear : float
        Maximum shear stirrup spacing (mm).
    s_actual : float
        Actual stirrup spacing (mm).
    av : float
        Total shear reinforcement area (mm^2).
    s : float
        Stirrup spacing for shear (mm).

    Returns
    -------
    dict
        Keys: checks (list of dicts), torsion_threshold, cracking_moment,
        at_per_s, smax, al_required, status
    """
    depth = height - cover - ds - db / 2.0

    # Core geometry
    x1 = width - 2 * (cover + db / 2.0)
    y1 = height - 2 * (cover + db / 2.0)
    aoh = x1 * y1
    ph = 2 * (x1 + y1)

    # Outer cross-section area and perimeter
    if tf > 0:
        acp = width * (height - tf) + beff * tf
        pcp = 2 * (beff + height)
    else:
        acp = width * height
        pcp = 2 * (width + height)

    checks = []

    # 1. Torsional threshold
    phi_tth = (phi_torsion * 0.083 * math.sqrt(fc) / 10.0
               * (acp ** 2) / pcp / 100.0)
    torsion_needed = tu >= phi_tth
    checks.append({
        "check": "torsional_threshold",
        "value": round(phi_tth, 4),
        "tu": tu,
        "required": torsion_needed,
        "status": "MUST CONSIDER" if torsion_needed else "MAY NEGLECT",
    })

    # 2. Cracking torsional moment
    phi_tcr = 4.0 * phi_tth
    min_at_s = max(0.031 * math.sqrt(fc), 0.175) * width / fy
    checks.append({
        "check": "cracking_torsion",
        "phi_tcr": round(phi_tcr, 4),
        "min_at_per_s": round(min_at_s, 6),
    })

    # 3. Combined stress / crushing limit
    crushing_limit = phi_torsion * (
        (vc / (width * depth) * 10.0) + 0.66 * math.sqrt(fc)
    )
    combined_stress = math.sqrt(
        (vc / (width * depth) * 10.0) ** 2
        + (tu * 100.0 * ph / (1.7 * (aoh ** 2) * 10.0)) ** 2
    )
    crushing_ok = combined_stress < crushing_limit
    checks.append({
        "check": "crushing_limit",
        "crushing_limit": round(crushing_limit, 4),
        "combined_stress": round(combined_stress, 4),
        "status": "OK" if crushing_ok else "FAIL",
    })

    # 4. Transverse reinforcement
    at_per_s = tu * 100.0 / (2.0 * phi_torsion * 0.85 * aoh * (fy / 10.0))
    total_transverse = (av / s) / 2.0 + at_per_s
    transverse_ok = total_transverse >= min_at_s
    checks.append({
        "check": "transverse_reinforcement",
        "at_per_s": round(at_per_s, 6),
        "total_transverse": round(total_transverse, 6),
        "min_at_per_s": round(min_at_s, 6),
        "status": "OK" if transverse_ok else "FAIL",
    })

    # 5. Stirrup spacing
    smax = min(smax_shear, ph / 8.0, 305.0)
    spacing_ok = s_actual <= smax
    checks.append({
        "check": "stirrup_spacing",
        "smax": round(smax, 2),
        "s_actual": s_actual,
        "status": "OK" if spacing_ok else "FAIL",
    })

    # 6. Longitudinal reinforcement
    min_al = (0.42 * math.sqrt(fc) * acp / fy
              - max(at_per_s, 0.175 * width / fy) * ph)
    req_al = at_per_s * ph
    al_required = max(min_al, req_al)
    bar_area = math.pi * (db / 2.0) ** 2
    bar_ok = bar_area > al_required / 4.0
    checks.append({
        "check": "longitudinal_torsion_steel",
        "al_required": round(al_required, 2),
        "min_al": round(min_al, 2),
        "bar_area": round(bar_area, 2),
        "status": "OK" if bar_ok else "FAIL",
    })

    all_ok = all(
        c.get("status") in ("OK", "MAY NEGLECT", None)
        for c in checks
    )

    return {
        "checks": checks,
        "torsion_threshold": round(phi_tth, 4),
        "cracking_moment": round(phi_tcr, 4),
        "at_per_s": round(at_per_s, 6),
        "smax": round(smax, 2),
        "al_required": round(al_required, 2),
        "status": "PASS" if all_ok else "FAIL",
        "geometry": {
            "x1": round(x1, 2), "y1": round(y1, 2),
            "aoh": round(aoh, 2), "ph": round(ph, 2),
            "acp": round(acp, 2), "pcp": round(pcp, 2),
        },
    }
