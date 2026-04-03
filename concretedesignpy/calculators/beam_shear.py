# MIT License
# Copyright (c) Albert Pamonag Engineering Consultancy

"""
Beam Shear Capacity Calculator
================================

Computes concrete and steel shear strength, and required stirrup spacing.

Reference: NSCP 2015 Section 422.5 / ACI 318-19
"""

import math


def compute_concrete_shear_strength(
    fc, b, d, lamda=1.0, vc_type="simple",
    vu=None, mu=None, nu=None, rho_w=None,
    h=None, db=None, ds=None, cover=None, n_bars=None,
):
    """
    Compute concrete shear strength Vc.

    Parameters
    ----------
    fc : float
        Concrete compressive strength (MPa).
    b : float
        Beam width (mm).
    d : float
        Effective depth (mm).
    lamda : float
        Lightweight concrete factor (1.0 for normal weight).
    vc_type : str
        One of: 'simple', 'detailed', 'axial_compression',
        'detailed_axial_compression', 'axial_tension'.
    vu, mu, nu : float or None
        Factored shear (N), moment (N-mm), axial force (N).
    rho_w : float or None
        Longitudinal reinforcement ratio. If None, computed from inputs.
    h, db, ds, cover, n_bars : float or None
        Section geometry for computing Ag and rho_w when needed.

    Returns
    -------
    dict
        Keys: vc (N), vc_kn (kN), type
    """
    if fc <= 0 or b <= 0 or d <= 0:
        raise ValueError("fc, b, d must all be positive.")

    # Compute Ag and rho_w if not given
    if h is None:
        h = d + (db / 2 if db else 0) + (ds if ds else 0) + (cover if cover else 40)
    ag = h * b

    if rho_w is None and n_bars and db:
        as_bar = n_bars * math.pi * (db / 2) ** 2
        rho_w = as_bar / (b * d)

    if vc_type == "simple":
        vc = (1.0 / 6.0) * lamda * math.sqrt(fc) * b * d

    elif vc_type == "detailed":
        if vu is None or mu is None:
            raise ValueError("vu and mu are required for 'detailed' type.")
        if rho_w is None:
            raise ValueError("rho_w is required for 'detailed' type.")
        ratio = min(vu * d / mu, 1.0) if mu != 0 else 1.0
        vc1 = (0.16 * lamda * math.sqrt(fc) + 17 * rho_w * ratio) * b * d
        vc2 = 0.29 * lamda * math.sqrt(fc) * b * d
        vc = min(vc1, vc2)

    elif vc_type == "axial_compression":
        if nu is None:
            raise ValueError("nu is required for 'axial_compression' type.")
        vc = (1.0 / 6.0) * (1 + nu / (14 * ag)) * lamda * math.sqrt(fc) * b * d

    elif vc_type == "detailed_axial_compression":
        if nu is None or vu is None or mu is None:
            raise ValueError("nu, vu, mu required for 'detailed_axial_compression'.")
        if rho_w is None:
            raise ValueError("rho_w is required.")
        mm = mu - nu * ((4 * h - d) / 8)
        vc2 = 0.29 * lamda * math.sqrt(fc) * b * d * math.sqrt(1 + 0.29 * nu / ag)
        if mm <= 0:
            vc = vc2
        else:
            vc1 = (0.16 * lamda * math.sqrt(fc) + 17 * rho_w * (vu * d / mm)) * b * d
            vc = min(vc1, vc2)

    elif vc_type == "axial_tension":
        if nu is None:
            raise ValueError("nu is required for 'axial_tension' type.")
        vc = (1.0 / 6.0) * (1 + 0.29 * nu / ag) * lamda * math.sqrt(fc) * b * d
        vc = max(vc, 0)

    else:
        raise ValueError(f"Unknown vc_type: {vc_type}")

    return {"vc": round(vc, 2), "vc_kn": round(vc / 1000, 2), "type": vc_type}


def compute_steel_shear_strength(av, fyt, d, s):
    """
    Compute steel shear reinforcement strength Vs = Av * fyt * d / s.

    Parameters
    ----------
    av : float
        Total area of shear reinforcement (mm^2).
    fyt : float
        Yield strength of transverse reinforcement (MPa).
    d : float
        Effective depth (mm).
    s : float
        Stirrup spacing (mm).

    Returns
    -------
    dict
        Keys: vs (N), vs_kn (kN)
    """
    if s <= 0:
        raise ValueError("Stirrup spacing must be positive.")
    vs = av * fyt * d / s
    return {"vs": round(vs, 2), "vs_kn": round(vs / 1000, 2)}


def compute_shear_spacing(fc, b, d, fyt, vu_required, phi, av,
                          lamda=1.0, vc_type="simple", **kwargs):
    """
    Compute required stirrup spacing given a required Vu.

    Parameters
    ----------
    fc, b, d, fyt : float
        Material and section properties.
    vu_required : float
        Required ultimate shear (N).
    phi : float
        Strength reduction factor for shear.
    av : float
        Total area of shear reinforcement (mm^2).
    lamda : float
        Lightweight concrete factor.
    vc_type : str
        Type of Vc calculation.

    Returns
    -------
    dict
        Keys: spacing (mm), vs_required (N), smax (mm)
    """
    vc_result = compute_concrete_shear_strength(fc, b, d, lamda, vc_type, **kwargs)
    vc = vc_result["vc"]
    vs_required = (vu_required / phi) - vc
    vs_limit = (1.0 / 3.0) * math.sqrt(fc) * b * d

    if vu_required > phi * vc:
        use_s = av * fyt * d / vs_required if vs_required > 0 else 600.0
    else:
        avmin1 = 0.062 * math.sqrt(fc) * b / fyt
        avmin2 = 0.35 * b / fyt
        avmin_per_s = max(avmin1, avmin2)
        use_s = av / avmin_per_s

    if vs_required <= vs_limit:
        smax = min(d / 2.0, 600.0)
    else:
        smax = min(d / 4.0, 300.0)

    spacing = min(use_s, smax)

    return {
        "spacing": round(spacing, 2),
        "vs_required": round(vs_required, 2),
        "vs_required_kn": round(vs_required / 1000, 2),
        "smax": round(smax, 2),
        "vc_kn": vc_result["vc_kn"],
    }


def shear_torsion_design(fc, fyv, fy, phi, bw, h, cc, c, d,
                         vu, tu, nu, s_chosen, n_legs, db_stirrup, db_long):
    """
    Combined shear and torsion design per ACI 318M-14.

    Follows the computation flow from EA Spreadsheet Suite.

    Parameters
    ----------
    fc : float       Concrete compressive strength (MPa).
    fyv : float      Yield strength of stirrups (MPa).
    fy : float       Yield strength of longitudinal rft (MPa).
    phi : float      Strength reduction factor for shear/torsion.
    bw : float       Width of beam (mm).
    h : float        Total height of beam (mm).
    cc : float       Clear cover to links (mm).
    c : float        Cover to centroid of reinforcement (mm).
    d : float        Effective depth (mm).
    vu : float       Factored shear force (kN).
    tu : float       Factored torsional moment (kN.m).
    nu : float       Factored axial force (kN), positive = compression.
    s_chosen : float Chosen stirrup spacing (mm).
    n_legs : int     Total number of stirrup legs.
    db_stirrup : float  Stirrup bar diameter (mm).
    db_long : float  Longitudinal bar diameter for torsion (mm).

    Returns
    -------
    dict with all intermediate values and step references.
    """
    steps = []
    ag = bw * h  # gross area mm^2

    # ── 1. Material checks ──
    fc_check = "OK" if fc <= 70 else "fc' > 70 MPa limit"
    fyv_check = "OK" if fyv <= 420 else "fyv > 420 MPa limit"
    steps.append({
        "title": "Material Properties",
        "items": [
            {"label": "f'c", "value": fc, "unit": "MPa", "check": fc_check,
             "ref": "ACI 318M-14 Cl. 22.5.3.1"},
            {"label": "fyv", "value": fyv, "unit": "MPa", "check": fyv_check,
             "ref": "ACI 318M-14 Cl. 20.2.2.4"},
            {"label": "fy", "value": fy, "unit": "MPa"},
            {"label": "\u03c6", "value": phi, "ref": "ACI 318M-14 Cl. 21.2.1(b)"},
        ],
    })

    # ── 2. Section dimensions ──
    steps.append({
        "title": "Section Dimensions",
        "items": [
            {"label": "bw", "value": bw, "unit": "mm"},
            {"label": "h", "value": h, "unit": "mm"},
            {"label": "Cc", "value": cc, "unit": "mm"},
            {"label": "c", "value": c, "unit": "mm"},
            {"label": "d", "value": d, "unit": "mm",
             "formula": "h - c = {} - {} = {}".format(h, c, d)},
        ],
    })

    # ── 3. Concrete shear strength Vc ──
    # ACI 318M-14 Cl. 22.5.5.1 with axial load modification
    if nu > 0:
        # Compression: Cl. 22.5.6.1
        vc_1 = (1 + nu * 1000 / (14 * ag)) * (math.sqrt(fc) / 6) * bw * d / 1000
        vc_2 = 0.3 * math.sqrt(fc) * bw * d * (1 + 0.3 * nu * 1000 / ag) / 1000
        vc = min(vc_1, vc_2)
        vc_note = "with axial compression"
    elif nu < 0:
        # Tension: Cl. 22.5.7.1
        vc = max(0, (1 + 0.29 * nu * 1000 / ag) * (math.sqrt(fc) / 6) * bw * d / 1000)
        vc_note = "with axial tension"
    else:
        vc = (math.sqrt(fc) / 6) * bw * d / 1000
        vc_note = "no axial load"

    steps.append({
        "title": "Concrete Shear Strength",
        "items": [
            {"label": "Vc", "value": round(vc, 2), "unit": "kN",
             "note": vc_note,
             "ref": "ACI 318M-14 Cl. 22.5.5.1\u20137"},
        ],
    })

    # ── 4. Required stirrup shear strength Vs ──
    vs_max = (2.0 / 3.0) * math.sqrt(fc) * bw * d / 1000
    vs_req_raw = (vu - phi * vc) / phi

    if vs_req_raw <= 0:
        vs = 0
        shear_status = "Provide Min. Rft"
    elif vs_req_raw > vs_max:
        vs = vs_req_raw
        shear_status = "UNSAFE - Vs exceeds limit"
    else:
        vs = vs_req_raw
        shear_status = "SAFE"

    # Check: Vu <= phi*(Vc + 0.66*sqrt(fc')*bw*d)
    vu_max = phi * (vc + (2.0 / 3.0) * math.sqrt(fc) * bw * d / 1000)
    overall_check = "SAFE" if vu <= vu_max else "UNSAFE"

    # Av/s required for shear
    av_s_shear = (vu * 1000 - phi * vc * 1000) / (phi * fyv * d) if vs > 0 else 0

    steps.append({
        "title": "Transverse Rft for Shear",
        "items": [
            {"label": "Vs,req", "value": round(max(vs, 0), 2), "unit": "kN",
             "formula": "(Vu - \u03c6Vc)/\u03c6 = ({} - {}x{})/{} = {}".format(
                 vu, phi, round(vc, 2), phi, round(vs_req_raw, 2)),
             "ref": "ACI 318M-14 Cl. 22.5.10.5.3"},
            {"label": "Vs,max", "value": round(vs_max, 2), "unit": "kN",
             "formula": "\u03c6(2/3)\u221afc'bwd",
             "ref": "ACI 318M-14 Cl. 22.5.1.2"},
            {"label": "Check", "value": overall_check,
             "formula": "Vu \u2264 \u03c6Vc + \u03c6(0.66\u221afc')bwd = {} kN".format(
                 round(vu_max, 2)),
             "ref": "ACI 318M-14 Cl. 22.5.1.2"},
            {"label": "Av/s (shear)", "value": round(av_s_shear, 4), "unit": "mm\u00b2/mm",
             "ref": "ACI 318M-14 R22.5.10.5"},
        ],
    })

    # ── 5. Torsion geometry ──
    # Assume stirrup diameter for Aoh calculation
    aoh = (h - 2 * cc - db_stirrup) * (bw - 2 * cc - db_stirrup)
    ph = 2 * (h - 2 * cc - db_stirrup) + 2 * (bw - 2 * cc - db_stirrup)

    steps.append({
        "title": "Torsion Geometry",
        "items": [
            {"label": "Aoh", "value": round(aoh, 0), "unit": "mm\u00b2",
             "ref": "ACI 318M-14 R22.7.6.1.1"},
            {"label": "Ph", "value": round(ph, 0), "unit": "mm"},
        ],
    })

    # ── 6. Torsional thresholds ──
    # fTcr = phi * sqrt(fc') * (bw*h)^2 / (2*(bw+h)) / 3  (in kN.m)
    tcr = phi * math.sqrt(fc) * (bw * h) ** 2 / (2 * (bw + h)) / 3 / 1e6
    # fTth = fTcr / 4
    tth = tcr / 4

    if tu < tth:
        torsion_action = "Neglect Torsion"
    else:
        torsion_action = "Design for Torsion"

    # Cross-section dimension check
    shear_stress = vu * 1000 / (bw * d)
    if tu >= tth:
        torsion_stress = tu * 1e6 * ph / (1.7 * aoh ** 2)
        combined = math.sqrt(shear_stress ** 2 + torsion_stress ** 2)
        limit = phi * (vc * 1000 / (bw * d) + (2.0 / 3.0) * math.sqrt(fc))
        dim_check = "SAFE" if combined <= limit else "UNSAFE"
    else:
        torsion_stress = 0
        combined = shear_stress
        limit = phi * (vc * 1000 / (bw * d) + (2.0 / 3.0) * math.sqrt(fc))
        dim_check = "SAFE"

    steps.append({
        "title": "Torsion Thresholds",
        "items": [
            {"label": "\u03c6Tcr", "value": round(tcr, 4), "unit": "kN.m",
             "ref": "ACI 318M-14 Cl. 22.7.5.1"},
            {"label": "\u03c6Tth", "value": round(tth, 4), "unit": "kN.m",
             "ref": "ACI 318M-14 Cl. 22.7.4.1"},
            {"label": "Action", "value": torsion_action},
            {"label": "Section Check", "value": dim_check,
             "ref": "ACI 318M-14 Cl. 22.7.7.1"},
        ],
    })

    # ── 7. Torsion reinforcement ──
    if torsion_action == "Neglect Torsion":
        at_s = 0
        al = 0
        al_min = 0
    else:
        # At/s = Tu / (phi * 2 * 1.0 * 0.85 * Aoh * fyv)  per leg
        at_s = tu * 1e6 / (phi * 2 * 0.85 * aoh * fyv)
        # Al = Tu * Ph / (2 * phi * Aoh * fy)
        al = tu * 1e6 * ph / (2 * phi * aoh * fy)
        # Al,min
        al_min = (5 * math.sqrt(fc) * bw * h / (12 * fy)
                  - max(at_s, 0.175 * bw / fyv) * ph * fyv / fy)

    steps.append({
        "title": "Torsion Reinforcement",
        "items": [
            {"label": "At/s", "value": round(at_s, 4), "unit": "mm\u00b2/mm/leg",
             "ref": "ACI 318M-14 Cl. 22.7.6.1a"},
            {"label": "A\u2113", "value": round(max(al, al_min, 0), 2), "unit": "mm\u00b2",
             "ref": "ACI 318M-14 Cl. 22.7.6.1b"},
            {"label": "A\u2113,min", "value": round(max(al_min, 0), 2), "unit": "mm\u00b2",
             "ref": "ACI 318M-14 Cl. 9.6.4.3"},
        ],
    })

    # ── 8. Minimum shear reinforcement ──
    av_min_1 = (1.0 / 16.0) * math.sqrt(fc) * bw / fyv
    av_min_2 = 0.35 * bw / fyv
    av_min = max(av_min_1, av_min_2)

    steps.append({
        "title": "Minimum Shear Reinforcement",
        "items": [
            {"label": "Av,min/s", "value": round(av_min, 4), "unit": "mm\u00b2/mm",
             "ref": "ACI 318M-14 Table 9.6.3.3"},
        ],
    })

    # ── 9. Total required ──
    total_req = max(av_s_shear + 2 * at_s, av_min)

    steps.append({
        "title": "Total Required Reinforcement",
        "items": [
            {"label": "(Av+2At)/s", "value": round(total_req, 4),
             "unit": "mm\u00b2/mm",
             "ref": "ACI 318M-14 Cl. 9.6.4.2"},
            {"label": "A\u2113 (torsion)", "value": round(max(al, al_min, 0), 2),
             "unit": "mm\u00b2"},
        ],
    })

    # ── 10. Chosen reinforcement check ──
    av_leg = math.pi * (db_stirrup / 2) ** 2
    av_provided_per_s = n_legs * av_leg / s_chosen

    # Max spacing
    if vs > (1.0 / 3.0) * math.sqrt(fc) * bw * d / 1000:
        smax = min(300, d / 4, ph / 8)
    else:
        smax = min(600, d / 2)
    if torsion_action != "Neglect Torsion":
        smax = min(smax, ph / 8, 300)

    spacing_check = "OK" if s_chosen <= smax else "NOT OK"
    rft_check = "OK" if av_provided_per_s >= total_req else "NOT OK"

    # Compute chosen bar sizes
    if total_req > av_min:
        ext_dia = math.ceil(math.sqrt(4 / math.pi * (at_s * s_chosen + av_s_shear * s_chosen / n_legs)) / 2) * 2
    else:
        ext_dia = math.ceil(math.sqrt(4 / math.pi * av_min * s_chosen / n_legs) / 2) * 2

    steps.append({
        "title": "Chosen Reinforcement",
        "items": [
            {"label": "Spacing", "value": s_chosen, "unit": "mm",
             "check": spacing_check},
            {"label": "Smax", "value": round(smax, 0), "unit": "mm",
             "ref": "ACI 318M-14 Cl. 9.7.6.2.2"},
            {"label": "Legs", "value": n_legs},
            {"label": "Provided (Av+2At)/s", "value": round(av_provided_per_s, 4),
             "unit": "mm\u00b2/mm", "check": rft_check},
            {"label": "Suggested ext. stirrup \u2300", "value": ext_dia, "unit": "mm"},
        ],
    })

    return {
        "steps": steps,
        "vc": round(vc, 2),
        "vs": round(max(vs, 0), 2),
        "vs_max": round(vs_max, 2),
        "vu_max": round(vu_max, 2),
        "av_s_shear": round(av_s_shear, 4),
        "at_s": round(at_s, 4),
        "al": round(max(al, al_min, 0), 2),
        "total_req": round(total_req, 4),
        "av_min": round(av_min, 4),
        "aoh": round(aoh, 0),
        "ph": round(ph, 0),
        "tcr": round(tcr, 4),
        "tth": round(tth, 4),
        "torsion_action": torsion_action,
        "shear_status": overall_check,
        "dim_check": dim_check,
        "smax": round(smax, 0),
        "spacing_check": spacing_check,
    }


def shear_design(fc, fyv, phi, bw, h, cc, c, d, vu, nu,
                 s_chosen, n_legs, db_stirrup):
    """
    Beam shear design per ACI 318M-14 (shear only, no torsion).

    Parameters
    ----------
    fc : float       Concrete compressive strength (MPa).
    fyv : float      Yield strength of stirrups (MPa).
    phi : float      Strength reduction factor for shear.
    bw : float       Width of beam (mm).
    h : float        Total height of beam (mm).
    cc : float       Clear cover to links (mm).
    c : float        Cover to centroid of reinforcement (mm).
    d : float        Effective depth (mm).
    vu : float       Factored shear force (kN).
    nu : float       Factored axial force (kN), positive = compression.
    s_chosen : float Chosen stirrup spacing (mm).
    n_legs : int     Total number of stirrup legs.
    db_stirrup : float  Stirrup bar diameter (mm).

    Returns
    -------
    dict with all intermediate values for shear design report.
    """
    ag = bw * h  # gross area mm^2

    # ── 1. Concrete shear strength Vc ──
    if nu > 0:
        vc_1 = (1 + nu * 1000 / (14 * ag)) * (math.sqrt(fc) / 6) * bw * d / 1000
        vc_2 = 0.3 * math.sqrt(fc) * bw * d * (1 + 0.3 * nu * 1000 / ag) / 1000
        vc = min(vc_1, vc_2)
        vc_note = "with axial compression"
    elif nu < 0:
        vc = max(0, (1 + 0.29 * nu * 1000 / ag) * (math.sqrt(fc) / 6) * bw * d / 1000)
        vc_note = "with axial tension"
    else:
        vc = (math.sqrt(fc) / 6) * bw * d / 1000
        vc_note = "no axial load"

    # ── 2. Required stirrup shear strength Vs ──
    vs_max = (2.0 / 3.0) * math.sqrt(fc) * bw * d / 1000
    vs_req = (vu - phi * vc) / phi

    if vs_req <= 0:
        vs = 0
    elif vs_req > vs_max:
        vs = vs_req
    else:
        vs = vs_req

    vu_max = phi * (vc + vs_max)
    shear_status = "SAFE" if vu <= vu_max else "UNSAFE"

    # Av/s required for shear
    av_s_req = (vu * 1000 - phi * vc * 1000) / (phi * fyv * d) if vs > 0 else 0

    # ── 3. Minimum shear reinforcement ──
    av_min_1 = (1.0 / 16.0) * math.sqrt(fc) * bw / fyv
    av_min_2 = 0.35 * bw / fyv
    av_min = max(av_min_1, av_min_2)

    # Governing Av/s
    av_s_govern = max(av_s_req, av_min)

    # ── 4. Maximum spacing ──
    if vs > (1.0 / 3.0) * math.sqrt(fc) * bw * d / 1000:
        smax = min(300, d / 4)
    else:
        smax = min(600, d / 2)

    # ── 5. Chosen reinforcement check ──
    av_leg = math.pi * (db_stirrup / 2) ** 2
    av_provided = n_legs * av_leg / s_chosen

    spacing_ok = "OK" if s_chosen <= smax else "NOT OK"
    rft_ok = "OK" if av_provided >= av_s_govern else "NOT OK"

    return {
        "vc": round(vc, 2),
        "vc_note": vc_note,
        "vs": round(max(vs, 0), 2),
        "vs_req": round(max(vs_req, 0), 2),
        "vs_max": round(vs_max, 2),
        "vu_max": round(vu_max, 2),
        "shear_status": shear_status,
        "av_s_req": round(max(av_s_req, 0), 4),
        "av_min": round(av_min, 4),
        "av_s_govern": round(av_s_govern, 4),
        "smax": round(smax, 0),
        "av_provided": round(av_provided, 4),
        "spacing_ok": spacing_ok,
        "rft_ok": rft_ok,
    }
