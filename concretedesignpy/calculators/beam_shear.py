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
