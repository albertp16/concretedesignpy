# MIT License
# Copyright (c) Albert Pamonag Engineering Consultancy

"""
Joint Shear Verification
=========================

Verifies joint shear requirements for special moment frames
per NSCP 2015 Section 422.7.

Reference: NSCP 2015, ACI 318-19
"""

import math


def _tension_force(as_bar, n_bars, fy, factor=1.25):
    """Probable tensile force in rebar group (kN)."""
    return factor * n_bars * as_bar * fy / 1000.0


def joint_shear_check(
    ve, as1, n_bars1, as2, n_bars2, fy, fc,
    beam_width, joint_depth, perpendicular_dist=0,
    joint_config=1, lamda=1.0, phi=0.85,
):
    """
    Joint shear verification for special moment frames.

    Parameters
    ----------
    ve : float
        Column shear (kN).
    as1 : float
        Area of one bar in group 1 (mm^2).
    n_bars1 : int
        Number of bars in group 1.
    as2 : float
        Area of one bar in group 2 (mm^2).
    n_bars2 : int
        Number of bars in group 2.
    fy : float
        Steel yield strength (MPa).
    fc : float
        Concrete compressive strength (MPa).
    beam_width : float
        Width of the beam (mm).
    joint_depth : float
        Depth of the joint / column dimension (mm).
    perpendicular_dist : float
        Perpendicular distance from column face to beam edge (mm).
    joint_config : int
        1 = confined on all 4 faces (factor 1.7),
        2 = confined on 3 or 2 opposite faces (factor 1.2),
        3 = other (factor 1.0).
    lamda : float
        Lightweight concrete factor (1.0 for normal weight).
    phi : float
        Strength reduction factor for shear.

    Returns
    -------
    dict
        Keys: t1, t2, v_joint, joint_width, aj, vn, phi_vn, status
    """
    t1 = _tension_force(as1, n_bars1, fy)
    t2 = _tension_force(as2, n_bars2, fy)
    v_joint = t1 + t2 - ve

    # Effective joint area
    option_a = beam_width + joint_depth
    option_b = beam_width + 2 * perpendicular_dist
    joint_width = min(option_a, option_b)
    aj = joint_width * joint_depth

    # Joint shear strength factor
    config_factors = {1: 1.7, 2: 1.2, 3: 1.0}
    factor = config_factors.get(joint_config, 1.0)

    vn = factor * lamda * math.sqrt(fc) * aj / 1000.0
    phi_vn = phi * vn

    compliant = phi_vn >= v_joint

    return {
        "t1": round(t1, 2),
        "t2": round(t2, 2),
        "v_joint": round(v_joint, 2),
        "joint_width": round(joint_width, 2),
        "aj": round(aj, 2),
        "vn": round(vn, 2),
        "phi_vn": round(phi_vn, 2),
        "ve": ve,
        "joint_config": joint_config,
        "status": "OK" if compliant else "REINFORCEMENT REQUIRED",
    }
