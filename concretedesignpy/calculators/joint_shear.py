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
    beam_width, joint_depth, *, column_width,
    perpendicular_dist=0,
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
        Width of the beam, b (mm).
    joint_depth : float
        Overall depth of the column, h, in the direction of the joint
        shear considered (mm).
    column_width : float
        Overall width of the column perpendicular to the joint shear
        considered (mm). Keyword-only and required: without it the
        governing cap on Aj cannot be applied.
    perpendicular_dist : float
        Perpendicular distance from the column side face to the nearest
        beam edge, x (mm). Zero when the beam is flush with, or wider
        than, the column.

        NOTE for the next reader: this is measured **beam face to column
        face**, and that is deliberate. Section 15.5.2.2(b) is worded
        "twice the perpendicular distance from longitudinal axis of beam
        to nearest side face of the column", i.e. 2 * (b/2 + x), which is
        identically b + 2x. Fig. R15.5.2.2 draws the same limit as
        "b + 2x" with x face-measured. Do not "fix" this back to 2 * x.
    joint_config : int
        1 = confined on all 4 faces (gamma 1.7),
        2 = confined on 3 or 2 opposite faces (gamma 1.2),
        3 = other (gamma 1.0).
    lamda : float
        Lightweight concrete factor (1.0 for normal weight).
    phi : float
        Strength reduction factor for shear.

    Returns
    -------
    dict
        Keys: t1, t2, v_joint, joint_width, aj, vn, phi_vn, status

    Notes
    -----
    This function checks joint shear only. It does not check joint
    transverse reinforcement, bar development into the joint, or the
    column-to-beam flexural strength ratio.
    """
    if column_width <= 0:
        raise ValueError("column_width must be positive.")

    t1 = _tension_force(as1, n_bars1, fy)
    t2 = _tension_force(as2, n_bars2, fy)
    v_joint = t1 + t2 - ve

    # Effective joint area, Section 15.5.2.2.
    # "Effective joint width shall be the overall width of the column
    #  where the beam is wider than the column. Where the column is wider
    #  than the beam, effective joint width shall not exceed the lesser of
    #  (a) and (b)."  R15.5.2.2: "In no case is Aj greater than the column
    #  cross-sectional area."
    # Both sentences are satisfied by taking the least of the three.
    option_a = beam_width + joint_depth
    option_b = beam_width + 2 * perpendicular_dist
    joint_width = min(column_width, option_a, option_b)
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
        "column_width": column_width,
        "aj": round(aj, 2),
        "vn": round(vn, 2),
        "phi_vn": round(phi_vn, 2),
        "ve": ve,
        "joint_config": joint_config,
        "status": "OK" if compliant else "REINFORCEMENT REQUIRED",
    }
