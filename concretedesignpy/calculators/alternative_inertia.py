# MIT License
# Copyright (c) Albert Pamonag Engineering Consultancy

"""
Alternative Moment of Inertia
===============================

Computes alternative (effective) moment of inertia
per NSCP 2015 Section 424.2.3.5.
"""


def alternative_inertia_column_wall(ig, ast, ag, mu, pu, h, po):
    """
    Alternative moment of inertia for columns and walls.

    NSCP 2015 Section 424.2.3.5

    Parameters
    ----------
    ig : float
        Gross moment of inertia (mm^4).
    ast : float
        Total area of longitudinal steel (mm^2).
    ag : float
        Gross cross-sectional area (mm^2).
    mu : float
        Factored moment (kN-m or consistent units).
    pu : float
        Factored axial load (kN or consistent units).
    h : float
        Section depth in direction of bending (mm).
    po : float
        Nominal axial strength at zero eccentricity (kN or consistent units).

    Returns
    -------
    dict
        Keys: ialt, ireq, imin, imax
    """
    if ag <= 0 or po <= 0 or pu <= 0:
        raise ValueError("ag, po, pu must be positive.")

    imin = 0.35 * ig
    ireq = (0.80 + 25.0 * (ast / ag)) * (1 - mu / (pu * h) - 0.5 * pu / po) * ig
    imax = 0.875 * ig

    if ireq < imin:
        ialt = imin
    elif ireq > imax:
        ialt = imax
    else:
        ialt = ireq

    return {
        "ialt": round(ialt, 2),
        "ireq": round(ireq, 2),
        "imin": round(imin, 2),
        "imax": round(imax, 2),
    }


def alternative_inertia_beam_slab(ig, rho, b, d):
    """
    Alternative moment of inertia for beams, flat plates, and flat slabs.

    NSCP 2015 Section 424.2.3.5

    Parameters
    ----------
    ig : float
        Gross moment of inertia (mm^4).
    rho : float
        Longitudinal reinforcement ratio (As / (b * d)).
    b : float
        Effective width (mm).
    d : float
        Effective depth (mm).

    Returns
    -------
    dict
        Keys: ialt, ireq, imin, imax
    """
    if d <= 0:
        raise ValueError("d must be positive.")

    imin = 0.25 * ig
    ireq = (0.10 + 25.0 * rho) * (1.2 - 0.2 * (b / d)) * ig
    imax = 0.5 * ig

    if ireq < imin:
        ialt = imin
    elif ireq > imax:
        ialt = imax
    else:
        ialt = ireq

    return {
        "ialt": round(ialt, 2),
        "ireq": round(ireq, 2),
        "imin": round(imin, 2),
        "imax": round(imax, 2),
    }
