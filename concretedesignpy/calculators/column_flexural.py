# MIT License
# Copyright (c) Albert Pamonag Engineering Consultancy

"""
Column Minimum Flexural Strength Check
========================================

Verifies minimum flexural strength of columns
per NSCP 2015 Section 418.7.3.
"""


def check_min_flexural_strength(mn_col_top, mn_col_bot, mn_beam, limit=1.2):
    """
    Check minimum flexural strength ratio of columns to beams.

    Per NSCP 2015 Section 418.7.3, the sum of column moment capacities
    at a joint must exceed the sum of beam moment capacities by a factor.

    Parameters
    ----------
    mn_col_top : float
        Flexural capacity of the top column (kN-m).
    mn_col_bot : float
        Flexural capacity of the bottom column (kN-m).
    mn_beam : float
        Total flexural capacity of beams at the joint (kN-m).
    limit : float
        Required minimum ratio (default 1.2).

    Returns
    -------
    dict
        Keys: ratio, mn_col_total, mn_beam, limit, compliant, recommendation
    """
    if mn_beam <= 0:
        raise ValueError("mn_beam must be positive.")

    mn_col = mn_col_top + mn_col_bot
    ratio = mn_col / mn_beam
    compliant = ratio >= limit

    recommendation = ""
    if not compliant:
        recommendation = (
            "Increase column reinforcement, reduce beam reinforcement, "
            "or do both to achieve compliance."
        )

    return {
        "ratio": round(ratio, 4),
        "mn_col_total": round(mn_col, 2),
        "mn_beam": round(mn_beam, 2),
        "limit": limit,
        "compliant": compliant,
        "status": "COMPLIANT" if compliant else "NONCOMPLIANT",
        "recommendation": recommendation,
    }
