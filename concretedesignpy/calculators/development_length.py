# MIT License
# Copyright (c) Albert Pamonag Engineering Consultancy

"""
Development Length and Hook Geometry
=====================================

Standard hook geometry for reinforcement detailing
per NSCP 2015 Section 425.

Reference: NSCP 2015 Table 425.3.1 and 425.3.2
"""


def bend_and_hook_deformed_bars(bar_size, angle):
    """
    Standard hook geometry for deformed bars in tension.

    NSCP 2015 Table 425.3.1

    Parameters
    ----------
    bar_size : float
        Bar diameter (mm). Valid: 10-25, 28-36, 40-50.
    angle : int
        Hook angle (degrees). Valid: 90, 180.

    Returns
    -------
    dict
        Keys: bend_diameter (mm), lext (mm)
    """
    if 10 <= bar_size <= 25:
        bend = 6 * bar_size
    elif 28 <= bar_size <= 36:
        bend = 8 * bar_size
    elif 40 <= bar_size <= 50:
        bend = 10 * bar_size
    else:
        raise ValueError(f"Invalid bar size {bar_size} mm. Must be 10-50.")

    if angle == 90:
        lext = 12 * bar_size
    elif angle == 180:
        lext = max(4 * bar_size, 65)
    else:
        raise ValueError(f"Invalid hook angle {angle}. Must be 90 or 180.")

    return {
        "bend_diameter": bend,
        "lext": lext,
        "ldh": lext + bend,
        "bar_size": bar_size,
        "angle": angle,
    }


def bend_and_hook_stirrups(bar_size, angle):
    """
    Hook geometry for stirrups, ties, and hoops.

    NSCP 2015 Table 425.3.2

    Parameters
    ----------
    bar_size : float
        Bar diameter (mm). Valid: 10-16, 20-25.
    angle : int
        Hook angle (degrees). Valid: 90, 135, 180.

    Returns
    -------
    dict
        Keys: bend_diameter (mm), lext (mm)
    """
    if 10 <= bar_size <= 16:
        bend = 4 * bar_size
    elif 20 <= bar_size <= 25:
        bend = 6 * bar_size
    else:
        raise ValueError(f"Invalid bar size {bar_size} mm for stirrups. Must be 10-25.")

    if angle == 90:
        if bar_size <= 16:
            lext = max(6 * bar_size, 75)
        else:
            lext = 12 * bar_size
    elif angle == 135:
        lext = max(6 * bar_size, 75)
    elif angle == 180:
        lext = max(4 * bar_size, 65)
    else:
        raise ValueError(f"Invalid hook angle {angle}. Must be 90, 135, or 180.")

    return {
        "bend_diameter": bend,
        "lext": lext,
        "ldh": lext + bend,
        "bar_size": bar_size,
        "angle": angle,
    }
