# MIT License
# Copyright (c) Albert Pamonag Engineering Consultancy

"""
Rebar Layout Generator
========================

Generates rebar coordinate layout for rectangular cross-sections.
Computes stirrup center-line corners and main bar positions.
"""


def section_generate_rect(width, height, cover, ds, db, nx, ny):
    """
    Generate rebar layout for a rectangular cross-section.

    Parameters
    ----------
    width : float
        Section width (mm).
    height : float
        Section height (mm).
    cover : float
        Clear concrete cover (mm).
    ds : float
        Stirrup diameter (mm).
    db : float
        Main bar diameter (mm).
    nx : int
        Number of bars along the width (top and bottom).
    ny : int
        Number of bars along the height (left and right).

    Returns
    -------
    dict
        Keys: section_corners, stirrup_corners, rebar_coords,
              total_bars, rebar_offset
    """
    if width <= 0 or height <= 0:
        raise ValueError("width and height must be positive.")
    if nx < 1 or ny < 1:
        raise ValueError("nx and ny must be at least 1.")

    # Outer rectangle corners
    section_corners = [
        (0.0, 0.0),
        (width, 0.0),
        (width, height),
        (0.0, height),
    ]

    # Stirrup center-line corners
    stirrup_offset = cover + 0.5 * ds
    stirrup_corners = [
        (stirrup_offset, stirrup_offset),
        (width - stirrup_offset, stirrup_offset),
        (width - stirrup_offset, height - stirrup_offset),
        (stirrup_offset, height - stirrup_offset),
    ]

    # Main bar offset from outer face
    rebar_offset = cover + ds + 0.5 * db
    x_left = rebar_offset
    x_right = width - rebar_offset
    y_bot = rebar_offset
    y_top = height - rebar_offset

    rebar_coords = []

    # Horizontal bars (top & bottom)
    if nx == 1:
        x_positions = [(x_left + x_right) / 2.0]
    else:
        dx = (x_right - x_left) / (nx - 1)
        x_positions = [x_left + i * dx for i in range(nx)]

    for x_pos in x_positions:
        rebar_coords.append((round(x_pos, 2), round(y_bot, 2)))
        rebar_coords.append((round(x_pos, 2), round(y_top, 2)))

    # Vertical bars (left & right)
    if ny == 1:
        y_positions = [(y_bot + y_top) / 2.0]
    else:
        dy = (y_top - y_bot) / (ny - 1)
        y_positions = [y_bot + i * dy for i in range(ny)]

    for y_pos in y_positions:
        rebar_coords.append((round(x_left, 2), round(y_pos, 2)))
        rebar_coords.append((round(x_right, 2), round(y_pos, 2)))

    # Remove duplicates and sort
    rebar_coords = list(set(rebar_coords))
    rebar_coords.sort(key=lambda pt: (pt[1], pt[0]))

    return {
        "section_corners": section_corners,
        "stirrup_corners": stirrup_corners,
        "rebar_coords": rebar_coords,
        "total_bars": len(rebar_coords),
        "rebar_offset": round(rebar_offset, 2),
        "x_left": round(x_left, 2),
        "x_right": round(x_right, 2),
        "y_bot": round(y_bot, 2),
        "y_top": round(y_top, 2),
    }
