# SPDX-License-Identifier: MIT
"""
area_calculations.py
--------------------
A collection of utility functions for calculating:
- The area of a circle given its diameter,
- The total steel area given the number of bars and the area of each bar,
- The ratio of steel area to concrete area.

This code is open-source and licensed under the MIT License.
"""

import math

def area_diam(diameter: float) -> float:
    """
    Return the cross-sectional area of a circle given its diameter.

    Args:
        diameter (float): The diameter of the circle (must be positive).

    Returns:
        float: The area of the circle.

    Raises:
        ValueError: If 'diameter' is not positive.
    """
    if diameter <= 0:
        raise ValueError("Diameter must be positive.")
    return (math.pi / 4) * (diameter ** 2)


def steel_area(num_bars: int, bar_area: float) -> float:
    """
    Return total steel area for a given number of bars.

    Args:
        num_bars (int): Number of bars (must be positive).
        bar_area (float): Cross-sectional area of one bar (must be positive).

    Returns:
        float: The total steel area.

    Raises:
        ValueError: If 'num_bars' <= 0 or 'bar_area' <= 0.
    """
    if num_bars <= 0:
        raise ValueError("Number of bars must be positive.")
    if bar_area <= 0:
        raise ValueError("Bar area must be positive.")
    return num_bars * bar_area


def area_ratio(steel_area_val: float, concrete_area_val: float) -> float:
    """
    Return the ratio of steel area to concrete area.

    Args:
        steel_area_val (float): Total steel area (must be positive).
        concrete_area_val (float): Concrete area (must be positive).

    Returns:
        float: The ratio of steel to concrete area.

    Raises:
        ValueError: If 'steel_area_val' <= 0 or 'concrete_area_val' <= 0.
    """
    if steel_area_val <= 0:
        raise ValueError("Steel area must be positive.")
    if concrete_area_val <= 0:
        raise ValueError("Concrete area must be positive.")
    return steel_area_val / concrete_area_val
