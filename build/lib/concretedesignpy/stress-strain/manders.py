#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MIT License

Copyright (c) 2025 Albert Pamonag

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights 
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell 
copies of the Software, and to permit persons to whom the Software is 
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in 
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, 
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN 
THE SOFTWARE.

-------------------------------------------------------------------------------
This module provides functions to compute the confined concrete strength and 
its associated factors following (primarily) Mander's model.
"""

import numpy as np
import matplotlib.pyplot as plt


def compute_confined_strength(fpco: float, fpl: float) -> float:
    """Computes confined concrete strength using Mander's model.

    This function calculates the maximum confined compressive strength of 
    concrete (f'cc) given the unconfined compressive strength (f'co) and the 
    effective lateral confining stress (fpl).

    Args:
        fpco (float): Unconfined concrete strength, f'co, in MPa.
        fpl (float): Effective lateral confining stress, fpl, in MPa.

    Returns:
        float: Confined concrete strength, f'cc, in MPa.
    """
    A = 2.254 * np.sqrt(1.0 + ((7.92 * fpl) / fpco))
    B = 2.0 * (fpl / fpco)
    fpcc = fpco * (-1.254 + A - B)
    return fpcc


def compute_confined_strength_ratio(ratio: float) -> float:
    """Computes the ratio of confined strength to unconfined strength.

    This ratio is given by (f'cc / f'co) as a function of the lateral
    confinement stress ratio (fpl / f'co).

    Args:
        ratio (float): Lateral confinement stress ratio (fpl / f'co).

    Returns:
        float: Confined strength ratio (f'cc / f'co).
    """
    A = 2.254 * np.sqrt(1.0 + 7.92 * ratio)
    B = 2.0 * ratio
    return -1.254 + A - B


def compute_equation_A(r: float) -> float:
    """Computes coefficient A based on a confinement stress ratio.

    The expression is drawn from a modified approach to Mander's model, which 
    accounts for varying lateral confinement stress ratios (fl1 / fl2).

    Args:
        r (float): The ratio fl1 / fl2 of lateral stresses.

    Returns:
        float: Coefficient A for subsequent use in calculating the
            confined strength factor.
    """
    return 6.8886 - (0.6069 + (17.275 * r)) * np.exp(-4.989 * r)


def compute_equation_B(A: float, r: float) -> float:
    """Computes coefficient B based on coefficient A and confinement ratio.

    Args:
        A (float): Coefficient A from `compute_equation_A`.
        r (float): The ratio fl1 / fl2 of lateral stresses.

    Returns:
        float: Coefficient B for subsequent use in calculating the
            confined strength factor.
    """
    B2 = ((5.0 / A) * (0.9849 - 0.6306 * np.exp(-3.8939 * r))) - 0.1
    return (4.5 / B2) - 5.0


def compute_confined_strength_factor(A: float, B: float, bar_x: float) -> float:
    """Computes the confined strength factor K.

    K is a factor that modifies the unconfined strength f'co to reflect 
    additional confinement. Typically K = (f'cc / f'co).

    Args:
        A (float): Coefficient A from `compute_equation_A`.
        B (float): Coefficient B from `compute_equation_B`.
        bar_x (float): Normalized confining stress (fl / f'co).

    Returns:
        float: Confined strength factor K = (f'cc / f'co).
    """
    return 1.0 + A * bar_x * (0.1 + (0.9 / (1.0 + B * bar_x)))


def solver(r: float, xi: float) -> float:
    """Computes the confined strength factor K for a given ratio r and xi.

    This function combines the above intermediate computations to solve for 
    the confined strength factor.

    Args:
        r (float): Confinement stress ratio (fl1 / fl2).
        xi (float): Normalized confinement stress (fl / f'co).

    Returns:
        float: The confined strength factor K, i.e., f'cc / f'co.
    """
    A_val = compute_equation_A(r)
    B_val = compute_equation_B(A_val, r)
    k1 = A_val * (0.1 + 0.9 / (1.0 + B_val * xi))
    return 1.0 + (k1 * xi)


def plot_confined_strength() -> None:
    """Generates a plot of confined strength ratio vs. confinement stress ratio.

    This function compares various scenarios by varying:
      - fl1/fl2 from 0 to 0.2
      - f'l2/f'co from 0 to 0.3

    It also draws a reference curve where fl1 = fl2 (i.e., the original Mander 
    model ratio).
    """
    ratio_values = np.arange(0, 0.31, 0.01)

    # Compute solver-based data
    value0 = [solver(0.0, r) for r in ratio_values]
    value1 = [solver(0.1, r) for r in ratio_values]
    value2 = [solver(0.2, r) for r in ratio_values]

    # Compute the original Mander's approach for fl1=fl2
    value_fl1_eq_fl2 = [compute_confined_strength_ratio(r) for r in ratio_values]

    plt.figure(figsize=(8, 6))
    plt.plot(value0, -ratio_values, label="fl1/fl2 = 0.0")
    plt.plot(value1, -ratio_values, label="fl1/fl2 = 0.1")
    plt.plot(value2, -ratio_values, label="fl1/fl2 = 0.2")
    plt.plot(value_fl1_eq_fl2, -ratio_values, label="fl1 = fl2", linestyle="dashed")

    plt.xlabel(r"Confined Strength Ratio, $f'_{cc}/f'_{co}$")
    plt.ylabel(r"Largest Confining Stress Ratio, $f'_{l2}/f'_{co}$")
    plt.title("Confined Strength Determination from Lateral Confining Stresses")
    plt.legend()
    plt.grid(True)
    plt.show()


def main() -> None:
    """Main entry point for executing the script.

    Demonstrates how to use the functions defined above.
    """
    # Example usage for a single pair of fpco and fpl
    fpco_example = 30.0  # MPa, unconfined strength
    fpl_example = 3.0    # MPa, lateral confinement

    fcc = compute_confined_strength(fpco_example, fpl_example)
    print(f"Unconfined strength (f'co) = {fpco_example} MPa")
    print(f"Effective lateral stress (fpl) = {fpl_example} MPa")
    print(f"Confined strength (f'cc) = {fcc:.2f} MPa")

    # Plot the relationship between confinement ratio and confined strength ratio
    plot_confined_strength()


if __name__ == "__main__":
    main()
