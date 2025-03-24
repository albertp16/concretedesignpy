# MIT License
#
# Copyright (c) 2025 of Albert Pamonag
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

"""
Open-source code for reinforced concrete design calculations,
including strength reduction factor (ϕ), steel yield strain, and
beta_1 calculation for rectangular beam or column sections.

References:
- ACI 318 and NSCP (National Structural Code of the Philippines) 2015 
  (Similar to ACI for fundamental provisions, but always consult official documents).
- This script includes guard clauses and additional checks to ensure valid input ranges.

Disclaimer:
- Use this code with caution and always cross-verify with official design codes
  (ACI, NSCP, etc.). The authors assume no responsibility for use or misuse of this code.
"""

def calculate_strength_reduction_factor(
    actual_epsilon: float,
    epsilon_ty: float,
    section_type: str
) -> float:
    """
    Calculate the strength reduction factor (ϕ) based on strain compatibility.

    Parameters:
    -----------
    actual_epsilon : float
        Actual (computed) tensile strain in the reinforcement.
    epsilon_ty     : float
        Tensile strain at steel yielding = fy / Es.
    section_type   : str
        Type of section. Must be either 'spiral' (for a spirally reinforced column)
        or 'others' (for tied columns or beams).

    Returns:
    --------
    phi : float
        Strength reduction factor.

    Raises:
    -------
    ValueError
        If actual_epsilon or epsilon_ty are not numbers, if section_type is invalid,
        or if inputs are out of expected ranges.
    """

    if not isinstance(actual_epsilon, (int, float)):
        raise ValueError("actual_epsilon must be a numerical type (int or float).")
    if not isinstance(epsilon_ty, (int, float)):
        raise ValueError("epsilon_ty must be a numerical type (int or float).")
    if section_type not in ["spiral", "others"]:
        raise ValueError("section_type must be either 'spiral' or 'others'.")
    if actual_epsilon < 0:
        raise ValueError("actual_epsilon cannot be negative.")
    if epsilon_ty <= 0:
        raise ValueError("epsilon_ty (yield strain) must be positive.")

    # Tension-controlled limit for mild steel (ACI / NSCP typical): 0.003 + epsilon_ty
    epsilon_t_limit = 0.003 + epsilon_ty

    # Compression-controlled if actual_epsilon <= epsilon_ty
    # Transition region if epsilon_ty < actual_epsilon < epsilon_t_limit
    # Tension-controlled if actual_epsilon >= epsilon_t_limit

    if actual_epsilon <= epsilon_ty:
        # compression-controlled
        phi = 0.75 if section_type == "spiral" else 0.65
    elif actual_epsilon < epsilon_t_limit:
        # transition
        factor = (actual_epsilon - epsilon_ty) / 0.003
        if section_type == "spiral":
            # transitions linearly from 0.75 to 0.90
            phi = 0.75 + 0.15 * factor
        else:
            # transitions linearly from 0.65 to 0.90
            phi = 0.65 + 0.25 * factor
    else:
        # tension-controlled
        phi = 0.90

    return phi


def steel_yield_strain(fy: float, es: float) -> dict:
    """
    Calculate and report the steel yield strain.

    Parameters:
    -----------
    fy : float
        Steel yield strength (MPa or psi).
    es : float
        Modulus of elasticity of steel (MPa or psi).

    Returns:
    --------
    result : dict
        A dictionary with:
          - 'value': the yield strain value (float)
          - 'report': a formatted string describing the calculation

    Raises:
    -------
    ValueError
        If fy <= 0, es <= 0, or fy or es are not numbers.
    """
    if not isinstance(fy, (int, float)):
        raise ValueError("fy must be a numerical type (int or float).")
    if not isinstance(es, (int, float)):
        raise ValueError("es must be a numerical type (int or float).")
    if fy <= 0:
        raise ValueError("fy (yield strength) must be positive.")
    if es <= 0:
        raise ValueError("es (modulus of elasticity) must be positive.")

    epsilon_s = fy / es
    return {
        "value": epsilon_s,
        "report": f"Yield steel strain, ε_s = fy / Es = {fy:.6f} / {es:.6f} = {epsilon_s:.6f}"
    }


def calculate_beta_one(fc: float) -> dict:
    """
    Calculate the concrete beta_1 factor for compression block depth calculation.

    Beta_1 is a factor reflecting the shape of the stress block in concrete design.
    This function is based on ACI (and commonly adopted in NSCP 2015 with similar approach).
    NOTE: Always consult the applicable code for local or updated provisions.

    Provisions:
    -----------
    1. For 17 MPa <= fc <= 28 MPa:
       beta_1 = 0.85
    2. For 28 < fc < 55 MPa:
       beta_1 linearly decreases from 0.85 to 0.65 as fc increases from 28 MPa to 55 MPa.
    3. For fc >= 55 MPa:
       beta_1 = 0.65

    Parameters:
    -----------
    fc : float
        Concrete compressive strength in MPa. Must be >= 17 MPa.

    Returns:
    --------
    result : dict
        A dictionary with:
          - 'value': the calculated beta_1 value (float)
          - 'report': a formatted string describing the calculation.

    Raises:
    -------
    ValueError
        If fc < 17 MPa or fc is not a valid positive number.
    """
    if not isinstance(fc, (int, float)):
        raise ValueError("fc must be a numerical type (int or float).")
    if fc < 17:
        raise ValueError("fc < 17 MPa is out of scope for this formula (NSCP/ACI).")

    if 17 <= fc <= 28:
        beta_one = 0.85
        rep = "beta_1 = 0.85"
    elif 28 < fc < 55:
        beta_one = 0.85 - ((0.05 * (fc - 28)) / 7.0)
        rep = (
            f"beta_1 = 0.85 - (0.05 * (fc - 28)) / 7\n"
            f"       = 0.85 - (0.05 * ({fc:.2f} - 28)) / 7\n"
            f"       = {beta_one:.2f}"
        )
    else:  # fc >= 55
        beta_one = 0.65
        rep = "beta_1 = 0.65"

    # Additional safety checks or provisions (if needed):
    # For NSCP 2015, typically identical approach is used as ACI 318 for normal strength concrete.
    # High-strength concretes might have additional modifications or limitations in certain codes.
    # Always consult code if fc > 55 MPa regarding confinement, deflection checks, etc.

    return {"value": beta_one, "report": rep}


