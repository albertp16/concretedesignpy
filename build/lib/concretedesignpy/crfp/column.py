# column_materials.py
# Computes material properties and geometric parameters for reinforced concrete columns.

import numpy as np

def compute_modulus_of_elasticity(f_c):
    """
    Computes the modulus of elasticity of concrete based on compressive strength.
    """
    return 57000 * np.sqrt(f_c) # MPa

def compute_gross_area(b, h):
    """
    Computes the gross cross-sectional area of the column.
    """
    return b * h # mm^2

def main():
    # Given data
    f_c = 4000 / 145 # Concrete compressive strength (MPa)
    b = 16 * 25.4 # Column width (mm)
    h = 16 * 25.4 # Column depth (mm)
    
    # Compute properties
    E_c = compute_modulus_of_elasticity(f_c)
    A_g = compute_gross_area(b, h)
    
    # Output results
    print(f"Modulus of elasticity of concrete: {E_c:.2f} MPa")
    print(f"Gross area of column: {A_g:.2f} mm^2")

if __name__ == "__main__":
    main()

# frp_properties.py
# Computes FRP material properties for confinement calculations.

import numpy as np

def compute_effective_stress(ffu, Ef):
    """
    Computes effective stress in FRP based on tensile stress and modulus of elasticity.
    """
    ffe1 = 0.75 * ffu
    ffe2 = 0.004 * Ef
    return min(ffe1, ffe2)

def compute_frp_ratio(n, tf, b, h):
    """
    Computes the FRP volumetric ratio.
    """
    return (2 * n * tf * (b + h)) / (b * h)

def main():
    # Given data
    ffu = 155 * 6895 # Ultimate tensile stress in FRP (kPa)
    Ef = 9326 * 6895 # Modulus of elasticity of FRP (kPa)
    n = 5  # Number of FRP plies
    tf = 0.023 * 25.4 # Thickness of one FRP ply (mm)
    b = 16 * 25.4 # Column width (mm)
    h = 16 * 25.4 # Column depth (mm)
    
    # Compute properties
    ffe = compute_effective_stress(ffu, Ef)
    rho_f = compute_frp_ratio(n, tf, b, h)
    
    # Output results
    print(f"Effective stress in FRP: {ffe:.2f} mPa")
    print(f"FRP volumetric ratio: {rho_f:.4f}")

if __name__ == "__main__":
    main()

# backbone_curve.py
# Computes and plots the force-displacement backbone curve for FRP-confined columns.

import numpy as np
import matplotlib.pyplot as plt

def generate_backbone_curve(theta_a, theta_y, theta_u, cnl, Vy):
    """
    Generates the backbone curve based on the given parameters.

    Parameters:
    theta_a - Rotation at elastic limit (radians)
    theta_y - Yield rotation (radians)
    theta_u - Ultimate rotation (radians)
    cnl     - Normalized capping strength factor
    Vy      - Yield shear force (kip)

    Returns:
    theta_values - List of rotation values (radians)
    v_values - List of corresponding shear force values (kip)
    """
    # Define backbone curve points based on Mathcad matrix
    theta_values = np.array([
        -theta_a, -theta_a, -theta_u, -theta_y, 0,
        theta_y, theta_u, theta_a, theta_a
    ])
    v_values = np.array([
        0, -cnl * Vy, -Vy, -Vy, 0,
        Vy, Vy, cnl * Vy, 0
    ])
    
    return theta_values, v_values

def plot_backbone_curve(theta_values, v_values):
    """
    Plots the backbone curve.
    """
    plt.figure(figsize=(8,6))
    plt.plot(theta_values, v_values, marker='o', linestyle='-', color='b', label="Backbone Curve (Python)")
    plt.axhline(0, color='black', linewidth=0.8)
    plt.axvline(0, color='black', linewidth=0.8)
    plt.xlabel(r'$\theta$ (radians)')
    plt.ylabel(r'$V$ (kN)')  # Labeled in kip
    plt.title('Backbone Curve (Matching Mathcad)')
    plt.legend()
    plt.grid(True)
    plt.show()

if __name__ == "__main__":
    # Given values from Mathcad (units in kip and radians)
    theta_a = 0.082   # Elastic limit rotation (radians)
    theta_y = 0.0037  # Yield rotation (radians)
    theta_u = 0.0414  # Ultimate rotation (radians)
    cnl = 0.4         # Normalized capping strength factor from document
    Vy = 25 * 4.4482  # Yield shear force (kN)

    # Generate backbone curve values
    theta_values, v_values = generate_backbone_curve(theta_a, theta_y, theta_u, cnl, Vy)

    # Plot the backbone curve
    plot_backbone_curve(theta_values, v_values)
