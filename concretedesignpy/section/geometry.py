import matplotlib.pyplot as plt

# =============================================================================
# Function: compute_rebar_coordinates
# =============================================================================
def compute_rebar_coordinates(width, height, cover, ds, db, nx, ny):
    """
    Compute the coordinates for the rebar layout in a beam cross-section.

    Parameters:
        width    : Section width  [mm]
        height   : Section height [mm]
        cover    : Concrete cover (to stirrup outer face) [mm]
        ds       : Diameter of stirrups   [mm]
        db       : Diameter of main bars  [mm]
        nx       : Number of bars along top & bottom
        ny       : Number of bars along left & right

    Returns:
        dict: {
            'section_corners': 4 points (x, y) for the OUTER rectangle (unconfined section),
            'stirrup_corners': 4 points (x, y) for the STIRRUP CENTER-LINE rectangle,
            'rebar_coords':    (x, y) locations for each main rebar
        }
    """
    # --- 1) Outer rectangle corners (unconfined section) ---
    section_corners = [
        (0.0,  0.0),       # bottom-left
        (width,  0.0),     # bottom-right
        (width,  height),  # top-right
        (0.0,  height)     # top-left
    ]

    # --- 2) Stirrup center-line corners (confined section) ---
    # Offset from the outer face by (cover + 0.5 * ds)
    stirrup_offset = cover + 0.5 * ds
    stirrup_corners = [
        (stirrup_offset,             stirrup_offset),
        (width - stirrup_offset,     stirrup_offset),
        (width - stirrup_offset,     height - stirrup_offset),
        (stirrup_offset,             height - stirrup_offset)
    ]

    # --- 3) Main bar offset from the outer face ---
    # = cover + ds + 0.5 * db
    rebar_offset = cover + ds + 0.5 * db
    x_left  = rebar_offset
    x_right = width - rebar_offset
    y_bot   = rebar_offset
    y_top   = height - rebar_offset

    # --- 4) Compute rebar coordinates ---
    rebar_coords = []

    # Horizontal bars (top & bottom)
    if nx == 1:
        x_positions = [(x_left + x_right) / 2.0]
    else:
        dx = (x_right - x_left) / (nx - 1)
        x_positions = [x_left + i * dx for i in range(nx)]
    for x_pos in x_positions:
        rebar_coords.append((x_pos, y_bot))  # bottom bar
        rebar_coords.append((x_pos, y_top))  # top bar

    # Vertical bars (left & right)
    if ny == 1:
        y_positions = [(y_bot + y_top) / 2.0]
    else:
        dy = (y_top - y_bot) / (ny - 1)
        y_positions = [y_bot + i * dy for i in range(ny)]
    for y_pos in y_positions:
        rebar_coords.append((x_left,  y_pos))  # left bar
        rebar_coords.append((x_right, y_pos))  # right bar

    # Remove duplicates (corner bars appear in both sets)
    rebar_coords = list(set(rebar_coords))
    rebar_coords.sort(key=lambda pt: (pt[1], pt[0]))

    return {
        "section_corners": section_corners,
        "stirrup_corners": stirrup_corners,
        "rebar_coords": rebar_coords
    }


# =============================================================================
# Function: plot_beam_rebar
# =============================================================================
def plot_beam_rebar(width, height, section_corners, stirrup_corners, rebar_coords, show_labels=True):
    """
    Plot the beam cross-section with the rebar layout.

    Parameters:
        width           : Section width [mm]
        height          : Section height [mm]
        section_corners : List of 4 points for the outer rectangle (unconfined section)
        stirrup_corners : List of 4 points for the stirrup center-line rectangle (confined section)
        rebar_coords    : List of (x, y) coordinates for each main rebar
        show_labels     : If True, labels are added to the plot
    """
    fig, ax = plt.subplots(figsize=(5, 6))

    # --- 1) Plot the outer rectangle (unconfined section) ---
    x_vals = [p[0] for p in section_corners] + [section_corners[0][0]]
    y_vals = [p[1] for p in section_corners] + [section_corners[0][1]]
    ax.plot(x_vals, y_vals, 'k-', linewidth=2, label="Unconfined Section")

    # --- 2) Plot the stirrup center-line rectangle ---
    sx_vals = [p[0] for p in stirrup_corners] + [stirrup_corners[0][0]]
    sy_vals = [p[1] for p in stirrup_corners] + [stirrup_corners[0][1]]
    ax.plot(sx_vals, sy_vals, 'b--', linewidth=1.5, label="Stirrup (Center-Line)")

    # --- 3) Plot the rebar points ---
    rx = [p[0] for p in rebar_coords]
    ry = [p[1] for p in rebar_coords]
    ax.scatter(rx, ry, s=60, c='red', marker='o', label="Main Bars")

    # --- 4) Add labels if required ---
    if show_labels:
        for i, pt in enumerate(section_corners, start=1):
            ax.text(pt[0], pt[1], f"C{i}", color='blue', fontsize=9, ha='left', va='bottom')
        for i, pt in enumerate(stirrup_corners, start=1):
            ax.text(pt[0], pt[1], f"S{i}", color='green', fontsize=9, ha='left', va='bottom')
        for i, (xr, yr) in enumerate(rebar_coords, start=1):
            ax.text(xr, yr, f"R{i}", color='red', fontsize=8, ha='left', va='bottom')

    # --- 5) Final plot formatting ---
    ax.set_aspect('equal', adjustable='box')
    ax.set_xlim(-10, width + 10)
    ax.set_ylim(-10, height + 10)
    ax.set_xlabel("X-axis (mm)")
    ax.set_ylabel("Y-axis (mm)")
    ax.set_title("Cross-Section with Rebar Layout")
    ax.legend()
    ax.grid(True)

    plt.show()



def compute_section_properties_rect(width, height):
    """
    Compute basic section properties for a rectangular cross-section.

    Parameters:
        width  (b): Width of the rectangle  [mm]
        height (h): Height of the rectangle [mm]

    Returns:
        dict: {
            "Area":                 A,
            "FirstMomentX":         Qx,
            "FirstMomentY":         Qy,
            "MomentInertiaX":       Ix,
            "MomentInertiaY":       Iy,
            "RadiusGyrationX":      rx,
            "RadiusGyrationY":      ry,
            "SectionModulusX":      Sx,
            "SectionModulusY":      Sy
        }

    Notes on axes:
    - X axis is assumed horizontal, Y axis is vertical.
    - For the first moments (Qx, Qy), we assume the reference axes are at
      the bottom (y=0) and left (x=0) edges of the rectangle.
    - Moments of inertia (Ix, Iy) are about the centroidal axes (X-X' and Y-Y').
    - Section moduli (Sx, Sy) are based on the centroidal moments of inertia.
    """

    # 1) Area
    area = width * height

    # 2) First moment of area about:
    #    - X axis at y=0 -> Qx = ∫y dA = A*(h/2) = b*h^2 / 2
    #    - Y axis at x=0 -> Qy = ∫x dA = A*(b/2) = b^2*h / 2
    Qx = width * (height**2) / 2.0
    Qy = (width**2) * height / 2.0

    # 3) Moment of inertia about centroidal axes:
    #    Ix (about horizontal axis through centroid) = b*h^3 / 12
    #    Iy (about vertical axis through centroid)   = b^3*h / 12
    Ix = (width * height**3) / 12.0
    Iy = (width**3 * height) / 12.0

    # 4) Radius of gyration: r = sqrt(I / A)
    rx = (Ix / area)**0.5
    ry = (Iy / area)**0.5

    # 5) Section modulus: S = I / (c), where c is distance from centroid to extreme fiber
    #    For a rectangle: c_x = h/2, c_y = b/2
    Sx = Ix / (height / 2.0)  # = 2*Ix / h
    Sy = Iy / (width / 2.0)   # = 2*Iy / b

    return {
        "area": area,
        "first_moment": {
            "x" : Qx,
            "y" : Qy
        },
        "moment_interia" : {
            "x" : Ix,
            "y" : Iy
        },
        "radius_gyration" : {
            "x" : rx,
            "y" : ry
        },
        "section_modulus" : {
            "x" : Sx,
            "y" : Sy
        }
    }


# =============================================================================
# Main Execution Block
# =============================================================================
if __name__ == "__main__":
    # Define section parameters
    width   = 300.0    # mm
    height  = 600.0    # mm
    cover   = 40.0     # mm (measured to stirrup outer face)
    ds      = 12.0     # mm (stirrup diameter)
    db      = 20.0     # mm (main bar diameter)
    nx      = 3        # number of bars along top & bottom
    ny      = 4        # number of bars along left & right

    # Compute the rebar coordinates and corners
    results = compute_rebar_coordinates(width, height, cover, ds, db, nx, ny)
    # print(results)
    section = compute_section_properties_rect(width, height)
    # print(section)
    # section_corners = results["section_corners"]
    # stirrup_corners = results["stirrup_corners"]
    # rebar_coords    = results["rebar_coords"]

    # # Print the computed values
    # print("Unconfined (outer) rectangle corners:")
    # for i, corner in enumerate(section_corners, start=1):
    #     print(f"  C{i} = {corner}")

    # print("\nConfined (inner) rectangle corners (center-line):")
    # for i, corner in enumerate(stirrup_corners, start=1):
    #     print(f"  S{i} = {corner}")

    # print("\nRebar coordinates:")
    # for i, coord in enumerate(rebar_coords, start=1):
    #     print(f"  R{i} = {coord}")

    # # Plot the cross-section with the rebar layout
    # plot_beam_rebar(width, height, section_corners, stirrup_corners, rebar_coords, show_labels=True)

import numpy as np
import matplotlib.pyplot as plt


def calculate_po(fc, fy, width, height, asteel):
    """
    Calculate the axial compression load (Po) for a reinforced concrete column.
    
    Parameters:
    fc : float : Compressive strength of concrete (MPa)
    fy : float : Yield strength of steel (MPa)
    width : float : Width of the column (mm)
    height : float : Height of the column (mm)
    asteel : list : List of steel reinforcement areas (mm^2)
    
    Returns:
    P0 : float : Axial compression load (N)
    """
    # Calculate concrete contribution
    cc = 0.85 * fc * width * height
    
    # Calculate steel contributions
    cs_total = sum(as_index * (fy - 0.85 * fc) for as_index in asteel)
    
    # Total axial compression load
    value = cc + cs_total
    
    return value

def calculate_pt(fy, asteel):
    """
    Calculate the axial tensile load (Pt) for a reinforced concrete column.
    
    Parameters:
    fy : float : Yield strength of steel (MPa)
    asteel : list : List of steel reinforcement areas (mm^2)
    
    Returns:
    Pt : float : Axial tensile load (N)
    """
    # Calculate tensile contributions
    pt = sum(as_index * fy for as_index in asteel)
    
    return pt

# Example usage
fc_prime = 30  # MPa
fy = 420  # MPa
b = 600  # mm
h = 600  # mm
As_list = [2580, 1290, 1290, 2580]  # mm^2

po = calculate_po(fc_prime, fy, b, h, As_list)
pt = calculate_pt(fy,As_list)

rebar_data = [[2580,1290,1290,2580],[75,200,400,525]] ## [[steel area],[depth from top]] 600x600mm rectangular

def compute_balanced_point(rebar_data, depth, epsilon_y):
    """
    Compute the balanced neutral axis depth (C) and strain values (ε_si) for a given set of rebar data.

    Parameters:
    - rebar_data: List of two lists, first for steel area and second for depth from the top.
    - d: Total depth of the section.
    - epsilon_y: Yield strain of the steel.

    Returns:
    - C: Balanced neutral axis depth in mm.
    - epsilon_si: List of strain values for each reinforcement layer.
    """
    # Extract depths
    depths = rebar_data[1]
    
    # Compute C (neutral axis depth)
    neutral_axis = (0.003 / (0.003 + epsilon_y)) * depth
    
    # Compute strain values for each reinforcement layer
    epsilon_si = [(neutral_axis - depth_index) / neutral_axis * 0.003 for depth_index in depths]
    
    return neutral_axis, epsilon_si

if __name__ == "__main__":
    d = 600  # Total depth of section in mm
    epsilon_y = 525 / 200000  # Assuming Young's modulus of steel = 200000 MPa and yield stress = 525 MPa

    # Compute results
    compute_balanced_point(rebar_data, d, epsilon_y)

# def plot_PM_interaction(): 
#     # Given constants
#     b = 600.0  # width (mm)
#     h = 600.0  # depth (mm)
#     fc = 30.0  # MPa (concrete compressive strength)
#     fy = 420.0  # MPa (steel yield strength)
#     Ec = 30000.0  # MPa (concrete elastic modulus)
#     Es = 200000.0  # MPa (steel elastic modulus)
#     eps_cu = 0.003    # ultimate concrete strain
#     eps_y = fy / Es   # steel yield strain (~0.00207)

#     # Rebar layout (depth from top in mm and area in mm^2)
#     rebar_depths = np.array([75.0, 200.0, 400.0, 525.0])
#     rebar_areas = np.array([2580.0, 1290.0, 1290.0, 2580.0])
    
#     # Define concrete stress-strain behavior (simplified)
#     def concrete_stress(eps_c):
#         # eps_c positive for compression
#         if eps_c <= 0:      # no tension
#             return 0.0
#         if eps_c <= 0.002:  # linear up to f'c at 0.002
#             return (eps_c / 0.002) * fc
#         elif eps_c <= eps_cu:  # flat top from 0.002 to 0.003
#             return fc
#         else:
#             return 0.0  # beyond 0.003, concrete has crushed (no capacity)
    
#     # Lists to store results
#     P_values = []
#     M_values = []
    
#     # Sweep neutral axis depth from near 0 to a large value
#     for c in np.linspace(50, 2000, 100):  # mm, 50 to 2000
#         # print(c)
#         # Compute concrete compression resultant
#         # Integrate concrete stress from y=0 to y=c (or h if c > h)
#         y_max = min(c, h)
#         # print(y_max)
#         # Use small slices for integration
#         n_steps = 100
#         ys = np.linspace(0, y_max, n_steps)
#         conc_force = 0.0 # Compression forces
#         conc_moment_top = 0.0
#         for y in ys:
#             eps = eps_cu * (1 - y/c)  # strain at depth y
#             stress = concrete_stress(eps)
#             dA = b * (y_max / n_steps)  # differential area (width * slice thickness)
#             conc_force += stress * dA       # force = stress * area (N, since MPa*mm^2)
#             conc_moment_top += stress * dA * y  # moment about top = F * y

#         print(conc_force)
#         # Steel forces
#         steel_force = 0.0
#         steel_moment_top = 0.0
#         for (d_i, A_i) in zip(rebar_depths, rebar_areas):
#             eps_s = eps_cu * (1 - d_i/c)
#             # Steel stress
#             if eps_s >= eps_y:       # yielding in compression
#                 sigma_s = fy
#             elif eps_s <= -eps_y:    # yielding in tension
#                 sigma_s = -fy
#             else:                    # elastic range
#                 sigma_s = Es * eps_s
#             F_si = sigma_s * A_i
#             steel_force += F_si
#             steel_moment_top += F_si * d_i
#         # Total axial force and moment about top
#         P_int = conc_force + steel_force        # (N)
#         M_int_top = conc_moment_top + steel_moment_top  # (N·mm)
#         # Convert moment about top to moment about mid-depth (300 mm)
#         M_int_mid = M_int_top - P_int * (h/2)
#         # Store (convert to kN and kN·m for plotting convenience)
#         P_values.append(P_int / 1000.0)           # kN
#         M_values.append(abs(M_int_mid) / 1e6)     # kN·m (take abs for positive moment)
    
#     # Plot P-M interaction curve
#     fig, axs = plt.subplots(1, 2, figsize=(10,4))
#     axs[0].plot(M_values, P_values, color='blue')
#     axs[0].invert_yaxis()  # invert Y so compression P is plotted downward (optional)
#     axs[0].set_xlabel("Moment $M$ (kN·m)")
#     axs[0].set_ylabel("Axial Load $P$ (kN)")
#     axs[0].set_title("P-M Interaction Diagram")
#     # Mark balanced, pure bending, pure compression points (optional annotations)
#     # (One could compute or find indices for these special points for annotation)
    
#     # Plot stress-strain curves for concrete and steel
#     # Concrete from 0 to 0.003, Steel from -0.003 to +0.003
#     eps_conc = np.linspace(0, 0.003, 100)
#     sigma_conc = [concrete_stress(eps) for eps in eps_conc]
#     eps_steel = np.linspace(-0.003, 0.003, 200)
#     sigma_steel = []
#     for eps in eps_steel:
#         if eps >= eps_y:
#             sigma_steel.append(fy)
#         elif eps <= -eps_y:
#             sigma_steel.append(-fy)
#         else:
#             sigma_steel.append(Es * eps)
#     axs[1].plot(eps_conc, sigma_conc, label="Concrete (comp)", color='red')
#     axs[1].plot(eps_steel, sigma_steel, label="Steel (tension/comp)", color='green')
#     axs[1].axvline(x=0, color='k', linewidth=0.8)
#     axs[1].axhline(y=0, color='k', linewidth=0.8)
#     axs[1].set_xlabel("Strain $\epsilon$")
#     axs[1].set_ylabel("Stress (MPa)")
#     axs[1].set_title("Material Stress-Strain Curves")
#     axs[1].legend()
#     axs[1].grid(True)
#     plt.tight_layout()
#     return fig, axs  # return figure and axes for further use if needed

# # Call the function to compute and get the plots
# fig, axes = plot_PM_interaction()
# plt.show()  # This will display the plots when run in an appropriate environment


# =============================================================================
# Function: plot_interaction_diagram
# =============================================================================
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline

# -----------------------------
# Pure Compression and Tension
# -----------------------------
def calculate_po(fcprime, fy, b, h, as_list):
    """
    Calculate the pure axial compression capacity (Po) of a rectangular RC section.
    
    1) Compute the concrete contribution:
       cc = 0.85 * fcprime * b * h
    2) For each layer of steel As_i, compute:
       cs_i = As_i * (fy - 0.85 * fcprime)
    3) Po = cc + sum(cs_i for all layers)
    """
    # Concrete contribution
    cc = 0.85 * fcprime * b * h
    
    # Steel contribution
    cs_list = [As_i * (fy - 0.85 * fcprime) for As_i in as_list]
    
    # Sum everything
    po = cc + sum(cs_list)
    return po

def calculate_pt(fy, as_list):
    """
    Calculate the pure axial tension capacity (Pt) of the same rectangular section.
    For each layer of steel As_i, compute:
       ts_i = As_i * fy
    and sum them.
    """
    ts_list = [As_i * fy for As_i in as_list]
    pt = sum(ts_list)
    return pt

def generate_as_list(asbar, n_list):
    """
    Given a 'base bar area' (asbar) and a list of bar counts (n_list),
    produce a list of steel areas [as1, as2, ...] = asbar*n1, asbar*n2, ...
    Automatically ignores any layer with n_i == 0.
    """
    return [asbar * n_i for n_i in n_list if n_i > 0]

def generate_depth_list(d_list):
    """
    Filter out depths that are 0 or unused.
    """
    return [d_i for d_i in d_list if d_i > 0]

# -----------------------------
# P–M Points Calculations
# -----------------------------
def points_on_pm_diagram(fcprime, fy, b, h, betafactor, epsgamma, Es, as_list, d_list):
    """
    Computes three special points on the P–M diagram (Balanced, Above-Balanced, Below-Balanced).
    """
    d_last = d_list[-1]

    # Define c-values
    c_bal   = (0.003 / (0.003 + epsgamma)) * d_last
    c_above = d_last
    c_below = (0.003 / (0.003 + 2.0 * epsgamma)) * d_last

    def compute_pm(c_value, force_last_layer_tension=False):
        """
        Computes P, M, and phi for a given neutral axis depth c_value.
        """
        # 1) Concrete block
        cc_new = 0.85 * betafactor * fcprime * b * c_value

        compression_list = []
        tension_list = []
        f_s_i_list = []

        # Define n_layers here
        n_layers = len(as_list)

        # Compute for each layer
        for i in range(n_layers):
            As_i = as_list[i]
            di   = d_list[i]

            # 2) Compute strain and fs
            if i == (n_layers - 1) and force_last_layer_tension:
                # Force last layer in full tension (-fy) for Balanced Case
                f_s_i = -fy
            elif abs(c_value - di) < 1e-9:
                # Exactly on the neutral axis
                eps_si = 0.0
                f_s_i = 0.0
            else:
                eps_si = ((c_value - di) / c_value) * 0.003
                f_s_i = Es * eps_si
                # Clamp fs based on its sign
                if f_s_i >= 0:
                    if f_s_i > fy:
                        f_s_i = fy
                else:
                    if abs(f_s_i) > fy:
                        f_s_i = -fy

            f_s_i_list.append(f_s_i)

            # 3) Compute c_s_i first, decide compression or tension
            c_s_i = As_i * (f_s_i - 0.85 * fcprime)
            if c_s_i > 0:
                compression_list.append((c_s_i, di))
            else:
                # If c_s_i is negative, disregard and compute tension steel instead
                t_s_i = As_i * f_s_i
                tension_list.append((t_s_i, di))

        # 5) Axial load
        sum_c = sum([comp[0] for comp in compression_list])
        sum_t = sum([ten[0] for ten in tension_list])
        p_val = cc_new + sum_c + sum_t

        # 6) Moment about mid-height
        M_conc    = cc_new * ((h / 2.0) - (betafactor * c_value) / 2.0)
        M_c_steel = sum([c_s_i * ((h / 2.0) - di) for (c_s_i, di) in compression_list])
        M_t_steel = sum([t_s_i * (di - (h / 2.0)) for (t_s_i, di) in tension_list])
        m_val = M_conc + M_c_steel - M_t_steel

        # 7) reduction factor
        phi_val = 0.003 / c_value if c_value != 0.0 else 0.0

        return (c_value, p_val, m_val, phi_val, f_s_i_list)

    # Compute for each special point
    cB, pB, mB, phiB, fB = compute_pm(c_bal,  force_last_layer_tension=True)
    cA, pA, mA, phiA, fA = compute_pm(c_above)
    cL, pL, mL, phiL, fL = compute_pm(c_below)

    return {
        "balanced": {"c": cB, "P": pB, "M": mB, "phi": phiB, "fs_list": fB},
        "above":    {"c": cA, "P": pA, "M": mA, "phi": phiA, "fs_list": fA},
        "below":    {"c": cL, "P": pL, "M": mL, "phi": phiL, "fs_list": fL},
    }

# -----------------------------
# Plotting the Diagram
# -----------------------------
def plot_pm_points_with_curve(po, pt, mB, pB, mA, pA, mL, pL, fcprime, asbar, fy):
    """
    The five points are:
      1) Pure Tension: (0, -Pt) [from Cell 1]
      2) Below-Balanced: (M_bb, P_bb) [from Cell 2]
      3) Balanced: (M_bal, P_bal) [from Cell 2]
      4) Above-Balanced: (M_ab, P_ab) [from Cell 2]
      5) Pure Compression: (0, +Po) [from Cell 1]
    
    The curve is obtained by interpolating these points with a parametric spline.
    """
    # Convert units
    Po_kN = po / 1e3        # Pure compression, positive
    Pt_kN = -(pt / 1e3)      # Pure tension, negative
    mB_kNm = mB / 1e6
    pB_kN  = pB / 1e3
    mA_kNm = mA / 1e6
    pA_kN  = pA / 1e3
    mL_kNm = mL / 1e6
    pL_kN  = pL / 1e3
    
    # Define the key points (in a parametric t-scale)
    t_vals = np.array([0, 0.25, 0.5, 0.75, 1])
    M_points = np.array([0, mL_kNm, mB_kNm, mA_kNm, 0])
    P_points = np.array([Pt_kN, pL_kN, pB_kN, pA_kN, Po_kN])

    # Create a fine parameter for smooth interpolation
    t_fine = np.linspace(0, 1, 300)
    spline_M = make_interp_spline(t_vals, M_points, k=3)
    spline_P = make_interp_spline(t_vals, P_points, k=3)
    M_fine = spline_M(t_fine)
    P_fine = spline_P(t_fine)
    
    # Plotting
    plt.figure(figsize=(7,7))
    plt.plot(M_fine, P_fine, lw=2)
    
    # Plot the five key points
    plt.plot(0, Pt_kN, 'bo', ms=8, label="Pure Tension (Pt)")
    plt.text(5, Pt_kN - 100, f"(0.0, {Pt_kN:.1f})", color='blue')

    plt.plot(mL_kNm, pL_kN, 'cs', ms=8, label="Below-Balanced")
    plt.text(mL_kNm + 20, pL_kN, f"({mL_kNm:.1f}, {pL_kN:.1f})", color='c')

    plt.plot(mB_kNm, pB_kN, 'ks', ms=8, label="Balanced")
    plt.text(mB_kNm + 20, pB_kN, f"({mB_kNm:.1f}, {pB_kN:.1f})", color='black')

    plt.plot(mA_kNm, pA_kN, 'ys', ms=8, label="Above-Balanced")
    plt.text(mA_kNm + 20, pA_kN, f"({mA_kNm:.1f}, {pA_kN:.1f})", color='goldenrod')

    plt.plot(0, Po_kN, 'ro', ms=8, label="Pure Compression (Po)")
    plt.text(5, Po_kN + 100, f"(0.0, {Po_kN:.1f})", color='red')

    # Axes, labels, legend
    plt.axhline(0, color='k', ls='--', lw=0.8)
    plt.axvline(0, color='k', ls='--', lw=0.8)
    plt.xlabel("Moment, M (kN·m)")
    plt.ylabel("Axial Load, P (kN)")
    plt.title("P–M Interaction Diagram")
    plt.grid(True)
    plt.legend(loc='best')

    # Cross Section Information
    text_str = (
        "Cross-section properties:\n"
        f"f'c = {fcprime} MPa\n"
        f"As  = {asbar} mm^2\n"
        f"fy  = {fy} MPa"
    )

    plt.text(0.5, 0.5, text_str,
             fontsize=10, ha='center', va='center',
             transform=plt.gca().transAxes,
             bbox=dict(facecolor='white', alpha=0.7, boxstyle="round"))

    plt.show()

if __name__ == "__main__":
    # Section Properties and Material Data
    fcprime    = 30.0       # MPa
    fy         = 420.0      # MPa
    Ec         = 30000.0    # MPa
    Es         = 200000.0   # MPa
    betafactor = 0.85
    epsgamma   = 0.00207
    # barsize = 28
    # asbar = math.pi * math.pow (barsize/2,2)
    asbar      = 645.0      # mm^2 (area of one bar)
    
    # Cross-section geometry
    b = 600.0  # mm
    h = 600.0  # mm
    
    # Reinforcement layout
    d_list = [75.0, 200.0, 400.0, 525.0]  # Depths from top fiber (mm)
    n_list = [4, 2, 2, 4]                 # Number of bars in each layer
    
    # Generate reinforcement and depth lists
    as_list = generate_as_list(asbar, n_list)
    d_list  = generate_depth_list(d_list)
    
    # Compute pure axial capacities (Cell 1)
    po = calculate_po(fcprime, fy, b, h, as_list)
    pt = calculate_pt(fy, as_list)
    
    
    # Compute P–M diagram key points (Cell 2)
    results = points_on_pm_diagram(fcprime, fy, b, h, betafactor, epsgamma, Es, as_list, d_list)
    pB = results["balanced"]["P"]
    mB = results["balanced"]["M"]
    pA = results["above"]["P"]
    mA = results["above"]["M"]
    pL = results["below"]["P"]
    mL = results["below"]["M"]
    
    
    # Plot the P–M interaction diagram (Cell 3)
    plot_pm_points_with_curve(po, pt, mB, pB, mA, pA, mL, pL, fcprime, asbar, fy)