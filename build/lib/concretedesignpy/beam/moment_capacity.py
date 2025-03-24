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
its associated factors following (primarily) RC Concrete Moment Capacity.
"""

import numpy as np
import concretedesignpy as gen  

def process_rebar_data(rebar_list):
    """
    Parameters
    ----------
    rebar_list : list of dict
        A list of dictionaries, where each dictionary must have the following keys:
            - "d"    (float): The distance of the rebar from a reference point in millimeters.
            - "diam" (float): The diameter of the rebar in mm.
            - "num"  (int)  : The number of rebars of that diameter and location.
        Example:
            [
                {"d": 58,  "diam": 16, "num": 2},
                {"d": 342, "diam": 16, "num": 3},
                ...
            ]
    
    Returns
    -------
    result : dict
        A dictionary containing:
        
        "value" : list of dict
            A list of dictionaries with the same fields as input plus an additional key:
                - "as" (float): The computed cross-sectional area of the rebars.
            
        "report" : str
            A human-readable multi-line string summarizing the computed values for each entry.
    """
    # Extract columns into NumPy arrays.
    ds    = np.array([r["d"]    for r in rebar_list], dtype=float)
    diams = np.array([r["diam"] for r in rebar_list], dtype=float)
    nums  = np.array([r["num"]  for r in rebar_list],  dtype=int)

    print("Distances:", ds)
    print("Diameters:", diams)
    print("Counts:   ", nums)

    area_diam_vectorized = np.vectorize(gen.area_diam)
    steel_area_vectorized = np.vectorize(gen.steel_area)
    
    area_each = area_diam_vectorized(diams)
    areas = steel_area_vectorized(area_each, nums)

    # Build output list and report.
    processed = []
    lines = ["rebar data (d, diam, num, as):"]

    for i in range(len(rebar_list)):
        item = {
            "d": ds[i],
            "diam": diams[i],
            "num": nums[i],
            "as": areas[i]
        }
        processed.append(item)
        lines.append(
            f"  d={ds[i]:7.2f}, "
            f"diam={diams[i]:5.2f}, "
            f"num={nums[i]:2d}, "
            f"as={areas[i]:10.2f}"
        )

    return {
        "value": processed,
        "report": "\n".join(lines)
    }


# def calculate_rebar_forces(
#     rebar_list, beam_height, beam_width,
#     fc, fy, es, ecu, beta_one
# ):
#     '''
#     rebar in list in mm and sq.mm 
#     beam height in mm 
#     beam width in mm 
#     fc in MPa 
#     fy in MPa
#     es in unitless
#     ecu in unitless
#     beta_one in unitless
#     '''
#     dt = beam_height
#     current_x = dt
#     iter_num = 50
#     iter_step = dt / iter_num
#     update_step = True
#     alpha_2 = 0.85

#     x_output = None
#     solution_found = False

#     iteration_data = []
#     i = 0
#     pass_count = 1

#     # final forces
#     fcc_final = 0.0
#     fcsi_final = 0.0
#     fsi_final = 0.0
#     ratio_final = None

#     while i < iter_num:
#         if i == 0:
#             x = current_x
#         else:
#             x = current_x - iter_step
#             current_x = x

#         # compute compression & tension rebar
#         fcsi_rebar = 0.0
#         fsi_rebar = 0.0
#         for bar in rebar_list:
#             d_i = bar["d"]
#             as_i = bar["as"]
#             if d_i < x:
#                 # compression zone
#                 esci = ((x - d_i)/x)*ecu
#                 f_sci = min(esci*es, fy)
#                 fi = as_i*f_sci
#                 fcsi_rebar += fi
#             else:
#                 # tension zone
#                 esti = ((d_i - x)/x)*ecu
#                 f_sti = min(esti*es, fy)
#                 fi = as_i*f_sti
#                 fsi_rebar += fi

#         # concrete block
#         a_comp = beta_one * x * beam_width
#         f_concrete = alpha_2 * fc * a_comp
#         fci = f_concrete + fcsi_rebar
#         fsi = fsi_rebar
#         ratio = fci/fsi if fsi != 0 else float('inf')

#         iteration_data.append({
#             "pass": pass_count,
#             "iteration": i+1,
#             "x": x,
#             "fc_concrete_kN": f_concrete/1000.0,
#             "fc_rebar_kN": fcsi_rebar/1000.0,
#             "fs_rebar_kN": fsi_rebar/1000.0,
#             "fci_kN": (f_concrete + fcsi_rebar)/1000.0,
#             "ratio": ratio
#         })

#         if not update_step:
#             if fci <= fsi:
#                 x_output = x
#                 fcc_final = f_concrete
#                 fcsi_final = fcsi_rebar
#                 fsi_final = fsi_rebar
#                 ratio_final = ratio
#                 solution_found = True
#                 break
#         else:
#             if fci <= fsi:
#                 update_step = False
#                 # store final
#                 x_output = x
#                 fcc_final = f_concrete
#                 fcsi_final = fcsi_rebar
#                 fsi_final = fsi_rebar
#                 ratio_final = ratio
#                 solution_found = True

#                 # revert x by adding step
#                 current_x += iter_step
#                 # refine step
#                 iter_step = iter_step / iter_num
#                 pass_count += 1
#                 i = -1  # so next loop iteration will be i=0
#         i += 1

#     return {
#         "iteration_data": iteration_data,
#         "neutral_axis": x_output,
#         "fc_concrete": fcc_final,
#         "fc_rebar": fcsi_final,
#         "fs_rebar": fsi_final,
#         "ratio": ratio_final
#     }

def calculate_rebar_forces(
    rebar_list,
    beam_height,
    beam_width,
    fc,
    fy,
    es,
    ecu,
    beta_one,
    alpha_2=0.85,
    max_outer_iterations=50
):
    """
    Calculate the neutral axis depth, concrete compressive force, steel compressive force,
    and steel tensile force in a rectangular RC section using an iterative approach.

    Parameters
    ----------
    rebar_list : list of dict
        Each dict must have:
          - d  : float, the distance from the top of the beam (mm).
          - as : float, the cross-sectional area of that rebar group (mm²).
        Example: [{"d": 50.0, "as": 804.25}, ...]

    beam_height : float
        Total height of the beam cross section (mm).
    beam_width : float
        Width (b) of the rectangular beam (mm).
    fc : float
        Concrete compressive strength (MPa).
    fy : float
        Steel yield strength (MPa).
    es : float
        Steel elastic modulus or slope factor to convert strain to stress (MPa).
        (If you are truly using E_s in MPa, note that es * strain yields stress in MPa.)
        Often E_s ~ 200,000 MPa in practice.
    ecu : float
        Ultimate concrete strain (unitless). Typical ~ 0.003 or 0.0035.
    beta_one : float
        The Whitney stress block parameter (unitless). Typical range ~ 0.65 - 0.85,
        depending on fc.
    alpha_2 : float, optional
        Factor for the equivalent concrete stress block (unitless). Default=0.85.
    max_outer_iterations : int, optional
        Maximum number of outer iterations for refining the step size. Default=50.

    Returns
    -------
    dict
        {
            "iteration_data": List of iteration records (each is a dict with pass, iteration,
                              x, fc_concrete_kN, fc_rebar_kN, fs_rebar_kN, fci_kN, ratio),
            "neutral_axis": float or None,
            "fc_concrete": float,
            "fc_rebar": float,
            "fs_rebar": float,
            "ratio": float or None
        }

    Notes
    -----
    - Forces are returned in N internally, then shown in kN in iteration_data.
    - The solver tries to find x such that total compression = total tension by
      iteratively decreasing x from beam_height down to near 0, then refining step
      when crossing from compression > tension to compression <= tension.
    - Because es, ecu, and beta_one are dimensionless, ensure es indeed produces
      stress in MPa when multiplied by strain.
    """

    def _compute_forces_at_x(x):
        """
        Compute total compression (f_comp) and tension (f_tens) at a given neutral-axis depth x.

        Returns
        -------
        f_comp : float (N)
            Sum of concrete compression + compression steel.
        f_tens : float (N)
            Sum of tension steel.
        """
        # 1) Steel forces
        fcsi_rebar = 0.0  # compression steel
        fsi_rebar = 0.0   # tension steel

        for bar in rebar_list:
            d_i = bar["d"]
            as_i = bar["as"]

            # Strain in steel at that layer
            if d_i < x:
                # Compression zone
                strain_s = ((x - d_i) / x) * ecu
                stress_s = min(strain_s * es, fy)  # limit to fy
                force_s = as_i * stress_s
                fcsi_rebar += force_s
            else:
                # Tension zone
                strain_s = ((d_i - x) / x) * ecu
                stress_s = min(strain_s * es, fy)  # limit to fy
                force_s = as_i * stress_s
                fsi_rebar += force_s

        # 2) Concrete block
        # a = beta_one * x
        # Compressive stress block area = a * beam_width
        # Compressive force in concrete (Whitney block)
        a_block = beta_one * x * beam_width
        f_concrete = alpha_2 * fc * a_block  # N (since fc is MPa, area is mm²)

        # Total compression
        f_comp = f_concrete + fcsi_rebar
        # Total tension
        f_tens = fsi_rebar

        return f_comp, f_tens, f_concrete, fcsi_rebar, fsi_rebar

    # Initial setup
    iteration_data = []
    pass_count = 1

    # We start from x close to the full beam height (dt)
    x_current = beam_height
    # We'll move in roughly equal steps
    step_size = beam_height / max_outer_iterations

    # Variables to store final results
    x_output = None
    fcc_final = 0.0
    fcsi_final = 0.0
    fsi_final = 0.0
    ratio_final = None

    # Controls whether we are still in the "coarse" search or in refinement mode
    update_step = True

    # Outer loop: up to max_outer_iterations
    # The approach:
    #   - Decrease x in steps
    #   - If we find that compression <= tension at some step,
    #     we revert x a step, refine the step size, and repeat.
    for i in range(max_outer_iterations):
        # For the first iteration, use x_current as is;
        # afterwards, keep decreasing by step_size
        if i > 0:
            x_current -= step_size

        (f_comp,
         f_tens,
         f_concrete,
         fc_rebar,
         fs_rebar) = _compute_forces_at_x(x_current)

        ratio = f_comp / f_tens if f_tens != 0 else float("inf")

        iteration_data.append({
            "pass": pass_count,
            "iteration": i + 1,
            "x": x_current,
            "fc_concrete_kN": f_concrete / 1000.0,
            "fc_rebar_kN": fc_rebar / 1000.0,
            "fs_rebar_kN": fs_rebar / 1000.0,
            "fci_kN": f_comp / 1000.0,  # total compression in kN
            "ratio": ratio
        })

        # If we're not refining the step anymore, check if we've found a solution
        if not update_step:
            # Once we see compression <= tension, we consider it "solved"
            if f_comp <= f_tens:
                x_output = x_current
                fcc_final = f_concrete
                fcsi_final = fc_rebar
                fsi_final = fs_rebar
                ratio_final = ratio
                break
        else:
            # If we are still coarse-searching...
            if f_comp <= f_tens:
                # We just crossed from compression > tension to compression <= tension
                # Save the solution at this crossing
                x_output = x_current
                fcc_final = f_concrete
                fcsi_final = fc_rebar
                fsi_final = fs_rebar
                ratio_final = ratio

                # Stop coarse searching; refine step
                update_step = False

                # Revert x by adding step (go back to the last iteration's x)
                x_current += step_size

                # Now make step_size smaller (refinement)
                step_size /= max_outer_iterations

                # We'll start again from i = -1 so that next loop iteration is i=0
                pass_count += 1
                i = -1
                continue  # jump to next iteration after resetting i

    return {
        "iteration_data": iteration_data,
        "neutral_axis": x_output,
        "fc_concrete": fcc_final,
        "fc_rebar": fcsi_final,
        "fs_rebar": fsi_final,
        "ratio": ratio_final
    }


# def example_usage():
#     """
#     Example usage of the calculate_rebar_forces function
#     with a hypothetical beam and reinforcement layout.
#     """

#     # Sample rebar layout (mm, mm²):
#     # Suppose we have:
#     #   - One layer of compression bars near top: d ~ 50 mm from top, total As=400 mm²
#     #   - Two layers of tension bars near bottom:
#     #       * First layer at d ~ 450 mm: As=804 mm²
#     #       * Second layer at d ~ 470 mm: As=804 mm²
#     #
#     # Distances (d) are measured from the top of the beam.
#     # Rebar area in mm² is the sum of the bars in that layer.
#     rebar_list = [
#         {"d": 50.0,  "as": 400.0},  # Compression bars
#         {"d": 450.0, "as": 804.0},  # Tension bars (layer 1)
#         {"d": 470.0, "as": 804.0},  # Tension bars (layer 2)
#     ]

#     # Beam dimensions (mm)
#     beam_height = 500.0
#     beam_width = 300.0

#     # Material properties
#     fc = 30.0       # MPa (concrete compressive strength)
#     fy = 400.0      # MPa (steel yield strength)
#     es = 200000.0   # MPa (steel modulus of elasticity; typical ~ 200,000 MPa)
#     ecu = 0.003     # ultimate concrete strain
#     beta_one = 0.85 # Whitney stress block parameter
#     alpha_2 = 0.85  # factor for the equivalent concrete stress block

#     result = calculate_rebar_forces(
#         rebar_list=rebar_list,
#         beam_height=beam_height,
#         beam_width=beam_width,
#         fc=fc,
#         fy=fy,
#         es=es,
#         ecu=ecu,
#         beta_one=beta_one,
#         alpha_2=alpha_2,
#         max_outer_iterations=50
#     )

#     # Print the final solution
#     print("=== FINAL SOLUTION ===")
#     print(f"Neutral axis depth (x)      : {result['neutral_axis']:.2f} mm")
#     print(f"Concrete compression (N)    : {result['fc_concrete']:.2f}")
#     print(f"Steel compression (N)       : {result['fc_rebar']:.2f}")
#     print(f"Steel tension (N)           : {result['fs_rebar']:.2f}")
#     print(f"Ratio (C/T)                 : {result['ratio']:.3f}\n")

#     # Optionally, inspect the iteration data
#     print("=== ITERATION DATA (showing only the last few) ===")
#     for row in result["iteration_data"][-5:]:  # last 5 entries
#         print(
#             f"Pass {row['pass']:2d} | "
#             f"Iter {row['iteration']:2d} | "
#             f"x={row['x']:.2f} mm | "
#             f"Fconcrete={row['fc_concrete_kN']:.3f} kN | "
#             f"FcompSteel={row['fc_rebar_kN']:.3f} kN | "
#             f"FtensSteel={row['fs_rebar_kN']:.3f} kN | "
#             f"Ratio={row['ratio']:.3f}"
#         )

# # Run the example
# if __name__ == "__main__":
#     example_usage()
