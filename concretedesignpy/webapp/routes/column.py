# MIT License
# Copyright (c) Albert Pamonag Engineering Consultancy

"""Column calculation API routes."""

from flask import Blueprint, request, jsonify

from concretedesignpy.calculators.column_interaction import (
    generate_interaction_diagram,
    check_capacity,
)
from concretedesignpy.calculators.column_flexural import (
    check_min_flexural_strength,
)
from concretedesignpy.calculators.rebar_layout import section_generate_rect
from concretedesignpy.calculators.diagrams import (
    matplotlib_interaction_diagram,
    svg_rebar_section,
)

column_bp = Blueprint("column", __name__)


@column_bp.route("/interaction", methods=["POST"])
def interaction_diagram():
    """Generate P-M interaction diagram data with optional rebar layout."""
    data = request.get_json()
    try:
        b = float(data["b"])
        h = float(data["h"])
        d_bar = float(data.get("d_bar", 25))
        cover = float(data.get("cover", 40))
        confinement = data.get("confinement", "tied")

        # Manual bar input: list of {area, y, z} objects
        # y = distance from center (mm), z = distance from bottom (mm)
        manual_bars = data.get("manual_bars")
        bar_coords = None
        bar_areas = None

        if manual_bars and len(manual_bars) > 0:
            # Convert z-from-bottom to distance-from-top (compression face)
            bar_coords = [h - float(bar["z"]) for bar in manual_bars]
            bar_areas = [float(bar["area"]) for bar in manual_bars]
            n_bars = len(manual_bars)
        else:
            n_bars = int(data["n_bars"])

        result = generate_interaction_diagram(
            fc=float(data["fc"]),
            fy=float(data["fy"]),
            b=b, h=h,
            n_bars=n_bars,
            d_bar=d_bar,
            n_bars_side=int(data.get("n_bars_side", 0)),
            cover=cover,
            confinement=confinement,
            n_points=int(data.get("n_points", 32)),
            bar_coords=bar_coords,
            bar_areas=bar_areas,
        )

        # Load combination checks
        load_combos = data.get("load_combos")
        pu_demand = data.get("pu_demand")
        mu_demand = data.get("mu_demand")
        demand_for_chart = None

        if load_combos and len(load_combos) > 0:
            combo_checks = []
            for lc in load_combos:
                pu = float(lc["pu"])
                mu = abs(float(lc["mu"]))
                cap = check_capacity(result, pu, mu)
                combo_checks.append({
                    "name": lc.get("name", ""),
                    "pu": pu,
                    "mu": mu,
                    "phi_mn_capacity": cap["phi_mn_capacity"],
                    "dc_ratio": cap["dc_ratio"],
                    "status": cap["status"],
                })
            result["load_combo_checks"] = combo_checks
            # Use the critical combo for chart demand point
            critical = max(combo_checks, key=lambda c: abs(c["mu"]))
            demand_for_chart = (critical["pu"], critical["mu"])
        elif pu_demand is not None and mu_demand is not None:
            capacity = check_capacity(
                result,
                float(pu_demand),
                float(mu_demand),
            )
            result["demand_check"] = capacity
            demand_for_chart = (float(pu_demand), float(mu_demand))

        # Matplotlib interaction diagram (base64 PNG)
        section_label = f"{b:.0f}\u00d7{h:.0f} mm | {n_bars}-\u03c6{d_bar:.0f} | f'c={float(data['fc'])}, fy={float(data['fy'])} MPa"
        result["chart_pm"] = matplotlib_interaction_diagram(
            result["points"],
            demand=demand_for_chart if not load_combos else None,
            demand_points=combo_checks if load_combos and len(load_combos) > 0 else None,
            section_label=section_label,
        )

        # SVG rebar layout (if nx/ny provided)
        nx = int(data.get("nx", 0))
        ny = int(data.get("ny", 0))
        if nx > 0 and ny > 0:
            tie = 10 if d_bar <= 32 else 12
            layout = section_generate_rect(
                width=b, height=h, cover=cover,
                ds=tie, db=d_bar, nx=nx, ny=ny,
            )
            result["svg_rebar"] = svg_rebar_section(layout, b, h)
            result["rebar_layout"] = layout

        # Strip per-point steel_forces from points to reduce payload
        # (keep only the balanced point detail sent separately)
        for pt in result["points"]:
            pt.pop("steel_forces", None)

        return jsonify({"status": "success", "result": result})
    except (KeyError, ValueError, TypeError) as e:
        return jsonify({"status": "error", "message": str(e)}), 400


@column_bp.route("/flexural", methods=["POST"])
def flexural_check():
    """Check column minimum flexural strength ratio."""
    data = request.get_json()
    try:
        result = check_min_flexural_strength(
            mn_col_top=float(data["mn_col_top"]),
            mn_col_bot=float(data["mn_col_bot"]),
            mn_beam=float(data["mn_beam"]),
            limit=float(data.get("limit", 1.2)),
        )
        return jsonify({"status": "success", "result": result})
    except (KeyError, ValueError, TypeError) as e:
        return jsonify({"status": "error", "message": str(e)}), 400
