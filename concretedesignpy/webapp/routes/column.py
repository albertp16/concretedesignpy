# MIT License
# Copyright (c) Albert Pamonag Engineering Consultancy

"""Column calculation API routes."""

from flask import Blueprint, request, jsonify

from concretedesignpy.calculators.column_interaction import (
    generate_interaction_diagram,
)
from concretedesignpy.calculators.column_flexural import (
    check_min_flexural_strength,
)

column_bp = Blueprint("column", __name__)


@column_bp.route("/interaction", methods=["POST"])
def interaction_diagram():
    """Generate P-M interaction diagram data."""
    data = request.get_json()
    try:
        result = generate_interaction_diagram(
            fc=float(data["fc"]),
            fy=float(data["fy"]),
            b=float(data["b"]),
            h=float(data["h"]),
            n_bars=int(data["n_bars"]),
            d_bar=float(data["d_bar"]),
            n_bars_side=int(data.get("n_bars_side", 0)),
            cover=float(data.get("cover", 40)),
            c_step=float(data.get("c_step", 5)),
        )
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
