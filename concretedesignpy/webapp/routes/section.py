# MIT License
# Copyright (c) Albert Pamonag Engineering Consultancy

"""Section, detailing, and analysis API routes."""

from flask import Blueprint, request, jsonify

from concretedesignpy.calculators.development_length import (
    bend_and_hook_deformed_bars,
    bend_and_hook_stirrups,
)
from concretedesignpy.calculators.moment_curvature import (
    moment_curvature_analysis,
)
from concretedesignpy.calculators.rebar_layout import section_generate_rect
from concretedesignpy.calculators.alternative_inertia import (
    alternative_inertia_column_wall,
    alternative_inertia_beam_slab,
)

section_bp = Blueprint("section", __name__)


@section_bp.route("/development-length", methods=["POST"])
def development_length():
    """Hook geometry calculation."""
    data = request.get_json()
    try:
        bar_type = data.get("bar_type", "deformed")
        if bar_type == "deformed":
            result = bend_and_hook_deformed_bars(
                bar_size=float(data["bar_size"]),
                angle=int(data["angle"]),
            )
        else:
            result = bend_and_hook_stirrups(
                bar_size=float(data["bar_size"]),
                angle=int(data["angle"]),
            )
        return jsonify({"status": "success", "result": result})
    except (KeyError, ValueError, TypeError) as e:
        return jsonify({"status": "error", "message": str(e)}), 400


@section_bp.route("/moment-curvature", methods=["POST"])
def moment_curvature():
    """Moment-curvature analysis."""
    data = request.get_json()
    try:
        result = moment_curvature_analysis(
            b=float(data["b"]),
            h=float(data["h"]),
            d=float(data["d"]),
            fc=float(data["fc"]),
            fy=float(data["fy"]),
            as_tension=float(data["as_tension"]),
            es=float(data.get("es", 200000)),
        )
        return jsonify({"status": "success", "result": result})
    except (KeyError, ValueError, TypeError) as e:
        return jsonify({"status": "error", "message": str(e)}), 400


@section_bp.route("/rebar-layout", methods=["POST"])
def rebar_layout():
    """Generate rebar layout coordinates."""
    data = request.get_json()
    try:
        result = section_generate_rect(
            width=float(data["width"]),
            height=float(data["height"]),
            cover=float(data["cover"]),
            ds=float(data["ds"]),
            db=float(data["db"]),
            nx=int(data["nx"]),
            ny=int(data["ny"]),
        )
        return jsonify({"status": "success", "result": result})
    except (KeyError, ValueError, TypeError) as e:
        return jsonify({"status": "error", "message": str(e)}), 400


@section_bp.route("/alternative-inertia", methods=["POST"])
def alternative_inertia():
    """Alternative moment of inertia calculation."""
    data = request.get_json()
    try:
        member_type = data.get("member_type", "column")
        if member_type in ("column", "wall"):
            result = alternative_inertia_column_wall(
                ig=float(data["ig"]),
                ast=float(data["ast"]),
                ag=float(data["ag"]),
                mu=float(data["mu"]),
                pu=float(data["pu"]),
                h=float(data["h"]),
                po=float(data["po"]),
            )
        else:
            result = alternative_inertia_beam_slab(
                ig=float(data["ig"]),
                rho=float(data["rho"]),
                b=float(data["b"]),
                d=float(data["d"]),
            )
        return jsonify({"status": "success", "result": result})
    except (KeyError, ValueError, TypeError) as e:
        return jsonify({"status": "error", "message": str(e)}), 400
