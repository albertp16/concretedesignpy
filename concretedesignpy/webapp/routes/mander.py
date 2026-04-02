# MIT License
# Copyright (c) Albert Pamonag Engineering Consultancy

"""Mander confined concrete API routes."""

from flask import Blueprint, request, jsonify

from concretedesignpy.calculators.mander import (
    confined_strength_ratio,
    confined_stress_strain,
)

mander_bp = Blueprint("mander", __name__)


@mander_bp.route("/confined", methods=["POST"])
def confined():
    """Compute confined concrete strength and strain."""
    data = request.get_json()
    try:
        result = confined_stress_strain(
            fc=float(data["fc"]),
            fy_transverse=float(data["fy_transverse"]),
            b=float(data["b"]),
            h=float(data["h"]),
            cover=float(data["cover"]),
            db_main=float(data["db_main"]),
            db_tie=float(data["db_tie"]),
            n_bars_x=int(data["n_bars_x"]),
            n_bars_y=int(data["n_bars_y"]),
            s_tie=float(data["s_tie"]),
            n_legs_x=int(data.get("n_legs_x", 2)),
            n_legs_y=int(data.get("n_legs_y", 2)),
        )
        return jsonify({"status": "success", "result": result})
    except (KeyError, ValueError, TypeError) as e:
        return jsonify({"status": "error", "message": str(e)}), 400


@mander_bp.route("/ratio", methods=["POST"])
def ratio():
    """Compute confined strength ratio only."""
    data = request.get_json()
    try:
        result = confined_strength_ratio(
            fpco=float(data["fpco"]),
            fpl=float(data["fpl"]),
        )
        return jsonify({"status": "success", "result": {"ratio": round(result, 4)}})
    except (KeyError, ValueError, TypeError) as e:
        return jsonify({"status": "error", "message": str(e)}), 400
