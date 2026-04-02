# MIT License
# Copyright (c) Albert Pamonag Engineering Consultancy

"""Joint shear API routes."""

from flask import Blueprint, request, jsonify

from concretedesignpy.calculators.joint_shear import joint_shear_check

joint_bp = Blueprint("joint", __name__)


@joint_bp.route("/shear", methods=["POST"])
def joint_shear():
    """Joint shear verification."""
    data = request.get_json()
    try:
        result = joint_shear_check(
            ve=float(data["ve"]),
            as1=float(data["as1"]),
            n_bars1=int(data["n_bars1"]),
            as2=float(data["as2"]),
            n_bars2=int(data["n_bars2"]),
            fy=float(data["fy"]),
            fc=float(data["fc"]),
            beam_width=float(data["beam_width"]),
            joint_depth=float(data["joint_depth"]),
            perpendicular_dist=float(data.get("perpendicular_dist", 0)),
            joint_config=int(data.get("joint_config", 1)),
            lamda=float(data.get("lamda", 1.0)),
            phi=float(data.get("phi", 0.85)),
        )
        return jsonify({"status": "success", "result": result})
    except (KeyError, ValueError, TypeError) as e:
        return jsonify({"status": "error", "message": str(e)}), 400
