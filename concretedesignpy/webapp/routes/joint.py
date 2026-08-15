# MIT License
# Copyright (c) Albert Pamonag Engineering Consultancy

"""Joint shear API routes."""

from flask import Blueprint, request, jsonify

from concretedesignpy.calculators.joint_shear import joint_shear_check
from concretedesignpy.calculators.beam_report import joint_shear_report

joint_bp = Blueprint("joint", __name__)


def _joint_kwargs(data):
    """Shared payload mapping. phi is deliberately absent: Section 21.2.4.4
    fixes it at 0.85 for special-moment-frame joints, so it is not a field a
    request can carry."""
    return dict(
        ve=float(data["ve"]),
        as1=float(data["as1"]),
        n_bars1=int(data["n_bars1"]),
        as2=float(data["as2"]),
        n_bars2=int(data["n_bars2"]),
        fy=float(data["fy"]),
        fc=float(data["fc"]),
        beam_width=float(data["beam_width"]),
        joint_depth=float(data["joint_depth"]),
        column_width=float(data["column_width"]),
        perpendicular_dist=float(data.get("perpendicular_dist", 0)),
        joint_config=int(data.get("joint_config", 1)),
        lamda=float(data.get("lamda", 1.0)),
    )


@joint_bp.route("/shear", methods=["POST"])
def joint_shear():
    """Joint shear verification."""
    data = request.get_json()
    try:
        return jsonify({"status": "success",
                        "result": joint_shear_check(**_joint_kwargs(data))})
    except (KeyError, ValueError, TypeError) as e:
        return jsonify({"status": "error", "message": str(e)}), 400


@joint_bp.route("/shear-report", methods=["POST"])
def joint_shear_report_route():
    """Printable calculation sheet for special-moment-frame joint shear."""
    data = request.get_json()
    try:
        return jsonify({"status": "success",
                        "result": joint_shear_report(**_joint_kwargs(data))})
    except (KeyError, ValueError, TypeError) as e:
        return jsonify({"status": "error", "message": str(e)}), 400
