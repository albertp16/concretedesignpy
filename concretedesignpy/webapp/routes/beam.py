# MIT License
# Copyright (c) Albert Pamonag Engineering Consultancy

"""Beam calculation API routes."""

from flask import Blueprint, request, jsonify

from concretedesignpy.calculators.beam_moment import calculate_beam_moment
from concretedesignpy.calculators.beam_shear import (
    compute_concrete_shear_strength,
    compute_steel_shear_strength,
    compute_shear_spacing,
)
from concretedesignpy.calculators.beam_torsion import torsion_design
from concretedesignpy.calculators.beam_deflection import deflection_computation

beam_bp = Blueprint("beam", __name__)


@beam_bp.route("/moment", methods=["POST"])
def beam_moment():
    """Calculate beam moment capacity."""
    data = request.get_json()
    try:
        result = calculate_beam_moment(
            rebar_list=data["rebar_list"],
            fc=float(data["fc"]),
            fy=float(data["fy"]),
            b=float(data["b"]),
            h=float(data["h"]),
            es=float(data.get("es", 200000)),
        )
        # Remove rebar_forces detail for JSON (too verbose)
        result.pop("rebar_forces", None)
        return jsonify({"status": "success", "result": result})
    except (KeyError, ValueError, TypeError) as e:
        return jsonify({"status": "error", "message": str(e)}), 400


@beam_bp.route("/shear", methods=["POST"])
def beam_shear():
    """Calculate beam shear capacity or required spacing."""
    data = request.get_json()
    try:
        mode = data.get("mode", "capacity")

        if mode == "capacity":
            vc_result = compute_concrete_shear_strength(
                fc=float(data["fc"]),
                b=float(data["b"]),
                d=float(data["d"]),
                lamda=float(data.get("lamda", 1.0)),
                vc_type=data.get("vc_type", "simple"),
            )
            vs_result = compute_steel_shear_strength(
                av=float(data["av"]),
                fyt=float(data["fyt"]),
                d=float(data["d"]),
                s=float(data["s"]),
            )
            vn = vc_result["vc"] + vs_result["vs"]
            phi = float(data.get("phi", 0.75))
            vu = phi * vn
            result = {
                **vc_result,
                **vs_result,
                "vn": round(vn, 2),
                "vn_kn": round(vn / 1000, 2),
                "vu": round(vu, 2),
                "vu_kn": round(vu / 1000, 2),
            }
        else:
            result = compute_shear_spacing(
                fc=float(data["fc"]),
                b=float(data["b"]),
                d=float(data["d"]),
                fyt=float(data["fyt"]),
                vu_required=float(data["vu_required"]),
                phi=float(data.get("phi", 0.75)),
                av=float(data["av"]),
                lamda=float(data.get("lamda", 1.0)),
                vc_type=data.get("vc_type", "simple"),
            )

        return jsonify({"status": "success", "result": result})
    except (KeyError, ValueError, TypeError) as e:
        return jsonify({"status": "error", "message": str(e)}), 400


@beam_bp.route("/torsion", methods=["POST"])
def beam_torsion():
    """Perform torsion design checks."""
    data = request.get_json()
    try:
        result = torsion_design(
            width=float(data["width"]),
            height=float(data["height"]),
            cover=float(data["cover"]),
            db=float(data["db"]),
            tf=float(data.get("tf", 0)),
            beff=float(data.get("beff", data["width"])),
            phi_torsion=float(data.get("phi_torsion", 0.75)),
            fc=float(data["fc"]),
            fy=float(data["fy"]),
            tu=float(data["tu"]),
            vc=float(data["vc"]),
            ds=float(data["ds"]),
            smax_shear=float(data["smax_shear"]),
            s_actual=float(data["s_actual"]),
            av=float(data["av"]),
            s=float(data["s"]),
        )
        return jsonify({"status": "success", "result": result})
    except (KeyError, ValueError, TypeError) as e:
        return jsonify({"status": "error", "message": str(e)}), 400


@beam_bp.route("/deflection", methods=["POST"])
def beam_deflection():
    """Calculate beam deflection."""
    data = request.get_json()
    try:
        result = deflection_computation(
            b=float(data["b"]),
            h=float(data["h"]),
            d=float(data["d"]),
            fc=float(data["fc"]),
            fy=float(data["fy"]),
            clearspan=float(data["clearspan"]),
            as_tension=float(data["as_tension"]),
            n_bars_comp=int(data.get("n_bars_comp", 0)),
            db_comp=float(data.get("db_comp", 0)),
            es=float(data.get("es", 200000)),
            ma=float(data["ma"]) if data.get("ma") else None,
            point_load=float(data.get("point_load", 0)),
            uniform_load=float(data.get("uniform_load", 0)),
            beam_type=data.get("beam_type", "simply_supported"),
            member_type=data.get("member_type", "floor"),
            sustained_duration=data.get("sustained_duration", "5_years_or_more"),
        )
        return jsonify({"status": "success", "result": result})
    except (KeyError, ValueError, TypeError) as e:
        return jsonify({"status": "error", "message": str(e)}), 400
