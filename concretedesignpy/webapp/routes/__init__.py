# MIT License
# Copyright (c) Albert Pamonag Engineering Consultancy

"""Route blueprints for the web application."""

from concretedesignpy.webapp.routes.main import main_bp
from concretedesignpy.webapp.routes.beam import beam_bp
from concretedesignpy.webapp.routes.column import column_bp
from concretedesignpy.webapp.routes.joint import joint_bp
from concretedesignpy.webapp.routes.mander import mander_bp
from concretedesignpy.webapp.routes.frp import frp_bp
from concretedesignpy.webapp.routes.section import section_bp

__all__ = [
    "main_bp", "beam_bp", "column_bp", "joint_bp",
    "mander_bp", "frp_bp", "section_bp",
]
