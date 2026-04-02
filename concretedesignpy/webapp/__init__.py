# MIT License
# Copyright (c) Albert Pamonag Engineering Consultancy

"""
concretedesignpy.webapp
========================

Flask web application for concrete design calculations.
Provides a browser-based UI and JSON API for all calculators.
"""

from concretedesignpy.webapp.app import create_app

__all__ = ["create_app"]
