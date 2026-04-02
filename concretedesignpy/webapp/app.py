# MIT License
# Copyright (c) Albert Pamonag Engineering Consultancy

"""
Flask Application Factory
===========================

Creates and configures the Flask application with all
calculator route blueprints.
"""

from flask import Flask

from concretedesignpy.webapp.routes import (
    main_bp,
    beam_bp,
    column_bp,
    joint_bp,
    mander_bp,
    section_bp,
)


def create_app():
    """Create and configure the Flask application."""
    app = Flask(
        __name__,
        template_folder="templates",
        static_folder="static",
    )
    app.config["SECRET_KEY"] = "concretedesignpy-dev"

    app.register_blueprint(main_bp)
    app.register_blueprint(beam_bp, url_prefix="/api/beam")
    app.register_blueprint(column_bp, url_prefix="/api/column")
    app.register_blueprint(joint_bp, url_prefix="/api/joint")
    app.register_blueprint(mander_bp, url_prefix="/api/mander")
    app.register_blueprint(section_bp, url_prefix="/api/section")

    return app


def main():
    """Entry point for the console script."""
    app = create_app()
    app.run(debug=True, host="127.0.0.1", port=5000)
