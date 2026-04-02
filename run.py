#!/usr/bin/env python3
# MIT License
# Copyright (c) Albert Pamonag Engineering Consultancy

"""
ConcreteDesignPy Web Application
==================================

Run this file to start the web interface:

    python run.py

Then open http://127.0.0.1:5000 in your browser.
"""

from concretedesignpy.webapp import create_app

app = create_app()

if __name__ == "__main__":
    print("=" * 50)
    print("  ConcreteDesignPy Web App")
    print("  http://127.0.0.1:5000")
    print("=" * 50)
    app.run(debug=True, host="127.0.0.1", port=5000)
