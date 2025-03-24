# MIT License
#
# Copyright (c) 2025 ...
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

def eRF(exposure: str, fiber_type: str) -> float:
    """
    Calculates the Environmental Reduction Factor (eRF) for Fiber-Reinforced Polymer (FRP)
    materials based on the exposure condition and fiber type, following ACI 440 guidelines.

    Parameters
    ----------
    exposure : str
        Exposure environment. Accepted values: 'interior', 'exterior', or 'aggressive'.
    fiber_type : str
        Type of fiber used in the FRP. Accepted values: 'carbon', 'glass', or 'aramid'.

    Returns
    -------
    float
        The environmental reduction factor for the specified exposure and fiber type.

    Raises
    ------
    ValueError
        If the exposure condition or fiber type is not recognized.
    """

    data = {
        "interior": {
            "carbon": 0.95,
            "glass": 0.75,
            "aramid": 0.85
        },
        "exterior": {
            "carbon": 0.85,
            "glass": 0.65,
            "aramid": 0.75
        },
        "aggressive": {
            "carbon": 0.85,
            "glass": 0.50,
            "aramid": 0.70
        }
    }

    # Guard clauses
    if exposure not in data:
        raise ValueError(
            f"Invalid exposure '{exposure}'. "
            f"Supported exposures: {list(data.keys())}"
        )
    if fiber_type not in data[exposure]:
        raise ValueError(
            f"Invalid fiber type '{fiber_type}' for exposure '{exposure}'. "
            f"Supported fiber types: {list(data[exposure].keys())}"
        )

    return data[exposure][fiber_type]
