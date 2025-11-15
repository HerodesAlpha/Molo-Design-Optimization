"""
Structural analysis utilities for converting pressures to forces.

This module provides functions for calculating forces from pressure distributions
on panel elements.
"""

__author__ = "Eivind Sonju"
__copyright__ = "Copyright (C) 2017-2019 Verbun AS. All rights reserved."
__version__ = "2.0"

from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from core.tool_box import PanelData

# TODO: Set up gravity section force


def p2f(pressure: np.ndarray, panel_data) -> np.ndarray:
    """
    Convert pressure distribution to force vectors on panels.
    
    Args:
        pressure: Pressure values for each panel [npanels]
        panel_data: PanelData object with panel geometry
        
    Returns:
        Force vectors [npanels, 3] in x, y, z directions
    """
    npanels = panel_data.ppanels.shape[0]
    f_normal = np.zeros((npanels), dtype=np.complex128)
    f = np.zeros((npanels,3), dtype=np.complex128)
    for i, panel in enumerate(panel_data.ppanels):
        f_normal[i] = pressure[i] * panel_data.ppanel_areas[i]
        for j in range(3):
            f[i,j] = -f_normal[i] * panel_data.ppanel_normals[i, j]
    return f