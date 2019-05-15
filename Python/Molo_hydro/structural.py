__author__ = "Eivind Sonju"
__copyright__ = "Copyright (C) 2017-2019 Verbun AS. All rights reserved."
__version__ = "2.0"

import numpy as np

# TODO: Set up gravity section force

def p2f(pressure, panel_data):
    npanels = panel_data.ppanels.shape[0]
    f_normal = np.zeros((npanels), dtype=np.complex)
    f = np.zeros((npanels,3), dtype=np.complex)
    for i, panel in enumerate(panel_data.ppanels):
        f_normal[i] = pressure[i] * panel_data.ppanel_areas[i]
        for j in range(3):
            f[i,j] = -f_normal[i] * panel_data.ppanel_normals[i, j]
    return f