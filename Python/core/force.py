"""
Force visualization and representation module.

This module provides classes for representing forces and visualizing
them on meshes using VTK.
"""

__author__ = "Eivind Sonju"
__copyright__ = "Copyright (C) 2017-2019 Verbun AS. All rights reserved."
__version__ = "2.0"

from typing import Tuple, Union

import numpy as np

from core.meshmagick.MMviewer import MMViewer


class Force:
    """Class to handle a linear force (resultant).

    Parameters
    ----------
    point : array_like, optional
        The force application point. Default is (0, 0, 0)
    value : array_like, optional
        The value of the force. Default is (0, 0, 0)
    name : str, optional
        The name of the force
    mode : str, ['relative', 'absolute']
        The mode of the force whether its direction follows the body motion ('relative') or remains fixed with
        respect to the reference frame ('absolute'). Default is 'relative'.

    """

    def __init__(
        self,
        point: Union[Tuple[float, float, float], np.ndarray] = (0, 0, 0),
        value: Union[Tuple[float, float, float], np.ndarray] = (0, 0, 0)
    ) -> None:
        """
        Initialize force with application point and value.
        
        Args:
            point: Force application point [x, y, z]
            value: Force vector [fx, fy, fz]
        """
        if len(point) != 3:
            raise ValueError(f"point must have length 3, got {len(point)}")
        if len(value) != 3:
            raise ValueError(f"value must have length 3, got {len(value)}")
        
        self._point = np.asarray(point, dtype=np.float64)
        self._value = np.asarray(value, dtype=np.float64)


class show_force:
    """
    Class for visualizing forces on a mesh using VTK.
    """
    
    def __init__(self, working_mesh, points: np.ndarray, forces: np.ndarray) -> None:
        """
        Initialize force visualization.
        
        Args:
            working_mesh: Mesh to display forces on
            points: Array of force application points [n, 3]
            forces: Array of force vectors [n, 3]
        """
        self._points = points
        self._forces = forces
        self.backup = dict()
        self.mesh = working_mesh.copy()

    def show(self) -> None:
        """Display forces on the mesh using VTK viewer."""
        vtk_polydata = self.mesh._vtk_polydata()
        self.viewer = MMViewer()
        self.viewer.add_polydata(vtk_polydata)

        # Normalize forces for visualization
        max_force = np.max(np.linalg.norm(self._forces, axis=1))
        if max_force > 0:
            self._forces = self._forces / max_force

        if self._points.shape[0] == self._forces.shape[0]:
            for i, p in enumerate(self._points):
                self.viewer.add_vector(p, self._forces[i] * 20, scale=3, color=[0, 0, 1])
        else:
            raise ValueError(
                f"Points and forces shape mismatch: {self._points.shape[0]} points, "
                f"{self._forces.shape[0]} forces"
            )

        # Showing the viewer
        self.viewer.show()
        self.viewer.finalize()

        return



