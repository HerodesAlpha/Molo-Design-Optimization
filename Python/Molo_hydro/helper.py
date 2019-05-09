import numpy as np

import helper
import nemoh
from meshmagick.MMviewer import MMViewer
from meshmagick.mesh import Mesh


class Force(object):
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

    def __init__(self, point=(0, 0, 0), value=(0, 0, 0)):
        assert len(point) == 3
        assert len(value) == 3
        # assert mode in ('relative', 'absolute')

        self._point = point
        self._value = value
        # self.name = str(name)
        # self.mode = mode


class show_force(object):

    def __init__(self, working_mesh, points, forces):
        self._points = points
        self._forces = forces
        self.backup = dict()
        self.mesh = working_mesh.copy()

    def show(self):

        vtk_polydata = self.mesh._vtk_polydata()
        self.viewer = MMViewer()
        self.viewer.add_polydata(vtk_polydata)

        self._forces = np.true_divide(self._forces, np.max(np.linalg.norm(self._forces, axis=1)[:, np.newaxis]))

        if self._points.shape[0] == self._forces.shape[0]:
            for i, p in enumerate(self._points):
                # self.viewer.add_point(p, color=[0, 0, 1])
                self.viewer.add_vector(p, self._forces[i] * 20, scale=3, color=[0, 0, 1])
                # print(self._forces[i])
        else:
            print(self._points.shape[0], self._forces.shape[0])
            exit()

        # Showing the viewer
        self.viewer.show_interactive()
        self.viewer.finalize()

        return



