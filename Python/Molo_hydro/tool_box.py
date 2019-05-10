import os
import pickle
import subprocess
from copy import deepcopy
from math import cos, sin

import numpy as np
from numpy import pi, sin, cos


class TotalMassMatrixClass(object):
    def __init__(self):
        self._cog = np.zeros((3), dtype='float')
        self._point = np.zeros((3), dtype='float')
        self._mass_matrix_local = np.zeros((6, 6), dtype='float')
        self._mass_matrix_global = np.zeros((6, 6), dtype='float')

    @property
    def cog(self):
        """The position of the center of gravity"""
        return self._cog

    # @cog.setter
    # def cog(self, point):
    #     """The position of the center of gravity"""
    #     self._cog = np.asarray(point, dtype=np.float)

    @property
    def mass(self):
        """The mass of the body"""
        return self._mass_matrix_local[0, 0]

    @mass.setter
    def mass(self, val):
        self._mass_matrix_local[:3, :3] = np.eye(3) * val

    @property
    def mass_matrix_local(self):
        return self._mass_matrix_local

    @property
    def ixx(self):
        return self._mass_matrix_local[3][3]

    @ixx.setter
    def ixx(self, val):
        self.assign_val_to_mass_matrix_local(3, 3, val)

    @property
    def iyy(self):
        return self._mass_matrix_local[4, 4]

    @iyy.setter
    def iyy(self, val):
        self.assign_val_to_mass_matrix_local(4, 4, val)

    @property
    def izz(self):
        return self._mass_matrix_local[5][5]

    @izz.setter
    def izz(self, val):
        self.assign_val_to_mass_matrix_local(5, 5, val)

    @property
    def inertia_matrix_global(self):
        return self._mass_matrix_global[3:, 3:]

    @property
    def mass_matrix_global(self):
        self.update_mass_matrix_global()
        return self._mass_matrix_global

    # Local matrix
    @property
    def vector_matrix_local(self):
        return self._mass_matrix_local / self.mass

    #

    def add_mass_matrix_local(self, mat):
        self._mass_matrix_local += mat

    def assign_val_to_mass_matrix_local(self, i, j, val):
        self._mass_matrix_local[i, j] = val
        self.update_mass_matrix_global()

    def rotate_z_inertia_mass_matrix_local(self, angle):
        rm = rotation_matrix([0, 0, angle])
        rm_t = np.transpose(rm)
        self._mass_matrix_local[3:, :3] = rm * self._mass_matrix_local[3:, :3] * rm_t
        self._mass_matrix_local[:3, 3:] = rm * self._mass_matrix_local[:3, 3:] * rm_t
        self._mass_matrix_local[3:, 3:] = rm * self._mass_matrix_local[3:, 3:] * rm_t
        self.update_mass_matrix_global()

    # Global matrix
    @property
    def vector_matrix_global(self):
        return self._mass_matrix_global / self.mass

    def rotate_z_inertia_mass_matrix_global(self, angle):
        rm = rotation_matrix([0, 0, angle])
        self._mass_matrix_global[3:, 3:] = rm * self._mass_matrix_global[3:, 3:] * np.transpose(rm)

    def reset(self):
        self.__init__()

    @property
    def reduction_point(self):
        """
        The reduction point of the inertia matrix

        Returns
        -------
        ndarray
        """
        return self._point

    @reduction_point.setter
    def reduction_point(self, point):
        """Set the reduction point"""
        assert len(point) == 3
        self._point = np.asarray(point, dtype=np.float)
        self.update_mass_matrix_global()

    def update_mass_matrix_global(self):
        self._mass_matrix_global = self._mass_matrix_local + self._huygens_transport() * self.mass

    @property
    def at_cog(self):
        """Returns a new inertia object that is expressed at cog.

        It makes a copy of itself.

        Returns
        -------
        ndarray
        """
        inertia = deepcopy(self)
        inertia.shift_at_cog()
        return inertia

    def shift_at_cog(self):
        """Shift the inertia matrix internally at cog.

        The reduction point is then cog.
        """
        self._mass_matrix_local -= self._huygens_transport() * self.mass
        self._point = self._cog

    def is_at_cog(self):
        """Returns whether the object is expressed at cog

        Returns
        -------
        bool
        """
        return np.all(self._point == self._cog)

    def _huygens_transport(self):
        p_g = self._cog - self._point
        x = p_g[0];
        y = p_g[1];
        z = p_g[2]
        # print('x={} y={} z={}'.format(x,y,z))
        A = np.asarray([[0, -z, y],
                        [z, 0, -x],
                        [-y, x, 0]], dtype=np.float)
        AT = np.transpose(A)
        m = np.zeros((6, 6), dtype=np.float)
        m[3:, :3] = A
        m[:3, 3:] = AT
        m[3:, 3:] = np.matmul(A, AT)
        # print(m[3:,3:])
        return m


def hollow_right_circular_cylinder(tmm, int_radius, ext_radius, length, density=1.):
    """Get the inertia of a hollow right circular cylinder

    Returns
    -------
    TotalMassMatrixClass instance
    """

    vol = pi * length * (ext_radius ** 2 - int_radius ** 2)
    tmm.mass = density * vol

    R2r2 = (ext_radius ** 2 + int_radius ** 2)
    tmm.ixx = tmm.iyy = ((3 * R2r2 + length ** 2) / 12.) * tmm.mass
    tmm.izz = (R2r2 / 2.) * tmm.mass


def rectangular_prism(tmm, a, b, h, density=1.):
    """Get the inertia of a rectangular prism

    Returns
    -------
    TotalMassMatrixClass instance

    Note
    ----
    * a is along x
    * b is along y
    * h is along z
    """
    vol = a * b * h
    tmm.mass = density * vol

    tmm.ixx = tmm.mass * (b ** 2 + h ** 2) / 12.
    tmm.iyy = tmm.mass * (a ** 2 + h ** 2) / 12.
    tmm.izz = tmm.mass * (a ** 2 + b ** 2) / 12.


def trig(angle):
    # r = radians(angle)
    r = angle
    return cos(r), sin(r)


def transformation_matrix(rotation, translation=None, scale=None):
    if translation is None:
        translation = [0, 0, 0]
    if scale is None:
        scale = [1, 1, 1]
    xC, xS = trig(rotation[0])
    yC, yS = trig(rotation[1])
    zC, zS = trig(rotation[2])
    dX = translation[0]
    dY = translation[1]
    dZ = translation[2]
    sX = scale[0]
    sY = scale[1]
    sZ = scale[2]
    Scale_matrix = np.array([[sX, 0, 0, 0],
                             [0, sY, 0, 0],
                             [0, 0, sZ, 0],
                             [0, 0, 0, 1]])
    Translate_matrix = np.array([[1, 0, 0, dX],
                                 [0, 1, 0, dY],
                                 [0, 0, 1, dZ],
                                 [0, 0, 0, 1]])
    Rotate_X_matrix = np.array([[1, 0, 0, 0],
                                [0, xC, -xS, 0],
                                [0, xS, xC, 0],
                                [0, 0, 0, 1]])
    Rotate_Y_matrix = np.array([[yC, 0, yS, 0],
                                [0, 1, 0, 0],
                                [-yS, 0, yC, 0],
                                [0, 0, 0, 1]])
    Rotate_Z_matrix = np.array([[zC, -zS, 0, 0],
                                [zS, zC, 0, 0],
                                [0, 0, 1, 0],
                                [0, 0, 0, 1]])
    return np.dot(Rotate_Z_matrix,
                  np.dot(Rotate_Y_matrix, np.dot(Rotate_X_matrix, np.dot(Translate_matrix, Scale_matrix))))


def rotation_matrix(rotation):
    xC, xS = trig(rotation[0])
    yC, yS = trig(rotation[1])
    zC, zS = trig(rotation[2])
    Rotate_X_matrix = np.array([[1, 0, 0],
                                [0, xC, -xS],
                                [0, xS, xC]])
    Rotate_Y_matrix = np.array([[yC, 0, yS],
                                [0, 1, 0],
                                [-yS, 0, yC], ])
    Rotate_Z_matrix = np.array([[zC, -zS, 0],
                                [zS, zC, 0],
                                [0, 0, 1]])
    return np.dot(Rotate_Z_matrix, np.dot(Rotate_Y_matrix, Rotate_X_matrix))


def write_gmsh(fio, floater_model, dens_t=1, dens_quarter_cirlce=4, dens_cylinder_height=14):
    with open(r'.\templates\MOLO_{}c.geo.template'.format(floater_model.nc), 'r') as file:
        filedata = file.read()

    # Replace the target string
    filedata = filedata.replace('#dia_rc#', '{}'.format(floater_model.dia_rc))
    filedata = filedata.replace('#dia_hc#', '{}'.format(floater_model.dia_hc))
    filedata = filedata.replace('#gap#', '{}'.format(floater_model.gap))
    filedata = filedata.replace('#hgt#', '{}'.format(floater_model.hgt))
    filedata = filedata.replace('#t_lf#', '{}'.format(floater_model.t_lf))

    # Set mesh density
    dens2 = dens_t + 1  # Thickness
    dens5 = dens_quarter_cirlce + 1  #
    dens9 = 2 * dens_quarter_cirlce + 1  #
    dens15 = dens_cylinder_height + 1  # Column height
    filedata = filedata.replace('#dens2#', '{}'.format(dens2))
    filedata = filedata.replace('#dens5#', '{}'.format(dens5))
    filedata = filedata.replace('#dens9#', '{}'.format(dens9))
    filedata = filedata.replace('#dens15#', '{}'.format(dens15))

    gmsh_geo_file = fio.gmsh_dir.joinpath('MOLO_{}c.geo'.format(floater_model.nc))
    gmsh_msh_file = fio.gmsh_dir.joinpath('MOLO_{}c.msh'.format(floater_model.nc))

    with open(gmsh_geo_file, 'w') as file:
        file.write(filedata)
    print('{}'.format(gmsh_geo_file))
    try:
        a = subprocess.check_output(
            [fio.gmsh_exe, '-2', '{}'.format(gmsh_geo_file), '-save_all', '-format', 'msh2', '-o',
             '{}'.format(gmsh_msh_file)])
    except:
        print(a)
        exit()
    finally:
        return gmsh_msh_file


def save_M_and_K(wdir, M, MMK):
    # From meshmagic 3x3 to 6x6
    K = np.zeros((6, 6), dtype='float')
    K[2, 2] = MMK[0, 0]
    K[2, 3] = MMK[0, 1]
    K[3, 2] = MMK[0, 1]
    K[2, 4] = MMK[0, 2]
    K[4, 2] = MMK[0, 2]
    K[3, 3] = MMK[1, 1]
    K[3, 4] = MMK[1, 2]
    K[4, 3] = MMK[1, 2]
    K[4, 4] = MMK[2, 2]

    pickle.dump(M, open(os.path.join(wdir, 'M.pkl'), "wb"))
    pickle.dump(K, open(os.path.join(wdir, 'K.pkl'), "wb"))


def load_M_and_K(wdir):
    M = pickle.load(open(os.path.join(wdir, 'M.pkl'), 'rb'))
    K = pickle.load(open(os.path.join(wdir, 'K.pkl'), 'rb'))

    return M, K


def get_dof_index(sel_dof, dof):
    dof_index = 0
    if dof[sel_dof - 1]:
        for i in dof:
            if i:
                dof_index += 1
                if dof_index == sel_dof:
                    break


    else:
        print('No such DOF')

    return dof_index - 1


#
def get_dir_index(sel_dir, dir):
    return sel_dir


def P2R(radii, angles):
    return radii * np.exp(1j * angles)


def R2P(x):
    return abs(x), np.angle(x)

