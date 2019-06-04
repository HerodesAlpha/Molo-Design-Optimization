__author__ = "Eivind Sonju"
__copyright__ = "Copyright (C) 2017-2019 Verbun AS. All rights reserved."
__version__ = "2.0"

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


def write_gmsh(fio, floater_model, dens_t=1, dens_quarter_cirlce=4,
               dens_cylinder_height=14):  # TODO: Delete this function
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


def msh_file(fio, settings, bool_thin):
    if bool_thin:
        str_thin = '_thin'
    else:
        str_thin = ''

    # str_thin = '_thin' if settings.job_data['analysis']['simulations']['default']['calculation'][
    #     'use_dipoles_implementation'] else ''

    if settings.mesh_name == None:
        settings.mesh_name = 'MOLO_{}c'.format(settings.job_data['floater']['Number of radial columns'])

    this_mesh_name = '{}{}'.format(settings.mesh_name, str_thin)

    with open(fio.templates_dir.joinpath('{}.geo.template'.format(this_mesh_name)), 'r') as file:
        filedata = file.read()

    # Replace the target string
    filedata = filedata.replace('#dia_rc#', '{}'.format(settings.job_data['floater']['Radial column diameter']))
    filedata = filedata.replace('#dia_hc#', '{}'.format(settings.job_data['floater']['Central column diameter']))
    filedata = filedata.replace('#gap#', '{}'.format(settings.job_data['floater']['Gap factor']))
    filedata = filedata.replace('#t_lf#', '{}'.format(settings.job_data['floater']['Lower flange thickness']))

    if bool_thin:
        filedata = filedata.replace('#hgt#', '{}'.format(settings.job_data['floater']['Draught']))
    else:
        filedata = filedata.replace('#hgt#', '{}'.format(settings.job_data['floater']['Radial height']))

    # Set mesh density thin templates
    filedata = filedata.replace('#nel#', '{}'.format(16))

    # Set mesh density old templates
    dens_t = 1
    dens_quarter_cirlce = 4
    dens_cylinder_height = 14
    dens2 = dens_t + 1  # Thickness
    dens5 = dens_quarter_cirlce + 1  #
    dens9 = 2 * dens_quarter_cirlce + 1  #
    dens15 = dens_cylinder_height + 1  # Column height
    filedata = filedata.replace('#dens2#', '{}'.format(dens2))
    filedata = filedata.replace('#dens5#', '{}'.format(dens5))
    filedata = filedata.replace('#dens9#', '{}'.format(dens9))
    filedata = filedata.replace('#dens15#', '{}'.format(dens15))

    # Write gmsh geo file to analysis directory
    gmsh_geo_file = fio.gmsh_dir.joinpath('{}.geo'.format(this_mesh_name))
    with open(gmsh_geo_file, 'w') as file:
        file.write(filedata)

    # Create mesh
    gmsh_msh_file = fio.gmsh_dir.joinpath('{}.msh'.format(this_mesh_name))
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


class PanelData(object):
    def __init__(self, vertices, faces):
        self._ppoints = vertices
        self._ppanels = faces
        self._npoints = len(vertices)
        self._npanel = len(faces)
        _a1 = np.linalg.norm(np.cross(self._ppoints[self._ppanels[:, 1]] - self._ppoints[self._ppanels[:, 0]],
                                      self._ppoints[self._ppanels[:, 2]] - self._ppoints[self._ppanels[:, 0]]),
                             axis=1) * 0.5
        _a2 = np.linalg.norm(np.cross(self._ppoints[self._ppanels[:, 3]] - self._ppoints[self._ppanels[:, 0]],
                                      self._ppoints[self._ppanels[:, 2]] - self._ppoints[self._ppanels[:, 0]]),
                             axis=1) * 0.5
        self._ppanel_areas = _a1 + _a2

        self._c1 = np.sum(self._ppoints[self._ppanels[:, :3]], axis=1) / 3.
        self._c2 = (np.sum(self._ppoints[self._ppanels[:, 2:4]], axis=1) + self._ppoints[self._ppanels[:, 0]]) / 3.

        self._ppanel_centers = (np.array(([_a1, ] * 3)).T * self._c1 + np.array(([_a2, ] * 3)).T * self._c2)
        self._ppanel_centers /= np.array(([self._ppanel_areas, ] * 3)).T

        self._ppanel_cross = np.cross(self._ppoints[self._ppanels[:, 2]] - self._ppoints[self._ppanels[:, 0]],
                                      self._ppoints[self._ppanels[:, 3]] - self._ppoints[self._ppanels[:, 1]])
        self._norm = np.linalg.norm(self._ppanel_cross, axis=1)
        self._ppanel_normals = np.true_divide(self._ppanel_cross, self._norm[:, np.newaxis])

    @property
    def ppoints(self):
        return self._ppoints

    @property
    def ppanels(self):
        return self._ppanels

    @property
    def ppanel_normals(self):
        return self._ppanel_normals

    @property
    def ppanel_areas(self):
        return self._ppanel_areas

    @property
    def ppanel_centers(self):
        return self._ppanel_centers

    @property
    def npoints(self):
        return self._ppoints.shape[0]

    @property
    def npanels(self):
        return self._ppanels.shape[0]


def prepare_dipol_mesh(vertices, faces, settings):
    def printv(string):
        if 0:
            print(string)

    class BreakIt(Exception):
        pass

    vtol = 0.01
    drc = settings.job_data['floater']['Radial column diameter']
    dcc = settings.job_data['floater']['Central column diameter']
    gaf = settings.job_data['floater']['Gap factor']
    ncol = settings.job_data['floater']['Number of radial columns']

    da = (1 + gaf) * drc
    dtheta = 2 * np.pi / 3
    theta = [i * dtheta for i in range(3)]
    limit = drc / 2 * (1 + vtol)
    pd = PanelData(vertices, faces)

    is_flange_element = False
    dipol_index = []
    not_dipol_index = []
    flipped = []
    #
    for i in range(pd.npanels):
        found_inside = False
        is_flange_element = False
        xp = pd.ppanel_centers[i][0]
        yp = pd.ppanel_centers[i][1]
        xc = yc = 0

        if abs(pd.ppanel_normals[i][2]) == 1:  # Flange element
            printv('Process flange element {}'.format(i))
            is_flange_element = True
            # Check normal and flip if positive up
            if pd.ppanel_normals[i][2] == 1:
                flipped.append(i)
                printv('   Flip normal')

                faces[i] = faces[i][::-1]
            # Next, find dipols, i.e. elements not inside the cylinders
            # Assume dipol if not found inside cylinders

            try:
                for irad in range(3):
                    dxc = np.cos(theta[irad]) * da
                    dyc = np.sin(theta[irad]) * da
                    for icol in range(ncol):
                        xc = dxc * (icol + 1)
                        yc = dyc * (icol + 1)
                        if np.sqrt((xp - xc) ** 2 + (yp - yc) ** 2) < drc / 2:
                            found_inside = True
                            raise BreakIt
            except BreakIt:
                pass
            if np.sqrt((xp) ** 2 + (yp) ** 2) < dcc / 2:
                found_inside = True

            if not found_inside:
                printv('   Is dipol')
                # print(i)
                dipol_index.append(i)
            else:
                not_dipol_index.append(i)
        else:  # Cylinder element
            # Find the cylinder x,y to which the element belongs
            printv('Process cylinder element {}'.format(i))
            not_found = True
            try:
                for irad in range(3):
                    dxc = np.cos(theta[irad]) * da
                    dyc = np.sin(theta[irad]) * da
                    for icol in range(ncol):
                        xc = dxc * (icol + 1)
                        yc = dyc * (icol + 1)
                        dx = np.abs(xp - xc)
                        dy = np.abs(yp - yc)
                        if dx < limit and dy < limit:
                            printv('   Found on radial {}, column {}'.format(irad + 1, icol + 1))
                            not_found = False
                            raise BreakIt
            except BreakIt:
                pass
            if not_found:
                xc = yc = 0  # Check center column
                if np.sqrt((xp) ** 2 + (yp) ** 2) < dcc / 2 * (1 + vtol):
                    printv('   Found on center column'.format(i))
                    not_found = False
            if not_found:
                printv('   Not found on any column'.format(i))
                exit()

                # Check if normal is outwards from cylinder center
            printv('   xp = {: 6.2f}, yp = {: 6.2f}'.format(xp, yp))
            printv('   xc = {: 6.2f}, yc = {: 6.2f}'.format(xc, yc))
            vec = np.asarray([xp - xc, yp - yc])
            pn = pd.ppanel_normals[i][0:2]
            if np.dot(vec, pn) < 0:  # dot product i positive for coordinates on the positive side of the plane
                printv('   Flip normal')
                faces[i] = faces[i][::-1]
                flipped.append(i)


            not_dipol_index.append(i)

    if not (len(not_dipol_index) + len(dipol_index)) == len(faces):  # TODO: Fix dipol filter
        print(' prepare_dipol_mesh failed\n\tnot_dipol - {}\n\tdipol     - {}\n\ttotal     - {}'.format(
            len(not_dipol_index), len(dipol_index), len(faces)))
        exit()

    settings.job_data['analysis']['simulations']['default']['calculation'][
        'thin_panels'] = dipol_index  # TODO: Check if index must start with 1
    return vertices, faces, not_dipol_index, dipol_index
