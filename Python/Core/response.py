__author__ = "Eivind Sonju"
__copyright__ = "Copyright (C) 2017-2019 Verbun AS. All rights reserved."
__version__ = "2.0"

import numpy as np
import scipy
import nemoh
import pickle
import tool_box as tb

nax = np.newaxis


class ResponseModel(object):
    def __init__(self, candidate, seastate = None):
        # super().__init__(settings)
        if not seastate is None:
            pass

        self._settings = candidate.settings
        self._loads = candidate.loads

        self._h = self._settings.job_data['floater']['Radial']['Heigth']
        self._d_rc = self._settings.job_data['floater']['Radial']['Column']['Diameter']
        self._d_cc = self._settings.job_data['floater']['Central column diameter']
        self._t_lf = self._settings.job_data['floater']['Radial']['Flange']['Lower']['Plate']['Thickness']
        self._t_uf = self._settings.job_data['floater']['Radial']['Flange']['Upper']['Plate']['Thickness']

        self._t_lfst = self._settings.job_data['floater']['Radial']['Flange']['Lower']['Stiffener']['Longitudinal'][
            'Thickness']
        self._h_lfst = self._settings.job_data['floater']['Radial']['Flange']['Lower']['Stiffener']['Longitudinal'][
            'Height']
        self._t_ufst = self._settings.job_data['floater']['Radial']['Flange']['Upper']['Stiffener']['Longitudinal'][
            'Thickness']
        self._h_ufst = self._settings.job_data['floater']['Radial']['Flange']['Upper']['Stiffener']['Longitudinal'][
            'Height']

        self._gaf = self._settings.job_data['floater']['Gap factor']
        self._nr = self._settings.job_data['floater']['Number of radials']
        self._nc = self._settings.job_data['floater']['Radial']['Number of columns']

        # Prepare discrete mass and hydro forces
        with open(self._settings.fio.data_io_dir.joinpath('unit_model.pkl'), 'rb') as f:
            self._unit_model = pickle.load(f)

        self._part_list = self._unit_model.get_parts_without_children()
        self._point_mass = np.asarray([part._inertias.mass_matrix_global for part in self._part_list])

        self._point_mass_centers = np.asarray([-part._inertias.reduction_point for part in self._part_list])
        self._panel_pressure_centers = self._loads._pd.ppanel_centers

        self._projected_panel_area = self._loads._an[nax, nax, :, :]

        # *****************************************************************************
        #                               S T A T I C
        # *****************************************************************************
        # w = self._loads.w
        self.w2 = self._loads.w ** 2
        # Gravity
        self._point_mass_gravity_force = self._point_mass[:, 2, 2][:, nax] * np.asarray(
                [0, 0, self._settings.gravity])[nax, :]  # Use m33

        # Hydro static / Buoyancy
        self._panel_pressure_buoyancy_force = self._loads._force['Buoyancy']

        # *****************************************************************************
        #                   D Y N A M I C  E Q U I L I B R I U M
        # *****************************************************************************

        # Set up RAOs
        self._rao = np.zeros([self._loads.nw, self._loads._nbeta, 6], dtype=complex)
        for ibeta in range(self._loads._nbeta):
            self._rao[:, ibeta, :] = self.calc_rao(ibeta)
        self.calc_dynamic_equilibrium()

    def calc_dynamic_equilibrium(self):
        # -----------------------------------------------------------------------------
        # Right hand side
        # -----------------------------------------------------------------------------

        # Froude-Krylof
        self._panel_pressure_froude_krylof_force = self._loads._force['Froude-Krylof']

        # Diffraction
        self._panel_pressure_diffraction_force = self._loads._force['Diffraction']

        # -----------------------------------------------------------------------------
        # Left hand side
        # -----------------------------------------------------------------------------

        # Radiation

        am = self._loads._added_mass
        rd = self._loads._radiation_damping
        # The RAO for each DOF is multiplied with each RAO dependent panel force (x,y,z)
        self._panel_added_mass_force_all_dof = -am[:, nax, :, :, :] * self._rao[:, :, :, nax,
                                                                      nax] * self.w2[:, nax,
                                                                             nax, nax,
                                                                             nax]
        self._panel_radiation_damping_force_all_dof = rd[:, nax, :, :, :] * self._rao[:, :, :, nax,
                                                                            nax] * self._loads.w[:,
                                                                                   nax,
                                                                                   nax, nax,
                                                                                   nax] * np.complex(0, 1)

        self._panel_added_mass_force = np.sum(self._panel_added_mass_force_all_dof, axis=2)
        self._panel_radiation_damping_force = np.sum(self._panel_radiation_damping_force_all_dof, axis=2)

        # RAO transformation matrix
        self._rao_tra_mat = np.zeros([self._loads._nw, self._loads._nbeta, 4, 4], dtype=complex)
        self._panel_pos = np.zeros([self._loads._nw, self._loads._nbeta, self._loads.pd.npanels, 3], dtype=complex)
        for ifreq in range(self._loads._nw):  # TODO: Vectorize
            for ibeta in range(self._loads._nbeta):
                rot = self._rao[ifreq, ibeta, 3:6]
                tra = self._rao[ifreq, ibeta, 0:3]
                self._rao_tra_mat[ifreq, ibeta, :, :] = tb.transformation_matrix(rot, tra)
        self._rao_transf_mat = self._rao_tra_mat[:, :, nax, :, :]

        # Stiffness
        # Gen dynamic position of panels and calc hydro static pressure
        # Append a 1 to the 3 dof vector to correspond with 4x4 tra_mat
        self._ppc = np.append(self._panel_pressure_centers, np.ones((self._loads.pd.npanels, 1)), 1)
        # Modify for broadcasting, add artificial dim to use matmul on stack of matrices
        self._ppc = self._ppc[nax, nax, :, :, nax]
        # Perform matmul and remove artificial dim and append 1. This code is fast ...
        self._panel_pos[:, :, :, :] = np.squeeze(np.matmul(self._rao_transf_mat, self._ppc), axis=4)[:, :, :, 0:3]
        self._panel_pos -= self._loads.pd.ppanel_centers[nax, nax, :, :]  # Subtract mean position
        self._panel_diff_buoyancy_pressure = (self._loads.rho_sw * abs(self._loads.gravity)) * self._panel_pos[:, :, :,
                                                                                               2]  # Change in pressure
        self._panel_diff_buoyancy_force = (-self._panel_diff_buoyancy_pressure[:, :, :,
                                            nax] * self._projected_panel_area)

        # Inertia
        # Get dynamic acceleration of part masses and calc inertia force
        #
        # self._pmc = np.append(self._point_mass_centers, np.ones((self._point_mass_centers.shape[0], 1)), 1)
        # self._pmc = self._pmc[nax, nax, :, :, nax]
        # self._dynamic_point_mass_pos = np.squeeze(np.matmul(self._rao_transf_mat, self._pmc), axis=4)[:, :, :, 0:3]
        # self._dynamic_point_mass_pos -= self._point_mass_centers[nax, nax, :]  # Subtract mean position
        #
        # self._part_dynacc = self._dynamic_point_mass_pos * -w2[:, nax, nax, nax]
        # self._point_mass_dynamic_inertia_force = self._part_dynacc * self._point_mass[nax, nax, :,
        #                                                              nax, 0, 0]

        self._point_mass_dynamic_inertia_force = np.squeeze(
                self._point_mass[nax, nax, :, :, :] @ self._rao[:, :, nax, :, nax],
                axis=4)
        # self._point_mass_dynamic_inertia_force *= -w2[:, nax, nax, nax]

    def calc_rao(self, idir):
        fe = self._loads._fe[:, idir, :]
        m = self._loads._m
        ma = self._loads._ma
        c_rad = self._loads._c_radiation * self._settings.radiaton_damping_factor
        k = self._loads._k
        vw = self._loads._w
        container = np.zeros([len(vw), 6], dtype=complex)
        for i, w in enumerate(vw):
            this_ma = ma[i, :, :]
            this_c = c_rad[i, :, :]
            denom = np.asarray(-w ** 2 * (m + this_ma) + 1j * w * this_c + k, dtype=complex)
            daf = scipy.linalg.inv(denom)
            x = daf @ fe[i, :]
            container[i, :] = x
        return container

    def get_section_index(self, section_point, section_normal):

        def mask1(coordinates):
            vec = coordinates - section_point
            dot = np.dot(vec,
                         section_normal)  # dot product i positive for coordinates on the positive side of the plane
            return dot >= 0

        return mask1(self._point_mass_centers), mask1(self._panel_pressure_centers)

    def get_flange_panel_index(self):
        def printv(string, do_print=None):
            if do_print == None:
                do_print = False
            if do_print:
                print(string)

        class BreakIt(Exception):
            pass

        def mask2(coordinates):
            n = coordinates.shape[0]
            da = (1 + self._gaf) * self._d_rc
            dtheta = 2 * np.pi / 3
            theta = [i * dtheta for i in range(3)]
            index_1 = np.zeros((n), dtype=bool)
            for i in range(n):
                found_inside = False
                is_flange_element = False
                xp = coordinates[i, 0]
                yp = coordinates[i, 1]
                xc = yc = 0

                if coordinates[i, 2] <= (np.min(
                        coordinates[:, 2]) + self._t_lf + 0.01):  # Flange element coordinate is at lowest position
                    # Next, find elements not inside the cylinders
                    try:
                        for irad in range(3):
                            dxc = np.cos(theta[irad]) * da
                            dyc = np.sin(theta[irad]) * da
                            for icol in range(self._nc):
                                xc = dxc * (icol + 1)
                                yc = dyc * (icol + 1)
                                if np.sqrt((xp - xc) ** 2 + (yp - yc) ** 2) < self._d_rc / 2:
                                    found_inside = True
                                    raise BreakIt
                    except BreakIt:
                        pass
                    if np.sqrt((xp) ** 2 + (yp) ** 2) < self._d_cc / 2:
                        found_inside = True

                    if not found_inside:
                        printv('   Is flange outside cylinders')
                        # print(i)
                        index_1[i] = True

            # Next find elements between cylinders
            # Hardcoded first radial, first bay
            # TODO: Make generic
            vec_a = coordinates - [0, 0, 0]  # First bay starts at origin
            dot_a = np.dot(vec_a, [1, 0, 0])  # dot product i positive for coordinates on the positive side of the plane
            index_2 = dot_a > 0

            vec_b = coordinates - [da, 0, 0]  # First bay ends at first radial column center
            dot_b = np.dot(vec_b, [1, 0, 0])  # dot product i positive for coordinates on the positive side of the plane
            index_3 = dot_b < 0  #

            return np.logical_and(np.logical_and(index_1, index_2), index_3)  # Return union of the three indexes

        return mask2(self._point_mass_centers), mask2(self._panel_pressure_centers)

    def sum_forces(self, index, forces, coordinates, moment_ref_point):

        vec = coordinates - moment_ref_point

        if 1 in index:  # At least one item is on the considered side of the section surface
            if forces.ndim == 4:  # Dynamic [freq, dir, panel, f]
                section_forces = forces[:, :, index, :]
                section_moments = np.cross(vec[nax, nax, index, :],
                                           section_forces)  # Calculate moment about section
                # Concatenate along 4th dimension contaning [fx, fy, fz] and [mx, mz, mz]
                # Then sum along 3rd dimension holding the panels or point mass indices
                return np.sum(np.concatenate((section_forces, section_moments), axis=3), axis=2)
            elif forces.ndim == 2:  # Static [panel, f]
                section_forces = forces[index, :]
                section_moments = np.cross(vec[index, :], section_forces)  # Calculate moment about section
                # Concatenate along 2nd dimension contaning [fx, fy, fz] and [mx, mz, mz]
                # Then sum along 1st dimension holding the panels or point mass indices
                return np.sum(np.concatenate((section_forces, section_moments), axis=1), axis=0)
        else:
            return np.zeros(6)

    def assemble_forces(self, imass, ipanel, moment_ref_point):

        # ---------------------------------------------------------------------------
        # STATIC FORCES
        # ---------------------------------------------------------------------------
        f_gravity = self.sum_forces(imass, self._point_mass_gravity_force, self._point_mass_centers, moment_ref_point)
        f_bouyancy = self.sum_forces(ipanel, self._panel_pressure_buoyancy_force, self._panel_pressure_centers,
                                     moment_ref_point)
        #
        # ---------------------------------------------------------------------------
        # DYNAMIC EXCITATION FORCES
        # ---------------------------------------------------------------------------
        f_fk = self.sum_forces(ipanel, self._panel_pressure_froude_krylof_force, self._panel_pressure_centers,
                               moment_ref_point)
        f_diff = self.sum_forces(ipanel, self._panel_pressure_diffraction_force, self._panel_pressure_centers,
                                 moment_ref_point)
        #
        # ---------------------------------------------------------------------------
        # DYNAMIC REACTION FORCES
        # ---------------------------------------------------------------------------
        # f_inertia = self.sum_forces(imass, self._point_mass_dynamic_inertia_force, self._point_mass_centers,
        #                             moment_ref_point)
        f_inertia = np.sum(self._point_mass_dynamic_inertia_force[:, :, imass, :], axis=2)
        f_inertia *= -self.w2[:, nax, nax, ]

        f_dz_s = self.sum_forces(ipanel, self._panel_diff_buoyancy_force, self._panel_pressure_centers,
                                 moment_ref_point)
        f_added_mass = self.sum_forces(ipanel, self._panel_added_mass_force, self._panel_pressure_centers,
                                       moment_ref_point)

        f_radiation_damping = self.sum_forces(ipanel, self._panel_radiation_damping_force, self._panel_pressure_centers,
                                              moment_ref_point)
        #
        force_out = dict()
        force_out['Static'] = dict()
        force_out['Dynamic'] = dict()

        force_out['Static']['Gravity'] = f_gravity
        force_out['Static']['Buoyancy'] = f_bouyancy
        force_out['Static']['SUM'] = f_gravity + f_bouyancy
        force_out['Dynamic']['Froude-Krylof'] = f_fk
        force_out['Dynamic']['Diffraction'] = f_diff
        force_out['Dynamic']['Mass'] = f_inertia
        force_out['Dynamic']['Added mass'] = f_added_mass
        force_out['Dynamic']['Radiation damping'] = f_radiation_damping
        force_out['Dynamic']['Buoyancy'] = f_dz_s
        force_out['Dynamic']['SUM'] = f_fk + f_diff - (f_added_mass + f_radiation_damping + f_dz_s + f_inertia)

        return force_out

    def force_comp(f):
        if f.ndim == 3:
            fx = f[:, :, 0]
            fy = f[:, :, 1]
            fz = f[:, :, 2]
            mx = f[:, :, 3]
            my = f[:, :, 4]
            mz = f[:, :, 5]
        elif f.ndim == 1:
            fx = f[0]
            fy = f[1]
            fz = f[2]
            mx = f[3]
            my = f[4]
            mz = f[5]

    def radial_cross_section(self):
        fdi = self._settings.job_data['floater']

        h = fdi['Radial height']
        d_rc = fdi['Radial column diameter']
        d_cc = fdi['Central column diameter']
        w = d_rc
        t_lf = fdi['Lower flange thickness']
        t_uf = fdi['Upper flange thickness']

        t_lfst = fdi['Lower flange stiffener thickness']
        h_lfst = fdi['Lower flange stiffener height']
        t_ufst = fdi['Upper flange stiffener thickness']
        h_ufst = fdi['Upper flange stiffener height']

    @property
    def rao(self):
        return self._rao
