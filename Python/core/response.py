__author__ = "Eivind Sonju"
__copyright__ = "Copyright (C) 2017-2019 Verbun AS. All rights reserved."
__version__ = "2.0"

import numpy as np
import scipy
import core.nemoh
import pickle
import core.tool_box as tb

nax = np.newaxis


class Viscous_Damper(object):
    def __init__(self, coordinate, rho, cd, d, l):
        # def __init__(self, pos, rho, cd, d, l, w, sea_spectrum, rao):

        self.coordinate = coordinate
        self.cd = cd
        self.d = d
        self.l = l
        self.alpha = np.asarray([0, 0, 0.5 * rho * cd * d * l,0,0,0])  # x,y,z,rx,ry,rz

    # self.x = np.squeeze((self.a[nax, nax, :, :] @ rao[:, :, :, nax]), axis=3) * sea_spectrum[:, nax, nax]
    # self.eq = 8 / 3 * self.alpha[nax, nax, :] * w[:, nax, nax] * self.x[:, :, :] / np.pi


class ResponseModel(object):
    def __init__(self, candidate, Short_Term_Wave_Conditions):

        self._stwc = Short_Term_Wave_Conditions

        self._settings = candidate.settings
        self._loads = candidate.loads

        self.fe = self._loads._fe
        self.m = self._loads._m
        self.ma = self._loads._ma
        self.c_rad = self._loads._c_radiation * self._settings.radiaton_damping_factor

        self.k = self._loads._k
        self.w = self._loads._w
        self.nbeta = self._loads._nbeta
        self.nw = self._loads._nw

        self._rao_tra_mat = np.zeros([self._loads._nw, self._loads._nbeta, 4, 4], dtype=complex)

        self._c_visc = np.zeros([self.nw, 6, 6], dtype=complex)
        self._viscous_damper_centers = None

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

        self.get_wave_forces()
        self.calc_rao_dependent_response()

    def get_wave_forces(self):
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



    def calc_rao_dependent_response(self):
        self._rao_init = np.zeros([self.nw, self.nbeta, 6], dtype=complex)
        for ib in range(self.nbeta):  # self.nbeta
            self._rao_init[:, ib, :] = self.calc_rao_linear(ib)

        self._rao = self._rao_init.copy()
        if self._settings.do_linearize:
            self.create_viscous_damping()

        # RAO transformation matrix must be updated after all coefficients are calculated, including linearized viscous damping
        self.update_rao_tra_mat()

        # --------------------------------------------------------------------------------------------------------------
        # RADIATION
        # --------------------------------------------------------------------------------------------------------------
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
                                                                                   nax] * 1j

        self._panel_added_mass_force = np.sum(self._panel_added_mass_force_all_dof, axis=2)
        self._panel_radiation_damping_force = np.sum(self._panel_radiation_damping_force_all_dof, axis=2)
        # --------------------------------------------------------------------------------------------------------------
        # STIFFNESS
        # --------------------------------------------------------------------------------------------------------------

        self._panel_pos = np.zeros([self._loads._nw, self._loads._nbeta, self._loads.pd.npanels, 3], dtype=complex)
        # Gen dynamic position of panels and calc hydro static pressure
        # Append a 1 to the 3 dof vector to correspond with 4x4 tra_mat
        self._ppc = np.append(self._panel_pressure_centers, np.ones((self._loads.pd.npanels, 1)), 1)
        # Modify for broadcasting, add artificial dim to use matmul on stack of matrices
        self._ppc = self._ppc[nax, nax, :, :, nax]
        # Perform matmul and remove artificial dim and append 1. This code is fast ...
        self._panel_pos[:, :, :, :] = np.squeeze(np.matmul(self._rao_tra_mat[:, :, nax, :, :], self._ppc), axis=4)[:, :,
                                      :, 0:3]
        self._panel_pos -= self._loads.pd.ppanel_centers[nax, nax, :, :]  # Subtract mean position
        self._panel_diff_buoyancy_pressure = (self._loads.rho_sw * abs(self._loads.gravity)) * self._panel_pos[:, :, :,
                                                                                               2]  # Change in pressure
        self._panel_diff_buoyancy_force = (-self._panel_diff_buoyancy_pressure[:, :, :,
                                            nax] * self._projected_panel_area)

        # --------------------------------------------------------------------------------------------------------------
        # Mass
        # --------------------------------------------------------------------------------------------------------------
        self._point_mass_dynamic_inertia_force = np.squeeze(
                self._point_mass[nax, nax, :, :, :] @ self._rao[:, :, nax, :, nax],
                axis=4)

    def update_rao_tra_mat(self):
        # --------------------------------------------------------------------------------------------------------------
        # RAO TRANSFORMATION MATRIX
        # --------------------------------------------------------------------------------------------------------------
        for ifreq in range(self._loads._nw):  # TODO: Vectorize
            for ibeta in range(self._loads._nbeta):
                rot = self._rao[ifreq, ibeta, 3:6]
                tra = self._rao[ifreq, ibeta, 0:3]
                self._rao_tra_mat[ifreq, ibeta, :, :] = tb.transformation_matrix(rot, tra)

    def get_section_index(self, section_point, section_normal):

        def mask1(coordinates):
            vec = coordinates - section_point
            dot = np.dot(vec,
                         section_normal)  # dot product i positive for coordinates on the positive side of the plane
            return dot >= 0

        if self._settings.do_linearize:
            return mask1(self._point_mass_centers), mask1(self._panel_pressure_centers), mask1(self._viscous_damper_centers)
        else:
            return mask1(self._point_mass_centers), mask1(self._panel_pressure_centers), None

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

    def assemble_forces(self, imass, ipanel, idamp, moment_ref_point):

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

        if self._settings.do_linearize:
            f_viscous_damping = self.sum_forces(idamp, self._strip_viscous_damping_force, self._viscous_damper_centers,
                                                moment_ref_point)
            f_viscous_damping *= self.w[:, nax, nax, ]*1j
        else:
            f_viscous_damping=f_radiation_damping*0

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
        force_out['Dynamic']['Viscous damping'] = f_viscous_damping
        force_out['Dynamic']['Buoyancy'] = f_dz_s
        force_out['Dynamic']['SUM'] = f_fk + f_diff - (
                f_added_mass + f_radiation_damping + f_viscous_damping + f_dz_s + f_inertia)

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

    def create_viscous_damping(self):

        # for ib in range(1):  # self.nbeta


        self._n_strips = 10
        self.vd_list = []
        da = self._gaf * self._d_rc
        dx = da / self._n_strips
        dtheta = 2 * np.pi / self._nr
        theta = [i * dtheta for i in range(self._nr)]
        coord = []
        for ir in range(self._nr):
            rot_mat = tb.rotation_matrix([0, 0, theta[ir]])
            for ibay in range(self._nc):
                for istrip in range(self._n_strips):
                    coord.append(rot_mat @ [self._d_cc / 2 + self._d_rc * ibay + dx * (istrip + 1 / 2), 0, 0])

        c = np.asarray(coord)
        rho = self._settings._rho_sw
        cd = 2
        dia = self._d_rc + 0.2
        l = dx
        for coord in c:
            self.vd_list.append(Viscous_Damper(coord, rho, cd, dia, l))

        alphas = np.asarray([vd.alpha for vd in self.vd_list])

        # Append a 1 to the 3 dof vector to correspond with 4x4 tra_mat
        self._viscous_damper_centers = c
        c_4 = np.append(c, np.ones((c.shape[0], 1)), 1)
        amp = self._stwc.hs / 2 # Damping target
        wp = self._stwc.wp

        ndof = 6

        self._c_visc = np.zeros([self.nbeta, ndof, ndof], dtype=complex)

        for i in range(12):

            self.update_rao_tra_mat()


            rao_wp = np.zeros([self.nbeta, ndof], dtype=complex)
            rao_tra_mat_wp = np.zeros([self.nbeta, 4, 4], dtype=complex)
            for ibeta in range(self.nbeta):
                for idof in range(ndof):
                    rao_wp[ibeta, idof] = np.interp(wp, self.w, self._rao[:, ibeta, idof])

                rot = rao_wp[ibeta, 3:6]
                tra = rao_wp[ibeta, 0:3]
                rao_tra_mat_wp[ibeta, :, :] = tb.transformation_matrix(rot, tra)

            # Gen dynamic position of dampers and calc damping coefficient

            pos = np.squeeze(np.matmul(rao_tra_mat_wp[:, nax, :, :], c_4[nax, :, :, nax]), axis=3)
            x = np.zeros([self.nbeta, alphas.shape[0],ndof], dtype=complex)
            x[:, :, 0:3] = (pos[:, :, 0:3] - c[nax, :, :]) * amp  # Subtract mean position and multiply with wave amp


            f_d_local = 8 / (3 * np.pi) * (alphas * (wp * x) ** 2) * complex(0, 1)  # Linearized damping force

            m_d_global = np.cross(c[nax, :, :], f_d_local[:, :, 0:3])
            f_d_global = np.concatenate((f_d_local[:,:,:3], m_d_global), axis=2)

            f_d_global_sum = np.sum(f_d_global, axis=1)

            vel = (rao_wp * amp * wp * complex(0, 1))

            c_visc_1d = f_d_global_sum / vel


            for ib in range(self.nbeta):  # self.nbeta
                self._c_visc[ib,:,:] = np.diag(c_visc_1d[ib,:])
                self._rao[:, ib, :] = self.calc_rao_linear(ib)

            c_visc_local = (f_d_local/ vel[:,nax,:])
            self._strip_viscous_damping_force = (c_visc_local[nax,:,:,:]*self._rao[:, :,nax,:])[:,:,:,:3]


    def calc_rao_linear(self, ibeta):
        out_rao = np.zeros([len(self.w), 6], dtype=complex)
        this_c_visc = self._c_visc[ibeta, :, :]
        for i, w in enumerate(self.w):
            #print('w = {}'.format(w))
            this_ma = self.ma[i, :, :]

            # print(this_c_visc)
            this_c = self.c_rad[i, :, :] + this_c_visc
            denom = np.asarray(-w ** 2 * (self.m + this_ma) + 1j * w * this_c + self.k, dtype=complex)
            daf = scipy.linalg.inv(denom)
            x = daf @ self.fe[i, ibeta, :]
            out_rao[i, :] = x
        return out_rao

    def point_rao(self,c):
        _c = np.append(c, np.ones((c.shape[0], 1)), 1)
        # Modify for broadcasting, add artificial dim to use matmul on stack of matrices
        _c = _c[nax, nax, :, :, nax]
        # Perform matmul and remove artificial dim and append 1. This code is fast ...
        return np.squeeze(np.matmul(self._rao_tra_mat[:, :, nax, :, :], _c), axis=4)[:, :,
                                      :, 0:3]



        





    @property
    def rao(self):
        return self._rao

    @property
    def rao_init(self):
        return self._rao_init