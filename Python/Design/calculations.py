__author__ = "Eivind Sonju"
__copyright__ = "Copyright (C) 2017-2019 Verbun AS. All rights reserved."
__version__ = "2.0"

import numpy as np
import scipy
import nemoh
import pickle
import tool_box as tb


class TransferFunctions(object):
    def __init__(self, settings, loads):
        # super().__init__(settings)
        self._settings = settings
        self._loads = loads

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

        # Set up RAOs
        self._rao = np.zeros([self._loads.nw, self._loads._nbeta, 6], dtype=complex)
        for ibeta in range(self._loads._nbeta):
            self._rao[:, ibeta, :] = self.get_rao(ibeta)

        # Prepare discrete mass and hydro forces
        with open(self._settings.fio.data_io_dir.joinpath('unit_model.pkl'), 'rb') as f:
            self._unit_model = pickle.load(f)

        self._part_list = self._unit_model.get_parts_without_children()
        self._point_mass = np.asarray([part._inertias.mass_matrix_global for part in self._part_list])

        self._point_mass_centers = np.asarray([-part._inertias.reduction_point for part in self._part_list])
        self._panel_pressure_centers = self._loads._pd.ppanel_centers

        self._projected_panel_area = self._loads._an[np.newaxis, np.newaxis, :, :]

        # Gravity
        self._point_mass_gravity_force = self._point_mass[:, 2, 2][:, np.newaxis] * np.asarray(
                [0, 0, self._settings.gravity])[
                                                                                    np.newaxis, :]  # Use m33

        # Hydro static / Buoyancy
        self._panel_pressure_buoyancy_force = self._loads._force['Buoyancy']

        # Froude-Krylof
        self._panel_pressure_froude_krylof_force = self._loads._force['Froude-Krylof']

        # Diffraction
        self._panel_pressure_diffraction_force = self._loads._force['Diffraction']

        # Radiation
        self._panel_pressure_radiation_unit_force = self._loads._force['Radiation'].copy()  # Dont mess with original
        # Here the RAO for each DOF is multiplied with each RAO dependent panel force (x,y,z)
        self._panel_pressure_radiation_force_all_dof = self._panel_pressure_radiation_unit_force[:, np.newaxis, :, :,
                                                       :] * self._rao[:, :, :, np.newaxis, np.newaxis]
        self._panel_pressure_radiation_force = np.sum(self._panel_pressure_radiation_force_all_dof, axis=2)

        # RAO transformation matrix
        self._rao_tra_mat = np.zeros([self._loads._nw, self._loads._nbeta, 4, 4], dtype=complex)
        self._panel_pos = np.zeros([self._loads._nw, self._loads._nbeta, self._loads.pd.npanels, 3], dtype=complex)
        for ifreq in range(self._loads._nw):
            for ibeta in range(self._loads._nbeta):
                rot = self._rao[ifreq, ibeta, 3:6]
                tra = self._rao[ifreq, ibeta, 0:3]
                self._rao_tra_mat[ifreq, ibeta, :, :] = tb.transformation_matrix(rot, tra)  # TODO: Vectorize
        self._rtm = self._rao_tra_mat[:, :, np.newaxis, :, :]

        # Gen dynamic position of panels and calc hydro static pressure
        # Append a 1 to the 3 dof vector to correspond with 4x4 tra_mat
        self._pc = np.append(self._panel_pressure_centers, np.ones((self._loads.pd.npanels, 1)), 1)
        # Modify for broadcasting, add artificial dim to use matmul on stack of matrices
        self._pc = self._pc[np.newaxis, np.newaxis, :, :, np.newaxis]
        # Perform matmul and remove artificial dim and appended 1. This code is fast ...
        self._panel_pos[:, :, :, :] = np.squeeze(np.matmul(self._rtm, self._pc), axis=4)[:, :, :, 0:3]
        self._panel_pos -= self._loads.pd.ppanel_centers[np.newaxis, np.newaxis, :]  # Subtract mean position
        # panel_pos += self._rao[:, :, np.newaxis, 0:3]
        self._pressures_diff_static = (self._loads.rho_sw * abs(self._loads.gravity)) * self._panel_pos[:, :, :,
                                                                                        2]  # Change in pressure
        self._panel_pressure_diff_static_force = -self._pressures_diff_static[:, :, :,
                                                  np.newaxis] * self._projected_panel_area

        # Get dynamic acceleration of part masses and calc inertia force
        # TODO: Include rotation/inertia moment from parts
        self._pmc = np.append(self._point_mass_centers, np.ones((self._point_mass_centers.shape[0], 1)), 1)
        self._pmc = self._pmc[np.newaxis, np.newaxis, :, :, np.newaxis]
        self._dynamic_point_mass_pos = np.squeeze(np.matmul(self._rtm, self._pmc), axis=4)[:, :, :, 0:3]
        self._dynamic_point_mass_pos -= self._point_mass_centers[np.newaxis, np.newaxis, :]  # Subtract mean position
        self._w2 = self._loads.w ** 2
        self._part_dynacc = self._dynamic_point_mass_pos * self._w2[:, np.newaxis, np.newaxis, np.newaxis]
        self._point_mass_dynamic_inertia_force = -self._part_dynacc * self._point_mass[np.newaxis, np.newaxis, :,
                                                                      np.newaxis, 0, 0]

    def get_rao(self, idir):
        fe = self._loads._fe[:, idir, :]
        m = self._loads._m
        ma = self._loads._ma
        c = self._loads._c_hyd
        k = self._loads._k
        return self.rao(fe, m, ma, c, k, self._loads._w)

    def rao(self, f, m, ma, c, k, vw):
        container = np.zeros([len(vw), 6], dtype=complex)
        for i, w in enumerate(vw):
            this_ma = ma[i, :, :]
            denom = np.asarray(-w ** 2 * (m + this_ma) + 1j * w * c[i, :, :] + k, dtype=complex)
            daf = scipy.linalg.inv(denom)
            # daf =np.asarray([[1/x for x in col]for col in denom])
            x = daf @ f[i, :]
            container[i, :] = x

        return container

    def get_section_index(self, section_point, section_normal):

        def mask(coordinates):
            vec = coordinates - section_point
            dot = np.dot(vec,
                         section_normal)  # dot product i positive for coordinates on the positive side of the plane
            return dot >= 0

        return mask(self._point_mass_centers), mask(self._panel_pressure_centers)

    def get_flange_panel_index(self):
        def printv(string):
            if 0:
                print(string)

        class BreakIt(Exception):
            pass

        def mask(coordinates):
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

        return mask(self._point_mass_centers), mask(self._panel_pressure_centers)

    def sum_forces(self, index, forces, coordinates, moment_ref_point):

        vec = coordinates - moment_ref_point

        if 1 in index:  # At least one item is on the considered side of the section surface
            if forces.ndim == 4:  # Dynamic [freq, dir, panel, f]
                section_forces = forces[:, :, index, :]
                section_moments = np.cross(vec[np.newaxis, np.newaxis, index, :],
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

        f_gravity = self.sum_forces(imass, self._point_mass_gravity_force, self._point_mass_centers, moment_ref_point)
        f_inertia = self.sum_forces(imass, self._point_mass_dynamic_inertia_force, self._point_mass_centers,
                                    moment_ref_point)
        f_bouyancy = self.sum_forces(ipanel, self._panel_pressure_buoyancy_force, self._panel_pressure_centers,
                                     moment_ref_point)
        f_fk = self.sum_forces(ipanel, self._panel_pressure_froude_krylof_force, self._panel_pressure_centers,
                               moment_ref_point)
        f_diff = self.sum_forces(ipanel, self._panel_pressure_diffraction_force, self._panel_pressure_centers,
                                 moment_ref_point)
        f_rad = self.sum_forces(ipanel, self._panel_pressure_radiation_force, self._panel_pressure_centers,
                                moment_ref_point)
        f_dz_s = self.sum_forces(ipanel, self._panel_pressure_diff_static_force, self._panel_pressure_centers,
                                 moment_ref_point)

        force_out = dict()
        force_out['Static'] = dict()
        force_out['Dynamic'] = dict()

        force_out['Static']['Total'] = f_gravity + f_bouyancy
        force_out['Static']['Gravity'] = f_gravity
        force_out['Static']['Buoyancy'] = f_bouyancy
        force_out['Dynamic']['Total'] = f_rad + f_fk + f_diff + f_dz_s + f_inertia
        force_out['Dynamic']['Radiation'] = f_rad
        force_out['Dynamic']['Froude-Krylof'] = f_fk
        force_out['Dynamic']['Diffraction'] = f_diff
        force_out['Dynamic']['Buoyancy'] = f_dz_s
        force_out['Dynamic']['Inertia'] = f_inertia

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

    def flange_normal_stress(self, f):

        # Get forces at section center
        # TODO: Change z to section center, now at waterline

        a_reinf = 2 * 4 * self._t_uf * 3 * self._t_uf
        a = self._d_rc * self._t_uf + a_reinf
        wz = self._t_uf * self._d_rc ** 2 / 6 + a_reinf * self._d_rc / 2

        def sig(f):
            if f.ndim == 3:
                sig_ax = f[:, :, 0] / (2 * a)
                sig_by = f[:, :, 4] / (self._h * a)
                sig_bz = f[:, :, 5] / (2 * wz)
            else:
                sig_ax = f[0] / (2 * a)
                sig_by = f[4] / (self._h * a)
                sig_bz = f[5] / (2 * wz)
            return sig_ax + sig_by + sig_bz

        return sig(f)



    def flange_lateral_pressure(self, f):
        return f/(self._d_rc *self._d_rc * self._gaf)



class Short_Term_Wave_Conditions():
    def __init__(self, hs, tp):
        self._hs = hs
        self._tp = tp
        self._gamma = self.gamma()
        self._tz = self.tp2tz()

    def gamma(self):
        x = self._tp / np.sqrt(self._hs)
        if x <= 3.6:
            return 5
        elif x < 5:
            return np.exp(5.75 - 1.15 * x)
        else:
            return 1

    def s_jonswap(self, w, gamma=None):
        sig_a = 0.07
        sig_b = 0.09
        delta_sig = sig_b - sig_a
        wp = 2 * np.pi / self._tp

        if gamma == None:
            gamma = self._gamma

            # print('Gamma {}'.format(gamma))

        a_gamma = 1 - 0.287 * np.log(gamma)

        def spec_pm(w):
            return (5 / 16) * (self._hs ** 2) * (wp ** 4) * (w ** (-5)) * np.exp(-(5 / 4) * ((w / wp) ** (-4)))

        def spec_j(w):
            def sig(w):
                return sig_a if w <= wp else sig_b

            sig_ab = np.array(list(map(sig, w)))

            return a_gamma * spec_pm(w) * gamma ** np.exp(-0.5 * ((w - wp) / sig_ab * wp))

        if gamma == 1:
            return spec_pm(w)
        else:
            return spec_j(w)

    def tp2tz(self):
        return (0.6673 + 0.05037 * self._gamma - 0.006230 * self._gamma ** 2 + 0.0003341 * self._gamma ** 3) * self._tp

    def expected_largest_maximum(self, x, w, gamma=None):
        x_r = np.abs(x ** 2) * self.s_jonswap(w, gamma)[:, np.newaxis]
        dw = w[1] - w[0]
        sig_r_m0 = sum(x_r) * dw
        nz = 3 * 3600 / self._tz
        return np.sqrt(sig_r_m0) * (np.sqrt(2 * np.log(nz)) + 0.5772 / np.sqrt(2 * np.log(nz)))


class DNVGL_RP_C201_Part1():
    def __init__(self, settings):
        self._settings = settings
        self._info = dict()
        self._info['Company'] = 'Det Norske Veritas'
        self._info['Type'] = 'Recommended Practice'
        self._info['Title'] = 'Buckling Strength Of Plated Structures'
        self._info['Version'] = 'October 2010'

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
        self._n_lfst = self._settings.job_data['floater']['Radial']['Flange']['Lower']['Stiffener']['Longitudinal'][
            'Number of']

        self._gaf = self._settings.job_data['floater']['Gap factor']

    def eq_7_1(self, sig_x_sd, tau_sd):
        # Return the equivalent axial force: N_Sd
        s = self._d_rc / (self._n_lfst - 1)  # spacing between stiffeners
        l = self._d_rc * self._gaf  # length of plate
        e = self._settings.emod_st  # modulus of elasticity for steel
        t = self._t_lf  # plate thickness
        l_g = self._d_rc  # Girder length, ref. fig 3.1 in C201
        gamma_m = self._settings.gamma_m  # material factor

        # eq. 7.7
        if l >= s:
            k_l = 5.34 + 4 * (s / l) ** 2
        else:
            k_l = 5.34 * (s / l) ** 2 + 4

        # eq. 7.6
        tau_crl = k_l * 0.904 * e * (t / s) ** 2

        # eq. 7.5
        if l <= l_g:
            k_g = 5.34 + 4 * (l / l_g) ** 2
        else:
            k_g = 5.34 * (l / l_g) ** 2 + 4

        # eq. 7.4
        tau_crg = k_g * 0.904 * e * (t / l) ** 2

        # eq. 7.2 and 7.3
        if tau_sd >= (tau_crl / gamma_m):
            tau_tf = tau_sd - tau_crg
        else:
            tau_tf = 0

        # eq. 7.1
        a_s = self._t_lfst * self._h_lfst
        n_sd = sig_x_sd * (a_s + s * t) + tau_tf * s * t

        return n_sd

    def eq_7_8(self, p_Sd, f_y=None, gamma_m=None):
        s = self._d_rc / (self._n_lfst - 1)  # spacing between stiffeners
        # Return the equivalent lateral line load: q_Sd
        return p_Sd*s

    def dummy(self, p_Sd, f_y=None, gamma_m=None):
        # Return the equivalent lateral line load: q_Sd

        s = self._d_rc / (self._n_lfst - 1)  # spacing between stiffeners
        l = self._d_rc * self._gaf  # length of plate
        e = self._settings.emod_st  # modulus of elasticity for steel
        t = self._t_lf  # plate thickness
        l_g = self._d_rc  # Girder length, ref. fig 3.1 in C201
        w_p = s  # Width of plate
        t_s = self._t_lfst  # Width of stiffener
        h_s = self._h_lfst  # Height of stiffener

        def i(b, h):
            return b * h ** 3 / 12

        def i_f(w_p, t, t_s, h_s):
            a1 = w_p * t
            a2 = t_s * h_s
            y_ = ((a1 * t / 2) + (a2 * (h_s / 2 + t))) / (a1 + a2)
            i1 = i(w_p, t)
            i2 = i(t_s, h_s)
            i_s = i1 + a1 * (y_ - t / 2) + i2 + a2 * (y_ - (t + h_s / 2))
            return y_, i_s



        # eq. 7.12
        y_, i_s = i_f(w_p, t, t_s, h_s)
        k_c = 2 * (1 + np.sqrt(1 + (10.9 * i_s) / (t ** 3 * s)))

        # eq. 7.11
        m_c = 8.9  # 13.3 for continuous stiffeners or 8.9 for simple supported stiffeners (sniped stiffeners)
        s_e = s * self.eq_7_13(f_y)  # Effective plate
        y_, i_es = i_f(s_e, t, t_s, h_s)
        w_es = i_es/y_

    def eq_7_13(self, f_y):
        # The effective plate width for a continuous stiffener subjected
        # to longitudinal and transverse stress and shear

        s = self._d_rc / (self._n_lfst - 1)  # spacing between stiffeners
        e = self._settings.emod_st  # modulus of elasticity for steel
        t = self._t_lf  # plate thickness

        # eq. 7.16
        c_ys = 1  # No reduction factor due to transverse stress

        # eq. 7.15
        lambda_p = 0.525 * (s / t) * np.sqrt(f_y * e)

        # eq. 7.14
        if lambda_p > 0.673:
            c_xs = (lambda_p - 0.22) / lambda_p ** 2
        else:
            c_xs = 1.0

        return c_ys * c_xs
