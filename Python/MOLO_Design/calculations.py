__author__ = "Eivind Sonju"
__copyright__ = "Copyright (C) 2017-2019 Verbun AS. All rights reserved."
__version__ = "2.0"

import numpy as np
import scipy
import nemoh
import pickle
import tool_box as tb




class TransferFunctions(object):
    def __init__(self, settings,loads):
        #super().__init__(settings)
        self._settings = settings
        self._loads = loads
        self_env = Environment(settings)
        self._rao = np.zeros([self._loads.nw, self._loads._nbeta, 6], dtype=complex)
        for ibeta in range(self._loads._nbeta):
            self._rao[:, ibeta, :] = self.get_rao(ibeta)

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

        # def H(self, f, m, ma, c, k, vw):
        #     container = np.zeros([len(vw), 6], dtype=complex)
        #     for iw, w in enumerate(vw):
        #         for i in [2,4]:
        #             this_ma=ma[iw, i, i]
        #             denom = np.asarray(-w ** 2 * (m[i,i] + ma[iw, i, i]) + 1j * w * c[iw, i, i] + k[i,i], dtype=complex)
        #             daf = 1/denom
        #             # daf =np.asarray([[1/x for x in col]for col in denom])
        #             x = daf * f[iw, i]
        #             container[iw, i] = x

        return container

    def section_forces(self, section_point=None, section_normal=None):

        if section_point is None:
            section_point = np.array([1, 0, 0])
        if section_normal is None:
            section_normal = np.array([1, 0, 0])

        section_point = np.asarray(section_point)
        section_normal = np.asarray(section_normal)

        print('\n--------------------------------------------------------------------------------------------')
        print('SECTION FORCES')
        print('Point: {}\nNormal: {}'.format(np.array2string(section_point), np.array2string(section_normal)))
        print('--------------------------------------------------------------------------------------------')

        with open(self._settings.fio.data_io_dir.joinpath('unit_model.pkl'), 'rb') as f:
            unit_model = pickle.load(f)
        unit_model.set_new_reduction_point(section_point)

        # Prepare RAO's for motion dependent response variables

        part_list = unit_model.get_parts_without_children()
        # [print(part._type) for part in part_list]

        point_mass = np.asarray([part._inertias.mass_matrix_global for part in part_list])

        # Gravity
        point_mass_gravity_force = point_mass[:, 2, 2][:, np.newaxis] * np.asarray([0, 0, self._settings.gravity])[
                                                                        np.newaxis, :]  # Use m33
        point_mass_centers = np.asarray([-part._inertias.reduction_point for part in part_list])
        f_gravity = nemoh.get_section_values(point_mass_gravity_force, point_mass_centers, section_point,
                                             section_normal)
        print('\nGravity force')
        tb.matprint(f_gravity)

        # Hydro static / Buoyancy
        f_bouyancy = nemoh.get_section_values(self._loads._force['Buoyancy'], self._loads._pd.ppanel_centers,
                                              section_point,
                                              section_normal)
        print('\nBuoyancy force')
        tb.matprint(f_bouyancy)

        # Collect all forces acting on the section
        # Froude-Krylof and diffraction
        f_fk = nemoh.get_section_values(self._loads._force['Froude-Krylof'], self._loads.pd.ppanel_centers, section_point,
                                        section_normal)
        f_diff = nemoh.get_section_values(self._loads._force['Diffraction'], self._loads.pd.ppanel_centers, section_point,
                                          section_normal)

        # Calculate radiation force transferfunctions R = H * eta

        # get_rao create complex motion at origin per freq in all dofs for given wave dir

        rad = self._loads._force['Radiation'].copy()  # Dont mess with original
        # rad = np.swapaxes(rad,1,2) # Modify radiation axes to align axes with RAO for broadcasting
        f_rad_all_panels_all_dof = rad[:, np.newaxis, :, :, :] * self._rao[:, :, :, np.newaxis, np.newaxis]
        f_rad_all_panels = np.sum(f_rad_all_panels_all_dof, axis=2)  # TODO: Move this work to init. Will be reused

        # Here the RAO for each DOF is multiplied with each RAO dependent panel force (x,y,z)
        f_rad = nemoh.get_section_values(f_rad_all_panels, self._loads.pd.ppanel_centers,
                                         section_point,
                                         section_normal)

        # Gen dynamic position of panels and calc hydro static pressure
        rao_tra_mat = np.zeros([self._loads._nw, self._loads._nbeta, 4, 4], dtype=complex)
        panel_pos = np.zeros([self._loads._nw, self._loads._nbeta, self._loads.pd.npanels, 3], dtype=complex)
        for ifreq in range(self._loads._nw):
            for ibeta in range(self._loads._nbeta):
                rot = self._rao[ifreq, ibeta, 3:6]
                tra = self._rao[ifreq, ibeta, 0:3]
                rao_tra_mat[ifreq, ibeta, :, :] = tb.transformation_matrix(rot, tra)  # TODO: Vectorize

        rtm = rao_tra_mat[:, :, np.newaxis, :, :]
        # Append a 1 to the 3 dof vector to correspond with 4x4 tra_mat
        pc = np.append(self._loads.pd.ppanel_centers, np.ones((self._loads.pd.npanels, 1)), 1)
        # Modify for broadcasting, add artificial dim to use matmul on stack of matrices
        pc = pc[np.newaxis, np.newaxis, :, :, np.newaxis]
        # Perform matmul and remove artificial dim and appended 1. This code is fast ...
        panel_pos[:, :, :, :] = np.squeeze(np.matmul(rtm, pc), axis=4)[:, :, :, 0:3]
        panel_pos -= self._loads.pd.ppanel_centers[np.newaxis, np.newaxis, :]  # Subtract mean position
        # panel_pos += self._rao[:, :, np.newaxis, 0:3]
        dp = (self._loads._rho_sw * abs(self._loads._grav)) * panel_pos[:, :, :, 2]  # Change in pressure
        f_dz = -dp[:, :, :, np.newaxis] * self._loads._an[np.newaxis, np.newaxis, :, :]

        f_dz_s = nemoh.get_section_values(f_dz, self._loads.pd.ppanel_centers, section_point,
                                          section_normal)

        # Get dynamic acceleration of part masses and calc inertia force
        # TODO: Include rotation/inertia moment from parts

        pmc = np.append(point_mass_centers, np.ones((point_mass_centers.shape[0], 1)), 1)
        pmc = pmc[np.newaxis, np.newaxis, :, :, np.newaxis]
        point_mass_pos = np.squeeze(np.matmul(rtm, pmc), axis=4)[:, :, :, 0:3]

        point_mass_pos -= point_mass_centers[np.newaxis, np.newaxis, :]
        w2 = self._loads.w ** 2
        part_dynacc = point_mass_pos * w2[:, np.newaxis, np.newaxis, np.newaxis]
        part_inertia_force = -part_dynacc * point_mass[np.newaxis, np.newaxis, :, np.newaxis, 0, 0]
        f_inertia = nemoh.get_section_values(part_inertia_force,
                                             point_mass_centers, section_point,
                                             section_normal)

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

        # Neutral axis relative to bottom of cylinder
        #z0 =




    def section_stress(self,f):


        # Get forces at section center
        # TODO: Change z to section center, now at waterline

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

        a_reinf = 2*4*t_uf*3*t_uf
        a = w * t_uf + a_reinf
        wz = t_uf * w ** 2 / 6  +  a_reinf * w/2

        def sig(f):
            if f.ndim == 3:
                sig_ax = f[:, :, 0] / (2*a)
                sig_by = f[:, :, 4] / (h*a)
                sig_bz = f[:, :, 5] / (2 * wz)
            else:
                sig_ax = f[0] / (2*a)
                sig_by = f[4] / (h*a)
                sig_bz = f[5] / (2 * wz)
            return sig_ax + sig_by + sig_bz

        return sig(f)


class Environment():
    def __init__(self, settings):
        pass

    def gamma(self,hs,tp):
        x = tp / np.sqrt(hs)
        if x <= 3.6:
            return  5
        elif x < 5:
            return np.exp(5.75 - 1.15 * x)
        else:
            return  1

    def s_jonswap(self, hs, wp, w, gamma=None):
        sig_a = 0.07
        sig_b = 0.09
        delta_sig = sig_b - sig_a

        tp = 2 * np.pi / wp
        if gamma == None:
            gamma = self.gamma(hs,tp)

            # print('Gamma {}'.format(gamma))

        a_gamma = 1 - 0.287 * np.log(gamma)

        def spec_pm(w):
            return (5 / 16) * (hs ** 2) * (wp ** 4) * (w ** (-5)) * np.exp(-(5 / 4) * ((w / wp) ** (-4)))

        def spec_j(w):
            def sig(w):
                return sig_a if w <= wp else sig_b

            sig_ab = np.array(list(map(sig, w)))

            return a_gamma * spec_pm(w) * gamma ** np.exp(-0.5 * ((w - wp) / sig_ab * wp))

        if gamma == 1:
            return spec_pm(w)
        else:
            return spec_j(w)

    def tp2tz(self,tp,gamma):
        return (0.6673 + 0.05037*gamma - 0.006230*gamma**2 + 0.0003341*gamma**3)*tp






class DNVGL_RP_C201():
    def __init__(self):
        self._info = dict()
        self._info['Company'] = 'Det Norske Veritas'
        self._info['Type'] = 'Recommended Practice'
        self._info['Title'] = 'Buckling Strength Of Plated Structures'
        self._info['Version'] = 'October 2010'

    def func_sx_rd(self, cx, fy, gm):  # Equation 6.1
        return cx * fy / gm

    def func_cx(self, l_p):  # Equation 6.2
        if l_p <= 0.673:
            return 1
        else:
            return (l_p - 0.22) / l_p ** 2

    def func_l_p(self, s, t, fy, e):  # Equation 6.3
        return 0.525 * (s / t) * np.sqrt(fy / e)

    # def func_
