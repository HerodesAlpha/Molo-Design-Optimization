__author__ = "Eivind Sonju"
__copyright__ = "Copyright (C) 2017-2019 Verbun AS. All rights reserved."
__version__ = "2.0"

import numpy as np
import scipy
import nemoh
import pickle
import tool_box as tb

from loads import Sea_and_Inertia_Loads


class TransferFunctions(Sea_and_Inertia_Loads, object):
    def __init__(self, settings):
        super().__init__(settings)

        # def rao(self, fe, m, ma, c, k, w):
        #     return np.absolute(fe / (-w ** 2 * (m + ma) + 1j * w * (c) + k))
        #
        #
        # def get_rao(self, idof, idir):
        #     fe = self._fe[:, idir, idof]
        #     m = self._m[idof][idof]
        #     ma = self._ma[:, idof, idof]
        #     c = self._c_hyd[:, idof, idof]
        #     k = self._k[idof][idof]
        #     return self.rao(fe, m, ma, c, k, self._w)

        self._rao = np.zeros([self.nw, self._nbeta, 6], dtype=complex)
        for ibeta in range(self._nbeta):
            self._rao[:, ibeta, :] = self.get_rao(ibeta)

    def get_rao(self, idir):
        fe = self._fe[:, idir, :]
        m = self._m
        ma = self._ma
        c = self._c_hyd
        k = self._k
        return self.rao(fe, m, ma, c, k, self._w)

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

        if section_point == None:
            section_point = np.array([1, 0, 0])
        if section_normal == None:
            section_normal = np.array([1, 0, 0])

        print('\n--------------------------------------------------------------------------------------------')
        print('SECTION FORCES')
        print('Point: {}\nNormal: {}'.format(np.array2string(section_point), np.array2string(section_normal)))
        print('--------------------------------------------------------------------------------------------')

        with open(self._settings.fio.data_io_dir.joinpath('unit_model.pkl'), 'rb') as f:
            unit_model = pickle.load(f)
        unit_model.set_new_reduction_point(section_point)
        # Prepare RAO's for motion dependent response variables

        part_list = unit_model.get_parts()
        # [print(part._type) for part in part_list]

        point_mass = np.asarray([part._inertias.mass_matrix_global for part in part_list])

        # Gravity
        point_mass_gravity_force = point_mass[:, 2, 2][:, np.newaxis] * np.asarray([0, 0, self._settings.grav])[
                                                                        np.newaxis, :]  # Use m33
        point_mass_centers = np.asarray([-part._inertias.reduction_point for part in part_list])
        f_gravity = nemoh.get_section_values(point_mass_gravity_force, point_mass_centers, section_point,
                                             section_normal)
        print('\nGravity force')
        tb.matprint(f_gravity)

        # Hydro static / Buoyancy
        f_bouyancy = nemoh.get_section_values(self._force['Buoyancy'], self._pd.ppanel_centers,
                                              section_point,
                                              section_normal)
        print('\nBuoyancy force')
        tb.matprint(f_bouyancy)

        # Collect all forces acting on the section
        # Froude-Krylof and diffraction
        f_fk = nemoh.get_section_values(self._force['Froude-Krylof'], self.pd.ppanel_centers, section_point,
                                        section_normal)
        f_diff = nemoh.get_section_values(self._force['Diffraction'], self.pd.ppanel_centers, section_point,
                                          section_normal)

        # Calculate radiation force transferfunctions R = H * eta

        # get_rao create complex motion at origin per freq in all dofs for given wave dir


        rad = self._force['Radiation'].copy()  # Dont mess with original
        # rad = np.swapaxes(rad,1,2) # Modify radiation axes to align axes with RAO for broadcasting
        f_rad_all_panels_all_dof = rad[:, np.newaxis, :, :, :] * self._rao[:, :, :, np.newaxis, np.newaxis]
        f_rad_all_panels = np.sum(f_rad_all_panels_all_dof, axis=2)  # TODO: Move this work to init. Will be reused

        # Here the RAO for each DOF is multiplied with each RAO dependent panel force (x,y,z)
        f_rad = nemoh.get_section_values(f_rad_all_panels, self.pd.ppanel_centers,
                                         section_point,
                                         section_normal)

        # Gen dynamic position of panels and calc hydro static pressure

        rao_rot_mat = np.zeros([self._nw, self._nbeta, 3, 3], dtype=complex)
        panel_pos = np.zeros([self._nw, self._nbeta, self.pd.npanels, 3], dtype=complex)
        for ifreq in range(self._nw):
            for ibeta in range(self._nbeta):
                rao_rot_mat[ifreq, ibeta, :, :] = tb.rotation_matrix(self._rao[ifreq, ibeta, 3:6])  # TODO: Vectorize
                """
                for ipanel in range(self.pd.npanels):
                    panel_pos[ifreq, ibeta, ipanel, :] = np.transpose(
                        np.dot(rao_rot_mat[ifreq, ibeta, :, :], self.pd.ppanel_centers[ipanel, :].T))



        """
        a = rao_rot_mat[:, :, np.newaxis, :, :]
        b = self.pd.ppanel_centers[np.newaxis, np.newaxis, :, :, np.newaxis] # Modify for broadcasting, add artificial dim to use matmul on stack of matrices
        panel_pos[:, :, :, :] = np.squeeze(np.matmul(a, b)) # Perform matmul and remove artificial dim. This code is fast ...

        # panel_pos = np.transpose(np.dot(rao_rot_mat[:,:,:,:], self.pd.ppanel_centers.T))
        # panel_pos = np.transpose(np.matmul(rao_rot_mat[:,:,np.newaxis,:,:], self.pd.ppanel_centers[np.newaxis,np.newaxis,:,:, np.newaxis]))
        panel_pos += self._rao[:, :, np.newaxis, 0:3]
        p_dz = (self._rho_sw * self._grav) * panel_pos[:, :, :, 2]
        f_dz = -p_dz[:, :, :, np.newaxis] * self._an[np.newaxis, np.newaxis, :, :]

        f_varying_buoyancy = nemoh.get_section_values(f_dz, self.pd.ppanel_centers, section_point,
                                                      section_normal)

        """
            # Get dynamic acceleration of part masses and calc inertia force
            f_inertia = np.zeros([self.nw, 6], dtype=complex)

            part_dynpos = np.transpose(np.dot(rao_rot_mat, point_mass_centers.T))
            part_dynpos += rao[ifreq, 0:3]
            part_dynacc = part_dynpos * self.w[ifreq] ** 2
            part_inertia_force = part_dynacc * point_mass[:, np.newaxis]
            f_inertia[ifreq, :] = nemoh.get_section_values(part_inertia_force,
                                                           point_mass_centers, section_point,
                                                           section_normal)

        f_tot_dyn[:, idir, :] = f_fk + f_diff + f_rad + f_varying_buoyancy + f_inertia
        """

        return f_rad + f_fk + f_diff + f_varying_buoyancy
