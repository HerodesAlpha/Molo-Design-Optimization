import numpy as np
import environmental_conditions as ec
from scipy.optimize import minimize
from scipy.optimize import Bounds
from scipy.optimize import broyden1,broyden2, newton_krylov

class BreakIt(Exception): pass


# Code check of t

# Tripping of stiffener, use stiffened plate criteria?

class Panel():
    def __init__(self, settings):
        self._settings = settings
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

        self._l = self._d_rc * self._gaf

        self._w_p = None
        self._a = None
        self._y_e = None
        self._i = None
        self._z_y = None
        self._w_z = None

        self.init_cross_section()

    def init_cross_section(self):
        def i(b, h):
            return b * h ** 3 / 12

        e = self._settings.emod_st  # modulus of elasticity for steel
        w_p = self._d_rc  # Width of plate
        t_p = self._t_lf  # Thickness of plate
        n_s = self._n_lfst  # Number of stiffeners
        t_s = self._t_lfst  # Width of stiffener
        h_s = self._h_lfst  # Height of stiffener

        # Section areas
        a_p = w_p * t_s
        a_s = h_s * t_s
        a_tot = a_p + n_s * a_s

        # Elastic neutral axis
        y_e = ((a_p * t_p / 2) + n_s * (a_s * (h_s / 2 + t_p))) / (a_p + n_s * a_p)

        # Moment of inertia
        i_p = i(w_p, t_p)
        i_s = i(t_s, h_s)
        i_panel = i_p + a_p * (y_e - t_p / 2) ** 2 + i_s + a_s * (y_e - (t_p + h_s / 2)) ** 2

        # Plastic section modulus about y (in plane, lateral to radial)
        if a_p / a_tot >= 0.5:  # Is neutral axis in plate?
            y_p = (t_p + n_s * a_s / w_p) / 2
            z_y = a_tot * y_p / 2
        else:  # or in stiffener
            y_p = t_p + (a_tot / 2 - a_p) / (n_s * t_s)
            z_y = a_tot * (t_p + h_s - y_p) / 2

        # Elastic section modulus about z (vertical)
        w_z = self._t_uf * self._d_rc ** 2 / 6
        dw = w_p / (n_s - 1)
        for i in range(int(n_s / 2)):  # Number of stiffeners cannot be odd
            # Add in pairs
            w_z += 2 * a_s * ((i + 1) * dw) ** 2

        self._y_e = y_e
        self._i = i_panel
        self._a = a_tot
        self._z_y = z_y
        self._w_p = w_p
        self._w_z = w_z

    def hold_minimize_panel_setion(self, sigma_y, bc, sigma_x, p_lat, contourline, freq, nwdir):

        try:
            for self._t_lf in np.arange(0, 0.01, 0.001):
                for self._t_lfst in np.arange(0, 0.01, 0.001):
                    for self._h_lfst in np.arange(0, 2.0, 0.01):
                        self.init_cross_section()
                        if self.dynamic_panel_utilization(sigma_y, bc, sigma_x, p_lat, contourline, freq, nwdir) < 1:
                            print(self._t_lf)
                            print(self._t_lfst)
                            print(self._h_lfst)
                            raise BreakIt
        except BreakIt:
            pass

        return self.dynamic_panel_utilization(sigma_y, bc, sigma_x, p_lat, contourline, freq, nwdir)

    def hold2_minimize_panel_setion(self, sigma_y, bc, sigma_x, p_lat, contourline, freq, nwdir):
        def objective_function(x):
            # x = np.zeros(3)

            self._t_lf = x[0]  # Thickness of plate
            self._t_lfst = x[1]  # Width of stiffener
            self._h_lfst = x[2]  # Height of stiffener
            self.init_cross_section()
            y = np.abs(self.dynamic_panel_utilization(sigma_y, bc, sigma_x, p_lat, contourline, freq, nwdir) - 1)
            # print('{x[0]: 8.5f} {x[1]: 8.5f} {x[2]: 8.5f} {y: 8.5f}'.format(x=x,y = y))
            return y

        x0 = np.array([self._t_lf, self._t_lfst, self._h_lfst], dtype=float)
        x0 = np.array([0.02, 0.02, 0.2], dtype=float)
        bounds = Bounds([0.02, 0.02, 0.2], [0.1, 0.1, 6])

        res = minimize(objective_function, x0, method='tnc', options={'gtol': 1e-2}, bounds=bounds)
        #print(res)
        self._t_lf = res.x[0]
        self._t_lfst = res.x[1]
        self._h_lfst = res.x[2]
        self.init_cross_section()

        return self.dynamic_panel_utilization(sigma_y, bc, sigma_x, p_lat, contourline, freq, nwdir)

    def minimize_panel_setion(self, sigma_y, bc, sigma_x, p_lat, contourline, freq, nwdir):
        ones=np.ones(3,dtype=float)

        def objective_function(x):
            # x = np.zeros(3)

            self._t_lf = x[0]  # Thickness of plate
            self._t_lfst = x[1]  # Width of stiffener
            self._h_lfst = x[2]  # Height of stiffener
            self.init_cross_section()
            y = np.abs(self.dynamic_panel_utilization(sigma_y, bc, sigma_x, p_lat, contourline, freq, nwdir) - 1)
            # print('{x[0]: 8.5f} {x[1]: 8.5f} {x[2]: 8.5f} {y: 8.5f}'.format(x=x,y = y))
            return [0.1*y,0.5*y,y]

        x0 = np.array([self._t_lf, self._t_lfst, self._h_lfst], dtype=float)
        x0 = np.array([0.02, 0.02, 0.2], dtype=float)
        bounds = Bounds([0.02, 0.02, 0.2], [0.1, 0.1, 6])

        res = newton_krylov(objective_function,xin=[0.02, 0.02, 0.2])
        #print(res)
        self._t_lf = res[0]
        self._t_lfst = res[1]
        self._h_lfst = res[2]
        self.init_cross_section()

        return self.dynamic_panel_utilization(sigma_y, bc, sigma_x, p_lat, contourline, freq, nwdir)


    def dynamic_panel_utilization(self, sigma_y, bc, sigma_x, p_lat, contourline, freq, nwdir):
        # From Ultimate Load Analysis of Marine Structures
        # Tore H. S�reide
        # Section 5.5 Beam-Columns with no torsional buckling

        nfreq = freq.shape[0]

        e = self._settings.emod_st
        l = self._l
        l_e = self._l  # TODO: Check setting effective buckling length equal system length
        a = self._a
        i = self._i

        lambda_0 = l_e / np.sqrt(i / a)
        lambda_y = np.pi * np.sqrt(e / sigma_y)
        lambda_red = lambda_0 / lambda_y

        # Using buckling curve from EN 1993.1.1:2005 section  6.3.1.2:
        # A     section area
        # f_y   yield strength
        # n_cr  elastic critical force for the relevant buckling mode based on the gross sectional properties
        # alpha imperfection factor
        alpha = {
            'a0': 0.13,
            'a': 0.21,
            'b': 0.34,
            'c': 0.49,
            'd': 0.76
        }
        curve_name = 'c'
        phi = 0.5 * (1 + alpha[curve_name] * (lambda_red - 0.2) + lambda_red ** 2)
        chi = 1 / (phi + np.sqrt(phi ** 2 - lambda_red ** 2))

        # eq. 5.147
        kappa = chi  # Different notation between EN 1993 and Tore
        sigma_k = sigma_y * kappa

        # --------------------------------------------------
        #   Loop trough all frequencies and compute the
        #   interaction function for each combination of
        #   sigma_x and p_lat
        # -------------------------------------------------
        int_for = np.zeros([nfreq, nwdir], dtype='complex')
        elm_contour =np.zeros([nwdir,len(contourline)], dtype='float')
        for iwdir in range(nwdir):
            for ifreq in range(nfreq):


                # eq. 5.153 or 5.154, also table 5.21 b) and e)
                # Equivalent moments due to evenly distributed loads
                q = p_lat[ifreq, iwdir] / self._w_p
                cxm = None
                if bc == 'fixed':
                    cxm = 0.85 * q * l ** 2 / 16
                elif bc == 'pinned':
                    cxm = q * l ** 2 / 8
                else:
                    print('No such boundary condition: {}'.format(bc))
                    exit()

                # eq. 5.155
                m_p = sigma_y * self._z_y

                # eq. 5.156
                p_k = sigma_k * a
                sigma_e = np.pi ** 2 * e / lambda_0 ** 2
                p_e = sigma_e * a
                p = sigma_x[ifreq, iwdir] * a

                int_for[ifreq,iwdir] = p / p_k + cxm / ((1 - p / p_e) * m_p)

            # Now get expected max for each direction
            elm_contour[iwdir,:] = np.array([stwcl.expected_largest_maximum(int_for[:,iwdir], freq) for stwcl in contourline])

        # for item in x:
        #    print(item)

        return elm_contour.flatten().max()

    def axial_stress(self, f, pos_y_side=None):
        if pos_y_side == None:
            pos_y_side = True

        # Compression is positive
        a = self._a
        wz = self._w_z
        h = self._h

        # Get forces at section center
        # TODO: Change z to section center, now at waterline

        def sig(f):
            if f.ndim == 3:
                sig_ax = -f[:, :, 0] / (2 * a)
                # Simplified, assuming neutral axis at center
                sig_by = f[:, :, 4] / (h * a)
                sig_bz = f[:, :, 5] / (2 * wz)
            else:
                sig_ax = -f[0] / (2 * a)
                sig_by = f[4] / (h * a)
                sig_bz = f[5] / (2 * wz)
            if pos_y_side:
                return sig_ax + sig_by + sig_bz
            else:
                return sig_ax + sig_by - sig_bz

        return sig(f)

    def lateral_pressure(self, f):
        a = self._l * self._w_p
        if f.ndim == 3:
            return f[:, :, 2] / a
        else:
            return f[2] / a


# class EN_1993_1_1_2005():
#     def __init__(self, settings):
#         self._h = self._settings.job_data['floater']['Radial']['Heigth']
#         self._d_rc = self._settings.job_data['floater']['Radial']['Column']['Diameter']
#         self._d_cc = self._settings.job_data['floater']['Central column diameter']
#         self._t_lf = self._settings.job_data['floater']['Radial']['Flange']['Lower']['Plate']['Thickness']
#         self._t_uf = self._settings.job_data['floater']['Radial']['Flange']['Upper']['Plate']['Thickness']
#
#         self._t_lfst = self._settings.job_data['floater']['Radial']['Flange']['Lower']['Stiffener']['Longitudinal'][
#             'Thickness']
#         self._h_lfst = self._settings.job_data['floater']['Radial']['Flange']['Lower']['Stiffener']['Longitudinal'][
#             'Height']
#         self._t_ufst = self._settings.job_data['floater']['Radial']['Flange']['Upper']['Stiffener']['Longitudinal'][
#             'Thickness']
#         self._h_ufst = self._settings.job_data['floater']['Radial']['Flange']['Upper']['Stiffener']['Longitudinal'][
#             'Height']
#         self._n_lfst = self._settings.job_data['floater']['Radial']['Flange']['Lower']['Stiffener']['Longitudinal'][
#             'Number of']
#         self._gaf = self._settings.job_data['floater']['Gap factor']
#
#     def interaction_factors(self, section_class):
#
#         my_y = (1 - n_ed / n_cr_y) / (1 - chi_y * n_ed / n_cr_y)
#         my_z = (1 - n_ed / n_cr_z) / (1 - chi_z * n_ed / n_cr_z)
#
#         w_y = w_pl_y / w_el_y
#         if w_y > 1.5:
#             w_y = 1.5
#
#
#         w_z = w_pl_z / w_el_z
#         if w_z > 1.5:
#             w_z = 1.5
#
#         eta_pl=n_ed/(n_rk/gamma_m1)
#
#         c_my =
#
#
#         kyy = c_my * c_mlt * my_y / (1 - n_ed / n_cr_y)
#         kyz = c_mz * my_y / (1 - n_ed / n_cr_z)
#         kzy = c_my * c_mlt * my_z / (1 - n_ed / n_cr_y)
#         kzz = c_mz * my_z / (1 - n_ed / n_cr_z)
#         if section_class == 3 or section_class == 4:
#             kyy = kyy / c_yy
#             kyz = kyz / c_yz * 0.6 * np.sqrt(wz / wy)
#             kzy = kzy / c_zy * 0.6 * np.sqrt(wy / wz)
#             kzz = kzz / c_zz
#
#         return k_yy, k_yz, k_zy, k_zz


class DNVGL_RP_C201_Part1():  # NOT APPLICABLE, KEPT FOR FUTURE REFERENCE
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

    # Section 7.2 Forces in the idealised stiffened plate
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
        # Return the equivalent lateral line load: q_Sd
        s = self._d_rc / (self._n_lfst - 1)  # spacing between stiffeners

        # eq. 7.9
        # No stress in transverse direction
        p0 = 0

        return (p_Sd + p0) * s

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
        w_es = i_es / y_

    # Section 7.3 Effective plate width
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

    # Section 7.4 Resistance of plate between stiffeners
    def eq_7_18(self, tau_Sd, f_y, gamma_m):
        return tau_Sd / (f_y / (1.73205080757 * gamma_m))

    # 7.5 Characteristic buckling strength of stiffeners
