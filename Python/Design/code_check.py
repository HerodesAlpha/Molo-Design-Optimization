import numpy as np


# Code check of t

# Tripping of stiffener, use stiffened plate criteria?

class panel():
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

        self._l = self._d_rc*self._gaf

        self._a = 0
        self._y = 0
        self._i = 0

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

        a_p = w_p * t_s
        a_s = h_s * t_s
        y_ = ((a_p * t_p / 2) + n_s * (a_s * (h_s / 2 + t_p))) / (a_p + n_s * a_p)
        i_p = i(w_p, t_p)
        i_s = i(t_s, h_s)
        i_panel = i_p + a_p * (y_ - t_p / 2)**2 + i_s + a_s * (y_ - (t_p + h_s / 2))**2
        self._y = y_
        self._i = i_panel
        self._a = a_p + n_s * a_s

    def utilization(self,sigma_y):
        # From Ultimate Load Analysis of Marine Structures
        # Tore H. Søreide
        # Section 5.5 Beam-Columns with no torsional buckling

        e = self._settings.emod_st
        l_e = self._l # TODO: Check setting effective buckling length equal system length
        lambda_0 = l_e/np.sqrt(self._i / self._a)
        lambda_y = np.pi*np.sqrt(e/sigma_y)

        lambda_red = lambda_0/lambda_y

        # Using buckling curve from EN 1993.1.1:2005 section  6.3.1.2:
        # A     section area
        # f_y   yield strength
        # n_cr  elastic critical force for the relevant buckling mode based on the gross sectional properties
        # alpha imperfection factor
        alpha = {
                'a0': 0.13,
                'a' : 0.21,
                'b' : 0.34,
                'c' : 0.49,
                'd' : 0.76
        }
        curve_name='c'
        phi = 0.5 * (1 + alpha[curve_name] * (lambda_red - 0.2) + lambda_red ** 2)
        chi = 1 / (phi + np.sqrt(phi ** 2 - lambda_red ** 2))

        # eq. 5.155
        m_p = sigma_y * w_p

        # eq. 5.156
        interaction_formula = p / p_k + c * m / ((1 - p / p_e) * m_p)
        return interaction_formula


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
