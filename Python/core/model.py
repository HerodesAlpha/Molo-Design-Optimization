__author__ = "Eivind Sonju"
__copyright__ = "Copyright (C) 2017-2019 Verbun AS. All rights reserved."
__version__ = "2.0"

from core.tool_box import *


class ModelClass(object):
    def __init__(self, type=None, red_point=None):
        self._inertias = TotalMassMatrixClass()
        self.verbose = 1
        self.parts_list = []
        if type is None:
            self._type = ''
        else:
            self._type = type

        if red_point is None:
            self._inertias._point = np.zeros(3)
        else:
            assert (len(red_point) == 3)
            self._inertias._point = np.asarray(red_point)

    # def set_reduction_point(self, vector):
    #     # Recursively update reduction point on myself and my children
    #     assert len(vector) == 3
    #     for part in self.get_all_parts():
    #         part._inertias.reduction_point = vector  # Will update global mass matrix also

    def move_reduction_point(self, vector):
        # Recursively update reduction point on myself and my children
        assert len(vector) == 3
        for part in self.get_all_parts():
            part._inertias.reduction_point = part._inertias._point + vector  # Will update global mass matrix also

    def get_parts_without_children(self, part_list=None):
        # Only get mass of parts that have no part
        if self.parts_list == []:
            part_list.append(self)
        if part_list == None:
            part_list = []
        for part in self.parts_list:
            part.get_parts_without_children(part_list)
        return part_list

    def get_all_parts(self, part_list=None):
        if part_list == None:
            part_list = []
        part_list.append(self)
        for part in self.parts_list:
            part.get_all_parts(part_list)
        return part_list

    @property
    def inertias(self):
        return self._inertias

    @property
    def mass(self):
        return self._inertias.mass

    @property
    def cog(self):
        return self._inertias.cog

    @cog.setter
    def cog(self, val):
        assert len(val) == 3
        self._inertias._cog = val

    # def __update_reference_point__(self):
    #    self._inertias.reduction_point = self._point
    #
    # def __set_cog_relative_to_point__(self):
    #     pass
    #     self._inertias._cog = self._inertias.cog - self._red_point

    def print_vector_matrix_global(self):
        if not self.verbose == 0:
            print(
                '{} {}\n\tMass = {m:6.1f} t\n\tReduction point is [{a[0]:6.2f}, {a[1]:6.2f}, {a[2]:6.2f}]\n\n\tVector matrix'.format(
                    self.__class__.__name__, self._type, m=self.inertias.mass / 1000,
                    a=self.inertias.reduction_point))
            m = self.inertias.mass_matrix_global / self.inertias.mass
            s = ''
            for i in range(6):
                for j in range(6):
                    s = s + '{: 10.3g}'.format(np.round(m[i, j], decimals=0))
                s = s + '\n'
            s = s + '\n'
            print(s)


class AssemblyClass(ModelClass, object):
    def __init__(self, type=None, red_point=None):
        ModelClass.__init__(self, type, red_point)

    def aggregate_inertias_from_parts(self, parts_list):
        # Check that mass matrix is zero
        sum_mass_matrix_global = np.zeros((6, 6), dtype='float')
        for part in parts_list:
            sum_mass_matrix_global += part.inertias.mass_matrix_global

        # Calculate origin relative to CoG
        self.inertias._point = np.zeros(3)
        for part in parts_list:  # Calculate reference point for this assembly
            self.inertias._point += part.inertias._point * part.mass
        self.inertias._point /= sum_mass_matrix_global[0, 0]

        self.inertias._mass_matrix_local = sum_mass_matrix_global - self.inertias._huygens_transport() * \
                                           sum_mass_matrix_global[0, 0]
        self._inertias.cog = np.zeros(3)  # Just to be sure
        # self.print_vector_matrix_global()


class UnitClass(AssemblyClass, object):
    # The mass matrix of all parts, parents and children are defined relative to the plane xy coinciding with the calm
    # water level plane. z - axis points upwards. Draught is zero.
    def __init__(self, models):
        ModelClass.__init__(self)

        self._models = models
        for model in self._models:
            self.parts_list.append(model)
        self.aggregate_inertias_from_parts(self.parts_list)


class WtgClass(AssemblyClass, object):
    def __init__(self, twr_data, rna_data, settings):
        ModelClass.__init__(self)
        #
        self._twr_data = twr_data
        self._rna_data = rna_data
        #
        self._rho_st = settings.rho_st

        self.parts_list.append(RNAClass(self._rna_data))
        self.parts_list.append(TowerClass(self._twr_data, self._rho_st))
        self.aggregate_inertias_from_parts(self.parts_list)


class FloaterClass(AssemblyClass, object):
    # CoGz will be set for BOS at z=0
    def __init__(self, floater_data, settings):
        ModelClass.__init__(self)
        self._nr = 3
        self._nc = floater_data.nc
        self._gap = floater_data.gap
        self._dia_rc = floater_data.dia_rc
        self._thi_rc = floater_data.thi_rc
        self._dia_hc = floater_data.dia_hc
        self._thi_hc = floater_data.thi_hc
        self._hgt = floater_data.hgt
        self._t_lf = floater_data.t_lf
        self._t_uf = floater_data.t_uf
        self._ballast_filling = floater_data._ballast_filling
        self._rho_st = settings.rho_st
        self._rho_bal = floater_data.rho_bal
        self._t_lfst = floater_data.t_lfst
        self._h_lfst = floater_data.h_lfst
        self._t_ufst = floater_data.t_ufst
        self._h_ufst = floater_data.h_ufst

        self._n_strips = 200  # Number of flange strips in longitudinal direction

        self._w_lf = floater_data.w_lf
        print('Width of lower flange {:1.2f}'.format(self._w_lf))
        self._w_uf = self._dia_rc
        print('Width of upper flange {:1.2f}'.format(self._w_uf))

        self._l_uf = (1 + self._gap) * self._dia_rc * self._nc + 0.5 * self._dia_rc
        print('Length of upper flange {:1.2f}'.format(self._l_uf))
        self._l_lf = self._l_uf + floater_data.l_lf_overlength
        print('Length of lower flange {:1.2f}'.format(self._l_lf))

        self._draught = settings.draught

        self._generate_parts()

        self.aggregate_inertias_from_parts(self.parts_list)

        pass
        # self.print_vector_matrix_global()

    def _generate_parts(self):
        # self._inertias.reset()

        da = (1 + self._gap) * self._dia_rc
        lower_flange_strip_width = self._l_lf / self._n_strips
        upper_flange_strip_width = self._l_uf / self._n_strips
        dtheta = 2 * pi / self._nr
        theta = [i * dtheta for i in range(self._nr)]
        zr = -(self._hgt / 2 + self._t_lf)
        hf_hc = self._ballast_filling
        zhc_bal = -1*(hf_hc * self._hgt / 2 + self._t_lf)

        def func_rpx(i):
            dx = - da * self._nc / self._n_strips

            def rpx(i):
                return dx * (i + 1 / 2)

            return rpx(i)

        rpx = func_rpx

        rpy = 0
        rpz_uf = -(self._t_lf + self._hgt + self._t_uf / 2)
        rpz_lf = -self._t_lf / 2
        dy_ufst = self._w_uf / 2 - self._t_ufst / 2
        dz_ufst = self._t_uf / 2 + self._h_ufst / 2
        dy_lfst = self._w_lf / 2 - self._t_lfst / 2
        dz_lfst = -self._t_lf / 2 - self._h_lfst / 2
        lr = ['left', 'right']

        self.parts_list.append(HubColumnClass(type='Hub column cylinder',
                                              dia_hc=self._dia_hc,
                                              thi_hc=self._thi_hc,
                                              hgt=self._hgt,
                                              rho_st=self._rho_st,
                                              red_point=[0, 0, zr + self._draught]))
        self.parts_list.append(BallastClass(type='Hub column ballast',
                                            irow=-1,
                                            icol=-1,
                                            dia_bal=self._dia_rc - 2 * self._thi_rc,
                                            rc_internal_hgt=self._hgt,
                                            filling=hf_hc,
                                            rho_bal=self._rho_bal,
                                            red_point=[0, 0, zhc_bal + self._draught]))

        for ir in range(self._nr):
            # dxc = cos(theta[ir]) * da
            # dyc = sin(theta[ir]) * da
            rot_mat = rotation_matrix([0, 0, theta[ir]])

            # Reduction point is set at center bottom of steel for all parts.
            # Flanges
            # Discretize flanges every meter or so in radial direction for better mass resolution
            for istrip in range(self._n_strips):
                self.parts_list.append(
                    FlangeClass(type='Radial {r:1.0f}, upper flange, strip {s:1.0f}'.format(r=ir + 1, s=istrip + 1),
                                irow=ir,
                                icol=None,
                                a=upper_flange_strip_width,
                                b=self._w_uf,
                                h=self._t_uf,
                                density=self._rho_st,
                                theta=theta[ir],
                                red_point=rot_mat @ [rpx(istrip), rpy, rpz_uf + self._draught]))

                # for ist in range(2):
                #     self.parts_list.append(FlangeClass(
                #         type='Radial {r:1.0f}, upper flange, {a} stiffener, strip {s:1.0f}'.format(r=ir + 1,
                #                                                                                    a=lr[ist],
                #                                                                                    s=istrip + 1),
                #         irow=ir,
                #         icol=None,
                #         a=upper_flange_strip_width,
                #         b=self._t_ufst,
                #         h=self._h_ufst,
                #         density=self._rho_st,
                #         theta=theta[ir],
                #         red_point=rot_mat @ [rpx(istrip), rpy + (-1) ** ist * dy_ufst,
                #                              rpz_uf + dz_ufst + self._draught]))

                self.parts_list.append(
                    FlangeClass(type='Radial {r:1.0f}, lower flange, strip {s:1.0f}'.format(r=ir + 1, s=istrip + 1),
                                irow=ir,
                                icol=None,
                                a=lower_flange_strip_width,
                                b=self._w_lf,
                                h=self._t_lf,
                                density=self._rho_st,
                                theta=theta[ir],
                                red_point=rot_mat @ [rpx(istrip), rpy, rpz_lf + self._draught]))

                # for ist in range(2):
                #     self.parts_list.append(FlangeClass(
                #         type='Radial {r:1.0f}, lower flange, {a} stiffener, strip {s:1.0f}'.format(r=ir + 1,
                #                                                                                    a=lr[ist],
                #                                                                                    s=istrip + 1),
                #         irow=ir,
                #         icol=None,
                #         a=lower_flange_strip_width,
                #         b=self._t_ufst,
                #         h=self._h_ufst,
                #         density=self._rho_st,
                #         theta=theta[ir],
                #         red_point=rot_mat @ [rpx(istrip), rpy + (-1) ** ist * dy_lfst,
                #                              rpz_lf + dz_lfst + self._draught]))

            # Radial columns
            for ic in range(self._nc):
                hf_rc = self._ballast_filling

                xr = -(ic + 1) * da
                yr = 0
                zrc_bal = -(hf_rc * self._hgt / 2 + self._t_lf)

                self.parts_list.append(RadialColumnClass(type='Radial column cylinder',
                                                         irow=ir,
                                                         icol=ic,
                                                         dia_rc=self._dia_rc,
                                                         thi_rc=self._thi_rc,
                                                         hgt=self._hgt,
                                                         rho_st=self._rho_st,
                                                         red_point=rot_mat @ [xr, yr, zr + self._draught]))

                self.parts_list.append(BallastClass(type='Radial column ballast',
                                                    irow=ir,
                                                    icol=ic,
                                                    dia_bal=self._dia_rc - 2 * self._thi_rc,
                                                    rc_internal_hgt=self._hgt,
                                                    filling=hf_rc,
                                                    rho_bal=self._rho_bal,
                                                    red_point=rot_mat @ [xr, yr, zrc_bal + self._draught]))

    @property
    def nc(self):
        return self._nc

    @property
    def gap(self):
        return self._gap

    @property
    def dia_rc(self):
        return self._dia_rc

    @property
    def dia_hc(self):
        return self._dia_hc

    @property
    def thi_hc(self):
        return self._thi_rc

    @property
    def hgt(self):
        return self._hgt

    @property
    def t_lf(self):
        return self._t_lf

    @property
    def t_uf(self):
        return self._t_uf


class RadialColumnClass(ModelClass, object):
    def __init__(self, type, irow, icol, dia_rc, thi_rc, hgt, rho_st, red_point=None):
        ModelClass.__init__(self, type, red_point)

        self._irow = irow
        self._icol = icol
        self._dia_rc = dia_rc
        self._thi_rc = thi_rc
        self._hgt = hgt
        self._rho_st = rho_st
        self.__update_inertias__()
        # self._inertias.reduction_point = self._red_point
        # print(self._red_point)
        # self.__set_cog_relative_to_point__()
        # self.print_vector_matrix_global()

    def __update_inertias__(self):
        hollow_right_circular_cylinder(self._inertias,
                                       int_radius=self._dia_rc / 2 - self._thi_rc,
                                       ext_radius=self._dia_rc / 2,
                                       length=self._hgt,
                                       density=self._rho_st)


class HubColumnClass(ModelClass, object):
    def __init__(self, type, dia_hc, thi_hc, hgt, rho_st, red_point=None):
        ModelClass.__init__(self, type, red_point)

        self._dia_rc = dia_hc
        self._thi_rc = thi_hc
        self._hgt = hgt
        self._rho_st = rho_st
        self.__update_inertias__()
        # self._inertias.reduction_point = self._red_point
        # print(self._red_point)
        # self.__set_cog_relative_to_point__()
        # self.print_vector_matrix_global()

    def __update_inertias__(self):
        hollow_right_circular_cylinder(self._inertias,
                                       int_radius=self._dia_rc / 2 - self._thi_rc,
                                       ext_radius=self._dia_rc / 2,
                                       length=self._hgt,
                                       density=self._rho_st)


class BallastClass(ModelClass, object):
    def __init__(self, type, irow, icol, dia_bal, rc_internal_hgt, filling, rho_bal, red_point=None):
        ModelClass.__init__(self, type, red_point)

        self._irow = irow
        self._icol = icol
        self._dia_bal = dia_bal
        self._rc_internal_hgt = rc_internal_hgt

        self._filling = filling
        self._rho_bal = rho_bal

        self.__update_inertias__()
        # self.print_vector_matrix_global()

    def __update_inertias__(self):
        d = self._dia_bal
        h = self._rc_internal_hgt
        vol = np.pi * d ** 2 / 4 * h * self._filling
        self.inertias.mass = vol * self._rho_bal


class FlangeClass(ModelClass, object):
    def __init__(self, type, irow, icol, a, b, h, density, theta, red_point=None):
        ModelClass.__init__(self, type, red_point)

        self._irow = irow
        self._icol = icol
        self._a = a
        self._b = b
        self._h = h
        self._theta = theta
        self._density = density

        self.__update_inertias__()
        # self._inertias.reduction_point = self._red_point
        # self.__set_cog_relative_to_point__()
        # self.print_vector_matrix_global()

    def __update_inertias__(self):
        rectangular_prism(self._inertias, a=self._a, b=self._b, h=self._h, density=self._density)
        self.inertias.rotate_z_inertia_mass_matrix_local(self._theta)


class RNAClass(ModelClass, object):
    def __init__(self, rna_data):
        ModelClass.__init__(self, rna_data.type, rna_data.p)
        self._rna_data = rna_data

        self.__update_inertias__()
        self._inertias.reduction_point = self._rna_data.p
        # self.__set_cog_relative_to_point__()
        # self.print_vector_matrix_global()

    def __update_inertias__(self):
        self.inertias.reset()
        self.inertias.mass = self._rna_data.m
        self.inertias.ixx = self._rna_data.ixx
        self.inertias.iyy = self._rna_data.iyy
        self.inertias.izz = self._rna_data.izz
        self.inertias.iyz = self._rna_data.iyz
        self.inertias.ixz = self._rna_data.ixz
        self.inertias.ixy = self._rna_data.ixy


class RNADataClass():
    def __init__(self, rdi):
        self.type = rdi['Type']
        self.m = rdi['Mass']
        self.d_rot = rdi['Rotor diameter']
        self.ixx = rdi['Ixx']
        self.iyy = rdi['Iyy']
        self.izz = rdi['Izz']
        self.iyz = rdi['Iyz']
        self.ixz = rdi['Iyz']
        self.ixy = rdi['Ixy']
        self.p = [0, 0, 0]


class TowerClass(ModelClass, object):
    def __init__(self, twr_data, rho_st):
        ModelClass.__init__(self, twr_data.type, twr_data.p)
        self._twr_data = twr_data
        self._rho_st = rho_st

        self.__update_inertias__()
        # print(self._inertias._3d_rotational_inertia)

        # self._inertias.reduction_point = self._twr_data.p
        # print(self._inertias._3d_rotational_inertia)
        # self.__set_cog_relative_to_point__()
        # self.print_vector_matrix_global()

    def __update_inertias__(self):
        hollow_right_circular_cylinder(self._inertias, int_radius=self._twr_data.d / 2 - self._twr_data.t,
                                       ext_radius=self._twr_data.d / 2,
                                       length=self._twr_data.h,
                                       density=self._rho_st)


class TowerDataClass():
    def __init__(self, tdi):
        self.type = tdi['Type']
        self.h = tdi['Height']
        self.t = tdi['Thickness']
        self.d = tdi['Diameter']
        self.p = [0, 0, 0]


class FloaterDataClass():
    def __init__(self, fdi):
        self.type = fdi['Type']
        self.nc = fdi['Radial']['Number of columns']  # 2, 3 or 4
        self.gap = fdi['Gap factor']
        self.dia_rc = fdi['Radial']['Column']['Diameter']
        self.thi_rc = fdi['Radial']['Column']['Thickness']
        self.dia_hc = fdi['Central column diameter']
        self.thi_hc = fdi['Central column thickness']
        self.hgt = fdi['Radial']['Heigth']
        self.t_lf = fdi['Radial']['Flange']['Lower']['eq_thick']
        self.w_lf = fdi['Radial']['Flange']['Lower']['Plate']['Width']
        self.t_uf = fdi['Radial']['Flange']['Upper']['eq_thick']
        self._ballast_filling = fdi['Ballast filling ratio']
        self.rho_bal = fdi['Ballast density']
        self.t_lfst = fdi['Radial']['Flange']['Lower']['Stiffener']['Longitudinal']['Thickness']
        self.h_lfst = fdi['Radial']['Flange']['Lower']['Stiffener']['Longitudinal']['Height']
        self.t_ufst = fdi['Radial']['Flange']['Upper']['Stiffener']['Longitudinal']['Thickness']
        self.h_ufst = fdi['Radial']['Flange']['Upper']['Stiffener']['Longitudinal']['Height']
        self.l_lf_overlength = fdi['Radial']['Flange']['Lower']['Overlength']
