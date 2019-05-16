__author__ = "Eivind Sonju"
__copyright__ = "Copyright (C) 2017-2019 Verbun AS. All rights reserved."
__version__ = "2.0"


import numpy as np

from tool_box import *


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

    def set_new_reduction_point(self, new_point):
        # Recursively update reduction point on myself and my children
        assert len(new_point) == 3
        self._inertias._point += new_point
        if self.parts_list:
            for part in self.parts_list:
                part.set_new_reduction_point(new_point)

    def get_parts(self, mass_list=None):
        # Only get mass of parts that have no part
        if self.parts_list == []:
            mass_list.append(self._inertias)
        if mass_list == None:
            mass_list=[]
        for part in self.parts_list:
            part.get_parts(mass_list)
        return mass_list




    #     def rec(x):
    #         if self.parts_list:
    #             for part in self.parts_list:
    #                 x.append(part.rec(x))
    #         else:
    #             return self._inertias
    #
    #     return rec(x)
    #
    # def pop_list(nodes=None, parent=None, node_list=None):
    #     if parent is None:
    #         return node_list
    #     node_list.append([])
    #     for node in nodes:
    #         if node['parent'] == parent:
    #             node_list[-1].append(node)
    #         if node['id'] == parent:
    #             next_parent = node['parent']
    #
    #     pop_list(nodes, next_parent, node_list)
    #     return node_list

    @property
    def inertias(self):
        return self._inertias

    @property
    def mass(self):
        return self._inertias.mass

    @property
    def inertia_matrix(self):
        return self._inertias.inertia_matrix_global

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
                    s = s + '{: 10.3g}'.format(np.round_(m[i, j], decimals=0))
                s = s + '\n'
            s = s + '\n'
            print(s)


class AssemblyClass(ModelClass, object):
    def __init__(self, type=None, red_point=None):
        ModelClass.__init__(self, type, red_point)

    def aggregate_inertias_from_parts(self, parts_list):
        for part in parts_list:
            # TODO: Revise mass and mass matrix
            self.inertias.add_mass_matrix_local(part.inertias.mass_matrix_global)

        for part in parts_list:
            self.inertias._point += part.inertias._point * part.inertias.mass / self.inertias.mass

        # self._inertias.update_vector_matrix_global()
        self._inertias.update_mass_matrix_global()
        self.print_vector_matrix_global()


class UnitClass(AssemblyClass, object):
    def __init__(self, models):
        ModelClass.__init__(self)

        self._models = models
        self.__update_global__()
        # self.print_vector_matrix_global()

        # self._inertias.reduction_point = self._red_point
        # self.__set_cog_relative_to_point__()

    def __update_global__(self):
        for model in self._models:
            self.parts_list.append(model)
        self.aggregate_inertias_from_parts(self.parts_list)


class WtgClass(AssemblyClass, object):
    def __init__(self, twr_data, rna_data, rho_st=7850):
        ModelClass.__init__(self)
        #
        self._twr_data = twr_data
        self._rna_data = rna_data
        #
        self._rho_st = rho_st
        # self.parts_list = []
        self.__update_global__()

        # self.print_vector_matrix_global()

    def __update_global__(self):
        self.parts_list.append(RNAClass(self._rna_data))
        self.parts_list.append(TowerClass(self._twr_data, self._rho_st))
        # for part in self.parts_list:
        #     print(part._inertias.inertia_matrix_global)
        self.aggregate_inertias_from_parts(self.parts_list)


class FloaterClass(AssemblyClass, object):
    # CoGz will be set for BOS at z=0
    def __init__(self, floater_data, rho_st):
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
        self._rho_st = rho_st
        self._rho_bal = floater_data.rho_bal
        # self._l_radial = 0
        # self._m_rc = 0
        # self._m_hc = 0
        # self._m_lf = 0
        # self._m_uf = 0
        # self._m_radial = 0
        # self._m_hub = 0
        # self._m_ballast = 0
        # self._m_global = 0
        # self._w_lf = 0
        # self._w_uf = 0
        # self.parts_list = []

        self.__update_global__()

        # self.print_vector_matrix_global()

        # self._inertias.reduction_point = self._red_point
        # self.__set_cog_relative_to_point__()

    def __update_global__(self):
        self._inertias.reset()

        da = (1 + self._gap) * self._dia_rc
        dtheta = 2 * pi / self._nr
        theta = [i * dtheta for i in range(self._nr)]
        zr = -(self._hgt / 2 + self._t_lf)
        hf_hc = self._ballast_filling[0]
        zhc_bal = -(hf_hc / 2 + self._t_lf)

        self.parts_list.append(HubColumnClass(type='Hub column cylinder',
                                              dia_hc=self._dia_hc,
                                              thi_hc=self._thi_hc,
                                              hgt=self._hgt,
                                              rho_st=self._rho_st,
                                              red_point=[0, 0, zr]))
        self.parts_list.append(BallastClass(type='Hub column ballast',
                                            irow=-1,
                                            icol=-1,
                                            dia_bal=self._dia_rc - 2 * self._thi_rc,
                                            rc_internal_hgt=self._hgt,
                                            filling=hf_hc,
                                            rho_bal=self._rho_bal,
                                            red_point=[0, 0, zhc_bal]))
        for ir in range(self._nr):
            dxc = cos(theta[ir]) * da
            dyc = sin(theta[ir]) * da
            # Reduction point is set at center bottom of steel for all parts.
            # Flanges
            # TODO: Discretize flanges every meter or so in radial direction for better mass resolution
            self.parts_list.append(FlangeClass(type='Upper flange',
                                               irow=ir,
                                               icol=None,
                                               a=da * self._nc,
                                               b=self._dia_rc,
                                               h=self._t_lf,
                                               density=self._rho_st,
                                               theta=theta[ir],
                                               red_point=[-dxc * self._nc / 2, -dyc * self._nc / 2,
                                                          -(self._t_lf + self._hgt + self._t_uf / 2)]))
            self.parts_list.append(FlangeClass(type='Lower flange',
                                               irow=ir,
                                               icol=None,
                                               a=da * self._nc,
                                               b=self._dia_rc,
                                               h=self._t_uf,
                                               density=self._rho_st,
                                               theta=theta[ir],
                                               red_point=[-dxc * self._nc / 2, -dyc * self._nc / 2, -self._t_lf / 2]))
            # Radial columns
            for ic in range(self._nc):
                hf_rc = self._ballast_filling[1][ic]

                xr = -(ic + 1) * dxc
                yr = -(ic + 1) * dyc
                zrc_bal = -(hf_rc / 2 + self._t_lf)

                self.parts_list.append(RadialColumnClass(type='Radial column cylinder',
                                                         irow=ir,
                                                         icol=ic,
                                                         dia_rc=self._dia_rc,
                                                         thi_rc=self._thi_rc,
                                                         hgt=self._hgt,
                                                         rho_st=self._rho_st,
                                                         red_point=[xr, yr, zr]))
                self.parts_list.append(BallastClass(type='Radial column ballast',
                                                    irow=ir,
                                                    icol=ic,
                                                    dia_bal=self._dia_rc - 2 * self._thi_rc,
                                                    rc_internal_hgt=self._hgt,
                                                    filling=hf_rc,
                                                    rho_bal=self._rho_bal,
                                                    red_point=[xr, yr, zrc_bal]))

        # Set inertias and CoG relative to bottom of tower
        self.aggregate_inertias_from_parts(self.parts_list)

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
        self.nc = fdi['Number of radial columns']  # 2, 3 or 4
        self.gap = fdi['Gap factor']
        self.dia_rc = fdi['Radial column diameter']
        self.thi_rc = fdi['Radial column thickness']
        self.dia_hc = fdi['Central column diameter']
        self.thi_hc = fdi['Central column thickness']
        self.hgt = fdi['Radial height']
        self.t_lf = fdi['Lower flange thickness']
        self.t_uf = fdi['Upper flange thickness']
        self._ballast_filling = fdi['Ballast filling ratio']
        self.rho_bal = fdi['Ballast density']
