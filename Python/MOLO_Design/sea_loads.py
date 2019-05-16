__author__ = "Eivind Sonju"
__copyright__ = "Copyright (C) 2017-2019 Verbun AS. All rights reserved."
__version__ = "2.0"

import pickle
from meshmagick.mesh import Mesh
from pylab import *
import helper
import nemoh
import tool_box as tb
from common import PhysicalQuantities


class HydroCoefficients(PhysicalQuantities, object):
    # TODO: Get added mass at zero and infinite frequency
    def __init__(self, dof, w, dir, fio, sym):
        super().__init__()

        self._w = w
        self._nw = len(w)
        self._dof = dof
        self._ndof = int(np.sum(dof))
        self._dir = dir
        self._ndir = len(dir)
        self._sym = sym

        step = 0
        sys.stdout.write("\nInit hydro:\n")
        # -------------------------------------------------------
        step += 1
        sys.stdout.write("\t({})Get panel data\n".format(step))
        self._pd = nemoh.PanelData(fio, self._sym)
        self._npanels = self._pd.ppanels.shape[0]
        # -------------------------------------------------------
        step += 1
        sys.stdout.write("\t({})Get forces\n".format(step))
        self._fe = nemoh.get_fe(fio, dof, w, dir)
        # -------------------------------------------------------
        step += 1
        sys.stdout.write("\t({})Get static mass and waterplane stiffness\n".format(step))
        self._m, self._k = tb.load_M_and_K(fio.data_io_dir)
        # -------------------------------------------------------
        step += 1
        sys.stdout.write("\t({})Get added mass and damping\n".format(step))
        self._ma, self._c_hyd = nemoh.get_ab(fio, dof, w, dir)
        # -------------------------------------------------------
        nproblems = (len(dir) + sum(dof)) * len(w)
        with open(fio.nemoh_root.joinpath('Normalvelocities.dat')) as f:
            if not int(f.readline()) == nproblems:
                print('nproblems mismatch')
                exit()
        # -------------------------------------------------------
        step += 1
        sys.stdout.write("\t({})Get pressures\n".format(step))
        pressure_file = fio.data_io_dir.joinpath('hydro_pressures.pkl')
        if pressure_file.is_file() and 1:
            self._p = pickle.load(open(pressure_file, 'rb'))
            sys.stdout.write("\t\tUnPickled from {}\n".format(pressure_file))
            for key in self._p.keys():
                print("\t\t\t{}".format(key))

        else:
            self._p = dict()
            sys.stdout.write("\t\tHydro static\n")
            self._p['Hydro static'] = (self._rho_sw * self._grav) * self._pd.ppanel_centers[:, 2]
            sys.stdout.write("\t\tFroude-Krylof \n")
            self._p['Froude-Krylof'] = nemoh.get_fk_pressure(fio, ndir=len(dir), nomega=self._nw,
                                                             npanels=self._pd.npanels)
            sys.stdout.write("\t\tDiffraction\n")
            self._p['Diffraction'] = np.zeros([self._nw, self._ndir, self._npanels], dtype=complex)
            for iw in range(self._nw):
                for ibeta in range(self._ndir):
                    pn = nemoh.diffraction_problem_number(iw, ibeta, self._ndir, self._ndof)
                    # print('\t\t\tProblem {}'.format(pn))
                    self._p['Diffraction'][iw, ibeta, :] = nemoh.get_poten_pressure(fio, npoints=self._pd.npoints,
                                                                                    ppanels=self._pd.ppanels,
                                                                                    problem_number=pn)
            sys.stdout.write("\t\tRadiation\n")
            self._p['Radiation'] = np.zeros([self._nw, self._ndof, self._npanels], dtype=complex)
            for iw in range(self._nw):
                for iradiation in range(self._ndof):
                    pn = nemoh.radiation_problem_number(iw, iradiation, self._ndir, self._ndof)
                    # print('\t\t\tProblem {}'.format(pn))
                    self._p['Radiation'][iw, iradiation, :] = nemoh.get_poten_pressure(fio, npoints=self._pd.npoints,
                                                                                       ppanels=self._pd.ppanels,
                                                                                       problem_number=pn)
            pickle.dump(self._p, open(pressure_file, "wb"))

        print('\n{} initialized\n'.format(self.__str__()))

    def show_pressure(self, ifreq, pressure_index, axis, pressure_type):
        xyz = np.zeros([self._pd._npanel, 3])
        xyz[:, axis] = 1
        nemoh_mesh = Mesh(self._pd.ppoints, self._pd.ppanels)
        vec = (self._pd.ppanel_normals * np.real(self._p[pressure_type][ifreq, pressure_index])[:, np.newaxis]) * xyz
        h = helper.show_force(nemoh_mesh, self._pd.ppanel_centers, vec)
        h.show()

    def p2f(self, ifreq, pressure_index, pressure_type):
        # Pressure to force
        f_normal = np.zeros((self._npanels), dtype=np.complex)
        f = np.zeros((self._npanels, 3), dtype=np.complex)
        for i, panel in enumerate(self._pd.ppanels):
            f_normal[i] = self._p[pressure_type][ifreq, pressure_index, i] * self._pd.ppanel_areas[i]
            for j in range(3):
                f[i, j] = -f_normal[i] * self._pd.ppanel_normals[i, j]
        return f

    @property
    def w(self):
        return self._w

    def sdof_val(self, m, dir_index, dof_index):
        return m[dir_index, :, dof_index]

    @property
    def available_dofs(self):
        return [i + 1 for i, e in enumerate(self._dof) if e != 0]

    @property
    def available_dirs(self):
        return [e for i, e in enumerate(self._dir)]

    def get_dof_index(self, sel_dof):
        dof_index = 0
        if self._dof[sel_dof - 1] == 1:
            return (sum(self._dof[:sel_dof]) - 1)
        else:
            print('No such DOF')

    #
    def get_dir_index(self, sel_dir):
        return int(sel_dir)

    @property
    def ma(self):
        return self._ma

    @property
    def m(self):
        return self._m

    @property
    def k(self):
        return self._k

    @property
    def fe(self):
        return self._fe

    @property
    def pd(self):
        return self._pd


def spec_wave(hs, wp, w, gamma=None):
    sig_a = 0.07
    sig_b = 0.09
    delta_sig = sig_b - sig_a

    tp = 2 * np.pi / wp
    x = tp / np.sqrt(hs)
    if gamma == None:
        if x <= 3.6:
            gamma = 5
        elif x < 5:
            gamma = np.exp(5.75 - 1.15 * x)
        else:
            gamma = 1
        # print('Gamma {}'.format(gamma))

    a_gamma = 1 - 0.287 * np.log(gamma)

    def spec_pm(w):
        return (5 / 16) * (hs ** 2) * (wp ** 4) * (w ** (-5)) * np.exp(-(5 / 4) * ((w / wp) ** (-4)))

    def spec_jonswap(w):
        def sig(w):
            return sig_a if w <= wp else sig_b

        sig_ab = np.array(list(map(sig, w)))

        return a_gamma * spec_pm(w) * gamma ** np.exp(-0.5 * ((w - wp) / sig_ab * wp))

    if not gamma == 1:
        spec = spec_jonswap
    else:
        spec = spec_pm

    return spec(w)


def spec_response(h, w, hs, wp):
    s = spec_wave(hs, wp, w)
    return h ** 2 * s
