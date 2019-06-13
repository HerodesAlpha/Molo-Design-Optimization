__author__ = "Eivind Sonju"
__copyright__ = "Copyright (C) 2017-2019 Verbun AS. All rights reserved."
__version__ = "2.0"

import h5py
from pylab import *
import pickle
import force
import nemoh
import tool_box as tb
from common import PhysicalQuantities
from meshmagick.mesh import Mesh
from pyNemoh.structure import BaseStructure


class HydroCoefficients(PhysicalQuantities, object):
    # TODO: Get added mass at zero and infinite frequency
    def __init__(self, settings):
        super().__init__()
        h5_bs = BaseStructure()
        with h5py.File(settings.fio.nemoh_root.joinpath('db.hdf5'), "r") as hdf5_db:

            #hdf5_db['results']['fk_pressure_raw'][0]

            self._w = hdf5_db[h5_bs.H5_RESULTS_CASE_W][:]
            self._nw = len(self._w)
            self._dof = [1,1,1,1,1,1]
            self._ndof = int(np.sum(self._dof))
            self._beta = hdf5_db[h5_bs.H5_RESULTS_CASE_BETA][:]
            self._nbeta = len(self._beta)
            #self._sym = sym

            step = 0
            sys.stdout.write("\nInit hydro:\n")
            # -------------------------------------------------------

            self._nemoh_mesh = pickle.load(open(settings.fio.data_io_dir.joinpath('nemoh_mesh.pkl'), 'rb'))
            self._pd = tb.PanelData(self._nemoh_mesh.vertices, self._nemoh_mesh.faces) # TODO: Get vertices and points

            # -------------------------------------------------------
            step += 1
            sys.stdout.write("\t({})Get forces\n".format(step))
            self._fe = hdf5_db[h5_bs.H5_RESULTS_EXCITATION_FORCES][:]
            # -------------------------------------------------------
            step += 1
            sys.stdout.write("\t({})Get static mass and waterplane stiffness\n".format(step))
            self._m, self._k = tb.load_M_and_K(settings.fio.data_io_dir)
            # -------------------------------------------------------
            step += 1
            sys.stdout.write("\t({})Get added mass and damping\n".format(step))
            self._ma = hdf5_db[h5_bs.H5_RESULTS_ADDED_MASS][:]
            self._c_hyd = hdf5_db[h5_bs.H5_RESULTS_RADIATION_DAMPING][:]
            #self._ma, self._c_hyd = nemoh.get_ab(fio, dof, w, dir)
            # -------------------------------------------------------
            # nproblems = (len(dir) + sum(dof)) * len(w)
            # with open(fio.nemoh_root.joinpath('Normalvelocities.dat')) as f:
            #     if not int(f.readline()) == nproblems:
            #         print('nproblems mismatch')
            #         exit()
            # -------------------------------------------------------
            step += 1
            sys.stdout.write("\t({})Get pressures\n".format(step))
            pressure_file = settings.fio.data_io_dir.joinpath('hydro_pressures.pkl')
            if pressure_file.is_file() and 0:
                self._p = pickle.load(open(pressure_file, 'rb'))
                sys.stdout.write("\t\tUnPickled from {}\n".format(pressure_file))
                for key in self._p.keys():
                    print("\t\t\t{}".format(key))

            else:
                self._p = dict()
                sys.stdout.write("\t\tHydro static\n")
                self._p['Hydro static'] = (self._rho_sw * self._grav) * self._pd.ppanel_centers[:, 2]
                sys.stdout.write("\t\tFroude-Krylof \n")
                self._p['Froude-Krylof'] = hdf5_db[h5_bs.H5_RESULTS_FK_PRESSURE_RAW][:]


                pressure = hdf5_db[h5_bs.H5_RESULTS_PRESSURE][:]


                sys.stdout.write("\t\tDiffraction\n")
                self._p['Diffraction'] = np.zeros([self._nw, self._nbeta, self._pd.npanels], dtype=complex)
                for iw in range(self._nw):
                    for ibeta in range(self._nbeta):
                        pn = nemoh.diffraction_problem_number(iw, ibeta, self._nbeta, self._ndof)
                        # print('\t\t\tProblem {}'.format(pn))
                        self._p['Diffraction'][iw, ibeta, :] = pressure[pn-1,:]

                sys.stdout.write("\t\tRadiation\n")
                self._p['Radiation'] = np.zeros([self._nw, self._ndof, self._pd.npanels], dtype=complex)
                for iw in range(self._nw):
                    for iradiation in range(self._ndof):
                        pn = nemoh.radiation_problem_number(iw, iradiation, self._nbeta, self._ndof)
                        # print('\t\t\tProblem {}'.format(pn))
                        self._p['Radiation'][iw, iradiation, :] = pressure[pn-1,:]

                pickle.dump(self._p, open(pressure_file, "wb"))

            print('\n{} initialized\n'.format(self.__str__()))

    def show_pressure(self, ifreq, pressure_index, axis, pressure_type):
        # Pressure index is either force degree of freedom or wave direction
        # Axis is the requested decomposed direction of the pressure
        xyz = np.zeros([self._pd._npanel, 3])
        xyz[:, axis] = 1
        nemoh_mesh = Mesh(self._pd.ppoints, self._pd.ppanels)
        p=self._p[pressure_type][ifreq, pressure_index,:]
        n=self._pd.ppanel_normals
        vec = (n * np.real(p)[:, np.newaxis]) * xyz
        h = force.show_force(nemoh_mesh, self._pd.ppanel_centers, vec)
        h.show()

    def p2f(self, ifreq, pressure_index, pressure_type):
        # Pressure to force
        f_normal = np.zeros((self._pd.npanels), dtype=np.complex)
        f = np.zeros((self._pd.npanels, 3), dtype=np.complex)
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
        return [e for i, e in enumerate(self._beta)]

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
