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


class Sea_and_Inertia_Loads(PhysicalQuantities, object):
    # TODO: Get added mass at zero and infinite frequency
    # TODO: Correct radiotion pressure due to artificial distance between upper and lower face of lower flange.
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

            with open(settings.fio.data_io_dir.joinpath('nemoh_mesh_vertices.pkl'), 'rb') as f:
                self._nemoh_mesh_vertices = pickle.load(f)
            with open(settings.fio.data_io_dir.joinpath('nemoh_mesh_faces.pkl'), 'rb') as f:
                self._nemoh_mesh_faces = pickle.load(f)

            # Reshape vertices and faces to full model
            if settings.use_symmmetri == True:
                a = self._nemoh_mesh_vertices
                b = a.copy(); b[:,1] = -b[:,1] # Create new set of nodes mirrored about xz
                self._nemoh_mesh_vertices = np.vstack((a,b)) # Append the new set to the old set
                a = self._nemoh_mesh_faces
                b = a.copy(); b = b + int(self._nemoh_mesh_vertices.shape[0]/2) # Create a new set of faces from the old
                # set and renumber by adding int(nvertices)
                b[:,:] = b[:,::-1] # Flip normals on mirrored faces (reverse nodes)
                self._nemoh_mesh_faces = np.vstack((a, b))
                del a; del b

                #print('Bottom {}'.format(np.min(self._nemoh_mesh_vertices[:,2])))
                if False:
                    ind = self._nemoh_mesh_vertices[:,2] <= np.min(self._nemoh_mesh_vertices[:,2])*0.99
                    self._nemoh_mesh_vertices[ind, 2] += settings.thin_panel_offset - settings.flange_thickness
                    #print('Bottom {}'.format(np.min(self._nemoh_mesh_vertices[:,2])))
                    del ind




            self._pd = tb.PanelData(self._nemoh_mesh_vertices, self._nemoh_mesh_faces) # TODO: Get vertices and points

            # -------------------------------------------------------
            step += 1
            sys.stdout.write("\t({})Get forces\n".format(step))
            self._fe = hdf5_db[h5_bs.H5_RESULTS_EXCITATION_FORCES][:]
            # self._fe = np.zeros([self._nw,self._nbeta,6,6],dtype=complex)
            # I = np.identity(6)
            # fe_diag=hdf5_db[h5_bs.H5_RESULTS_EXCITATION_FORCES][:]
            # for i in range(self._nw):
            #     for j in range(self._nbeta):
            #         self._fe[i,j,:,:] = fe_diag[i,j,:]*I
            # -------------------------------------------------------
            step += 1
            sys.stdout.write("\t({})Get static mass and waterplane stiffness\n".format(step))
            self._m, self._k = tb.load_M_and_K(settings.fio.data_io_dir)
            # -------------------------------------------------------
            step += 1
            sys.stdout.write("\t({})Get added mass and damping\n".format(step))
            self._ma = hdf5_db[h5_bs.H5_RESULTS_ADDED_MASS][:]
            # self._ma_inf = hdf5_db[h5_bs.H5_RESULTS_ADDED_MASS_INFINITE][:] TODO: Calc added mass inf
            # self._ma_zero = hdf5_db[h5_bs.H5_RESULTS_ADDED_MASS_ZERO][:]
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
                self._p['Hydro_static'] = (self._rho_sw * self._grav) * self._pd.ppanel_centers[:, 2] # TODO: z coordinate of lower face of flange is artificially low to avoid num. instab.. Dont use for hydro stat. pressure
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

                with open(pressure_file, "wb") as f:
                    pickle.dump(self._p,f )

            print('\n{} initialized\n'.format(self.__str__()))

    def show_pressure(self, ifreq, pressure_index, axis, pressure_type):
        # Pressure index is either force degree of freedom or wave direction
        # Axis is the requested decomposed direction of the pressure
        xyz = np.zeros([self._pd._npanel, 3])
        xyz[:, axis] = 1
        nemoh_mesh = Mesh(self._pd.ppoints, self._pd.ppanels)
        p=self._p[pressure_type][ifreq, pressure_index,:]
        n=self._pd.ppanel_normals
        vec = (n * np.imag(p)[:, np.newaxis]) * xyz
        h = force.show_force(nemoh_mesh, self._pd.ppanel_centers, vec)
        h.show()



    def p2f(self,pressure_type, ifreq=None, idir = None):

        if pressure_type=='Hydro_static':
            p_cmplx=self._p[pressure_type]
        else:
            p_cmplx=self._p[pressure_type][ifreq, idir, :]

        npanels = self.pd.ppanels.shape[0]
        f_normal = np.zeros((npanels), dtype=np.complex)
        f = np.zeros((npanels, 3), dtype=np.complex)
        for i, panel in enumerate(self.pd.ppanels):
            f_normal[i] = p_cmplx[i] * self.pd.ppanel_areas[i]
            for j in range(3):
                f[i, j] = -f_normal[i] * self.pd.ppanel_normals[i, j]
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

    # @property
    # def ma_inf(self):
    #     return self._ma_inf

    # @property
    # def ma_zero(self):
    #     return self._ma_zero

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
    def c_hyd(self):
        return self._c_hyd

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
