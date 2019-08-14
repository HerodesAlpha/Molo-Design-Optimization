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
import warnings

class Sea_and_Inertia_Loads(PhysicalQuantities, object):
    # TODO: Get added mass at zero and infinite frequency
    # TODO: Correct radiation pressure due to artificial distance between upper and lower face of lower flange.
    def __init__(self, settings):
        super().__init__()
        h5_bs = BaseStructure()
        self._settings = settings
        with h5py.File(settings.fio.nemoh_root.joinpath('db.hdf5'), "r") as hdf5_db:

            # hdf5_db['results']['fk_pressure_raw'][0]

            self._w = hdf5_db[h5_bs.H5_RESULTS_CASE_W][:]
            self._nw = len(self._w)
            self._dof = [1, 1, 1, 1, 1, 1]
            self._ndof = int(np.sum(self._dof))
            self._beta = hdf5_db[h5_bs.H5_RESULTS_CASE_BETA][:]
            self._nbeta = len(self._beta)
            # self._sym = sym

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
                b = a.copy()
                b[:, 1] = -b[:, 1]  # Create new set of nodes mirrored about xz
                self._nemoh_mesh_vertices = np.vstack((a, b))  # Append the new set to the old set
                a = self._nemoh_mesh_faces
                b = a.copy()
                b = b + int(self._nemoh_mesh_vertices.shape[0] / 2)  # Create a new set of faces from the old
                # set and renumber by adding int(nvertices)
                b[:, :] = b[:, ::-1]  # Flip normals on mirrored faces (reverse nodes)
                self._nemoh_mesh_faces = np.vstack((a, b))
                del a
                del b

                # print('Bottom {}'.format(np.min(self._nemoh_mesh_vertices[:,2])))
                if True: # Move lower faces to correct position
                    ind = self._nemoh_mesh_vertices[:, 2] <= np.min(self._nemoh_mesh_vertices[:, 2]) * 0.99
                    self._nemoh_mesh_vertices[ind, 2] += settings.thin_panel_offset - settings.flange_thickness
                    # print('Bottom {}'.format(np.min(self._nemoh_mesh_vertices[:,2])))
                    del ind
                    lower_face_corrected_z_pos = True
                else:
                    lower_face_corrected_z_pos = False

            self._pd = tb.PanelData(self._nemoh_mesh_vertices, self._nemoh_mesh_faces)  # TODO: Get vertices and points
            self._an = self.pd.ppanel_areas[:, np.newaxis] * self.pd.ppanel_normals
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
            self._ma = hdf5_db[h5_bs.H5_RESULTS_ADDED_MASS][:] #TODO: Calc added mass inf
            self._c_hyd = hdf5_db[h5_bs.H5_RESULTS_RADIATION_DAMPING][:]
            # -------------------------------------------------------

            step += 1
            sys.stdout.write("\t({})Get pressures and forces\n".format(step))
            pressure_file = settings.fio.data_io_dir.joinpath('hydro_pressures.pkl')
            force_file = settings.fio.data_io_dir.joinpath('hydro_forces.pkl')

            load_files = False
            save_files = False

            if pressure_file.is_file() and force_file.is_file() and load_files:
                self._pressure = pickle.load(open(pressure_file, 'rb'))
                self._force = pickle.load(open(force_file, 'rb'))
                sys.stdout.write("\t\tUnPickled from {}\n".format(pressure_file))
                for key in self._pressure.keys():
                    print("\t\t\t{}".format(key))
                sys.stdout.write("\t\tUnPickled from {}\n".format(force_file))
                for key in self._force.keys():
                    print("\t\t\t{}".format(key))

            else:

                self._pressure = dict()
                self._force = dict()

                sys.stdout.write("\t\tStatic buoyancy pressure\n")
                # TODO: z coordinate of lower face of flange is artificially low to avoid num. instab. Dont use for hydro stat. pressure
                assert lower_face_corrected_z_pos
                self._pressure['Buoyancy'] = (self._rho_sw * self._gravity) * self._pd.ppanel_centers[:, 2]

                sys.stdout.write("\t\tStatic buoyancy force\n")
                self._force['Buoyancy'] = -self._pressure['Buoyancy'][:,np.newaxis]*self._an


                sys.stdout.write("\t\tFroude-Krylof pressure\n")
                # TODO: z coordinate of lower face of flange is artificially low to avoid num. instab. Dont use for FK
                self._pressure['Froude-Krylof'] = hdf5_db[h5_bs.H5_RESULTS_FK_PRESSURE_RAW][:]

                sys.stdout.write("\t\tFroude-Krylof force\n")
                self._force['Froude-Krylof'] = -self._pressure['Froude-Krylof'][:,:,:,np.newaxis]*self._an[np.newaxis,np.newaxis,:,:]

                nemoh_pressure = hdf5_db[h5_bs.H5_RESULTS_PRESSURE][:]

                sys.stdout.write("\t\tDiffraction pressure\n")
                self._pressure['Diffraction'] = np.zeros([self._nw, self._nbeta, self._pd.npanels], dtype=complex)

                sys.stdout.write("\t\tDiffraction force\n")
                for iw in range(self._nw):
                    for ibeta in range(self._nbeta):
                        pn = nemoh.diffraction_problem_number(iw, ibeta, self._nbeta, self._ndof)
                        # print('\t\t\tProblem {}'.format(pn))
                        self._pressure['Diffraction'][iw, ibeta, :] = nemoh_pressure[pn - 1, :]
                self._force['Diffraction'] = -self._pressure['Diffraction'][:,:,:,np.newaxis]*self._an[np.newaxis,np.newaxis,:,:]

                sys.stdout.write("\t\tRadiation pressure\n")
                self._pressure['Radiation'] = np.zeros([self._nw, self._ndof, self._pd.npanels], dtype=complex)

                sys.stdout.write("\t\tRadiation force\n")
                for iw in range(self._nw):
                    for iradiation in range(self._ndof):
                        pn = nemoh.radiation_problem_number(iw, iradiation, self._nbeta, self._ndof)
                        # print('\t\t\tProblem {}'.format(pn))
                        self._pressure['Radiation'][iw, iradiation, :] = nemoh_pressure[pn - 1, :]
                self._force['Radiation'] = -self._pressure['Radiation'][:,:,:,np.newaxis]*self._an[np.newaxis,np.newaxis,:,:]

            if save_files:
                with open(pressure_file, "wb") as f:
                    pickle.dump(self._pressure, f)
                with open(force_file, "wb") as f:
                    pickle.dump(self._force, f)

            print('\n{} initialized\n'.format(self.__str__()))

    def show_pressure(self, ifreq, pressure_index, axis, pressure_type):
        # Pressure index is either force degree of freedom or wave direction
        # Axis is the requested decomposed direction of the pressure
        xyz = np.zeros([self._pd._npanel, 3])
        xyz[:, axis] = 1
        nemoh_mesh = Mesh(self._pd.ppoints, self._pd.ppanels)
        p = self._pressure[pressure_type][ifreq, pressure_index, :]
        n = self._pd.ppanel_normals
        vec = (n * np.imag(p)[:, np.newaxis]) * xyz
        h = force.show_force(nemoh_mesh, self._pd.ppanel_centers, vec)
        h.show()

    def p2f(self, p, ifreq=None, idir=None):
        if isinstance(p, str):
            if p == 'Hydro_static':
                p_cmplx = self._pressure[p]
            elif p in self._pressure.keys():
                p_cmplx = self._pressure[p][ifreq, idir, :]
        else:
            p_cmplx = p
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

    @property
    def nw(self):
        return self._nw


    @property
    def nbeta(self):
        return self._nbeta


