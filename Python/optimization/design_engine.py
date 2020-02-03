__author__ = "Eivind Sonju"
__copyright__ = "Copyright (C) 2017-2019 Verbun AS. All rights reserved."
__version__ = "2.0"

import logging
import multiprocessing
import pickle
import matplotlib.pyplot as plt
import numpy as np
from logutils.queue import QueueListener
import core.tool_box as tb

import core.model_set_up as msu
from core import stability
from nemoh_frontend import nemoh_frontend as nf
from core.common import SettingsClass
from pyNemoh.structure import BaseStructure

from core.meshmagick.mesh import Mesh
import core.code_check as cc
import core.environmental_conditions as ec
import json
import os
from pathlib import Path
import h5py
from pyNemoh.structure import BaseStructure
from core.response import ResponseModel
from core.loads import Sea_and_Inertia_Loads
from multiprocessing import freeze_support
import time

class Parameter_Space():
    def __init__(self, templates_dir):
        self._json_list = ['park', 'wtg', 'floater', 'analysis', 'design_basis']
        self._job_data = dict()

        self._templates_dir = templates_dir

        # Collect template data
        for item in self._json_list:
            with open(self._templates_dir.joinpath('{}_template.json'.format(item)), 'r') as f:
                self._job_data[item] = json.loads(f.read())

        # Save updated template to template dir
        for item in self._json_list:
            with open(self._templates_dir.joinpath('{}_template.json'.format(item)), 'w') as f:
                f.write(json.dumps(self._job_data[item], indent=4, sort_keys=True))

    @property
    def nradial(self):
        return self._job_data['floater']['Number of radials']

    @property
    def templates_dir(self):
        return self._templates_dir

    @nradial.setter
    def nradial(self, val):
        self._job_data['floater']['Number of radials'] = val

    @property
    def radial_column_diameter(self):
        return self._job_data['floater']['Radial']['Column']['Diameter']

    @radial_column_diameter.setter
    def radial_column_diameter(self, val):
        self._job_data['floater']['Radial']['Column']['Diameter'] = val

    @property
    def radial_column_thickness(self):
        return self._job_data['floater']['Radial']['Column']['Thickness']

    @radial_column_thickness.setter
    def radial_column_thickness(self, val):
        self._job_data['floater']['Radial']['Column']['Thickness'] = val


    @property
    def ncol(self):
        return self._job_data['floater']['Radial']['Number of columns']

    @ncol.setter
    def ncol(self, val):
        self._job_data['floater']['Radial']['Number of columns'] = val

    @property
    def gap(self):
        return self._job_data['floater']['Gap factor']

    @gap.setter
    def gap(self, val):
        self._job_data['floater']['Gap factor'] = val

    @property
    def height(self):
        return self._job_data['floater']['Radial']['Heigth']

    @height.setter
    def height(self, val):
        self._job_data['floater']['Radial']['Heigth'] = val

    @property
    def job_data(self):
        return self._job_data

    @property
    def stiffener_height(self):
        return self._job_data['floater']['Radial']['Flange']['Lower']['Stiffener']['Longitudinal']['Height']

    @stiffener_height.setter
    def stiffener_height(self, val):
        self._job_data['floater']['Radial']['Flange']['Lower']['Stiffener']['Longitudinal']['Height'] = val

    @property
    def filling_ratio(self):
        return self._job_data['floater']['Ballast filling ratio'][1]

    @filling_ratio.setter
    def filling_ratio(self, val):
        self._job_data['floater']['Ballast filling ratio'][1] = val

    @property
    def lower_plate_width(self):
        return self._job_data['floater']['Radial']['Flange']['Lower']['Plate']['Width']

    @lower_plate_width.setter
    def lower_plate_width(self, val):
        self._job_data['floater']['Radial']['Flange']['Lower']['Plate']['Width'] = val

    @property
    def lower_flange_overlength(self):
        return self._job_data['floater']['Radial']['Flange']['Lower']['Overlength']

    @lower_flange_overlength.setter
    def lower_flange_overlength(self, val):
        self._job_data['floater']['Radial']['Flange']['Lower']['Overlength'] = val

    @property
    def lower_flange_thickness(self):
        return self._job_data['floater']['Radial']['Flange']['Lower']['eq_thick']

    @lower_flange_thickness.setter
    def lower_flange_thickness(self, val):
        self._job_data['floater']['Radial']['Flange']['Lower']['eq_thick'] = val


    @property
    def upper_flange_thickness(self):
        return self._job_data['floater']['Radial']['Flange']['Upper']['eq_thick']

    @upper_flange_thickness.setter
    def upper_flange_thickness(self, val):
        self._job_data['floater']['Radial']['Flange']['Upper']['eq_thick'] = val




class Candidate():
    def __init__(self, parameter_space, fio):

        self.settings = SettingsClass(parameter_space, fio)
        self.h5_bs = BaseStructure()

        # TODO: Allow for none equidistant frequencies

        self.settings.set_file_structure_and_report()
        self.settings.simulation_dir = str(self.settings.fio.nemoh_root)
        self.settings.save_job_settings()

        self.settings.mesh_name = None
        self.settings.use_dipols = False
        self.settings.use_symmmetri = True
        self.settings.do_equilibriate = True
        self.settings.thin_panel_offset = 0.5
        self.response = None
        self.loads = None

    def init_model(self, state=None):
        if state == 'New':
            start_time = time.time()
            print('\n--------------------\nCREATE MODEL\n--------------------')

            self.unit_model, self.hs_floater = msu.init_models(self.settings)
            with open(self.settings.fio.data_io_dir.joinpath('unit_model.pkl'), 'wb') as f:
                pickle.dump(self.unit_model, f)
            with open(self.settings.fio.data_io_dir.joinpath('hs_floater.pkl'), 'wb') as f:
                pickle.dump(self.hs_floater, f)

            self.settings.thin_panels = []

            self.m, self.k = tb.load_M_and_K(self.settings.fio.data_io_dir)

            print('\nEigenvalue sollution WITHOUT added mass (given as lambda^0.5)')
            print('\nMass matrix:')
            tb.matprint(self.m)
            print('\nStiffness matrix:')
            tb.matprint(self.k)
            print('')
            tb.eigenvalprint(self.m, self.k)
            print("CREATE MODEL took {:1.2f} seconds ".format(time.time() - start_time))
        elif state == 'Old':
            with open(self.settings.fio.data_io_dir.joinpath('unit_model.pkl'), 'rb') as f:
                self.unit_model = pickle.load(f)
            with open(self.settings.fio.data_io_dir.joinpath('hs_floater.pkl'), 'rb') as f:
                self.hs_floater = pickle.load(f)
            for item in self.settings._job_data:
                with open(self.settings.fio.data_io_dir.joinpath('{}.json'.format(item)), 'r') as f:
                    self.settings._job_data[item] = json.loads(f.read())
        else:
            print('init_model state is either New or Old')
            exit()

    def init_load(self):
        self.loads = Sea_and_Inertia_Loads(self.settings)


    def init_response(self, short_term_wave_condition):

        self.response = ResponseModel(self, short_term_wave_condition)

    def intact_stability_ratio(self):
        start_time = time.time()
        print('\n--------------------\nSTABILITY ANALYSIS\n--------------------')
        res = stability.intact_stability(self.settings, self.hs_floater)

        print("STABILITY ANALYSIS took {:1.2f} seconds ".format(time.time() - start_time))
        return res

    def hydrodynamic_analysis(self):
        start_time = time.time()
        print('\n--------------------\nNEMOH ANALYSIS\n--------------------')
        self.remove_old_db()
        self.queue = multiprocessing.Queue(-1)
        self.ql = QueueListener(self.queue, *logging.getLogger().handlers)
        self.ql.start()
        nf.run(self.settings._job_data['analysis'], self.queue)
        self.ql.stop()
        print("NEMOH ANALYSIS took {:1.2f} seconds ".format(time.time() - start_time))

    def structural_analysis(self):

        self.panel_cc = cc.Panel(self.settings)

        imass, ipanel = self.response.get_flange_panel_index()
        self.f_part1 = self.response.assemble_forces(imass, ipanel, moment_ref_point=[0, 0, 0])
        del imass, ipanel

        # Consider first radial

        # Get section forces
        sp_x = self.settings.floater_data['Central column diameter'] * (0.5)  # + 0.8 + 1 + 0.8 + 1)
        sp_z = self.settings.floater_data['Radial']['Heigth'] / 2 - self.hs_floater.hs_data['draught']
        section_point = [sp_x, 0, sp_z]  # Used for moment reference
        section_normal = [1, 0, 0]

        imass, ipanel = self.response.get_section_index(section_point, section_normal)
        self.f_sec1 = self.response.assemble_forces(imass, ipanel, moment_ref_point=section_point)
        del imass, ipanel

        yr = self.settings.park_data['Design Basis']['ULS']['Return period']
        area = self.settings.park_data['Design Basis']['Area']
        self.ltwc1 = ec.Long_Term_Wave_Conditions(area=area)
        self.cl = self.ltwc1.contour_line(yr)

        # Create list of short terms from contour line
        self.contourline = []
        for hs, tz in self.cl:
            self.contourline.append(ec.Short_Term_Wave_Conditions(hs=hs, tz=tz))

        dpu = np.zeros([2], dtype=float)

        gamma_m = 1.15
        load_factor = 1.3
        sigma_y = 235000000 / gamma_m
        # Check lower, inner panel
        bc = 'pinned'

        f_sec = self.f_sec1['Dynamic']['SUM']
        f_part = self.f_part1['Dynamic']['SUM']
        sigma_x_p = self.panel_cc.axial_stress(f_sec, pos_y_side=True)
        sigma_x_n = self.panel_cc.axial_stress(f_sec, pos_y_side=False)
        p_lat = self.panel_cc.lateral_pressure(f_part)

        dpu[0] = self.panel_cc.minimize_panel_setion(sigma_y, bc, sigma_x_p, p_lat, self.contourline,
                                                     freq=self.loads.w, nwdir=self.loads.nbeta)
        dpu[1] = self.panel_cc.minimize_panel_setion(sigma_y, bc, sigma_x_n, p_lat, self.contourline,
                                                     freq=self.loads.w, nwdir=self.loads.nbeta)

        def f_max(f):
            def f_elm(f):
                return np.array([stwcl.expected_largest_maximum(f, self.loads.w) for stwcl in self.contourline])

            f_out = np.zeros(6)
            for idof in range(6):
                f_out[idof] = np.array([f_elm(f[:, ibeta, idof]) for ibeta in range(self.loads.nbeta)]).max()
            return f_out

        with h5py.File(self.settings.fio.structural_dir.joinpath('structural.hdf5'), "w") as hdf5_structural_db:
            hdf5_structural_db.create_dataset('max utilization ratio', data=dpu.max())
            # hdf5_structural_db.create_dataset('panel code check', data=self.panel_cc)
            hdf5_structural_db.create_dataset('section force', data=f_max(f_sec))
            hdf5_structural_db.create_dataset('panel force', data=f_max(f_part))
            # hdf5_structural_db.create_dataset('contour line', data=self.contourline)
            hdf5_structural_db.create_dataset('gamma m', data=gamma_m)

        return {
                'Max UR'       : dpu.max(),
                'panel_cc'     : self.panel_cc,
                'section force': f_max(f_sec),
                'panel force'  : f_max(f_part)
        }

        # plt.plot(self.cl[:, 1], self.cl[:, 0], 'tab:orange')
        # plt.title('{} yr contourlines'.format(yr))
        # plt.ylabel('Hs')
        # plt.xlabel('Tz')
        # plt.show()

    def print_report(self):
        self.settings._report._doc.generate_pdf(clean_tex=False)

    def has_complete_hydrodynamic_db(self):
        db_file = self.settings._fio.nemoh_root.joinpath('db.hdf5')
        if db_file.is_file():
            with h5py.File(db_file, "a") as hdf5_hydro_db:
                if self.h5_bs.H5_RESULTS_EXCITATION_FORCES in hdf5_hydro_db.keys():
                    return True
        else:
            return False

    def has_gz(self):
        gz_file = self.settings._fio.stability_dir.joinpath('gz.txt')
        if gz_file.is_file():
            return True
        else:
            return False

    def remove_old_db(self):
        db_file = self.settings._fio.nemoh_root.joinpath('db.hdf5')
        if db_file.is_file():
            db_file.unlink()
            print('\ndb.hdf5 deleted from {}\n'.format(str(self.settings.fio.nemoh_root)))

    def has_stability_db(self):
        f = self.settings.fio.stability_dir.joinpath('stability.hdf5')
        if f.exists():
            return True
        else:
            return False

    def is_stable(self):
        f = self.settings.fio.stability_dir.joinpath('stability.hdf5')
        key = 'intact_stability_area_ratio'
        if f.exists():
            with h5py.File(self.settings.fio.stability_dir.joinpath('stability.hdf5'), "a") as hdf5_stability_db:
                if key in hdf5_stability_db.keys():
                    if hdf5_stability_db.get('intact_stability_area_ratio')[()] > 1.4:

                        return True
                    else:
                        return False

        else:
            return False

    def has_model(self):
        unit_model = self.settings.fio.data_io_dir.joinpath('unit_model.pkl')
        hs_floater = self.settings.fio.data_io_dir.joinpath('hs_floater.pkl')
        nemoh_mesh = self.settings.fio.nemoh_dir.joinpath(
                'MOLO_{}c_nemoh.dat'.format(self.settings.job_data['floater']['Radial']['Number of columns']))
        if unit_model.exists():
            if hs_floater.exists():
                if  nemoh_mesh.exists():

                    return True
                else:
                    print('nemoh_mesh does not exist')
                    return False

            else:
                print('hs_floater does not exist')
                return False
        else:
            print('unit_model does not exist')
            return False

    def has_structural_db(self):
        f = self.settings.fio.stability_dir.joinpath('structural.hdf5')
        if f.exists():
            return True
        else:
            return False
