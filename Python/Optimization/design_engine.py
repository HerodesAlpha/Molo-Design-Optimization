__author__ = "Eivind Sonju"
__copyright__ = "Copyright (C) 2017-2019 Verbun AS. All rights reserved."
__version__ = "2.0"

import logging
import multiprocessing
import pickle
import matplotlib.pyplot as plt
import numpy as np
from logutils.queue import QueueListener
import tool_box as tb
import calculations
import model_set_up as msu
import stability
from Nemoh_Frontend import nemoh_frontend as nf
from common import SettingsClass
from pyNemoh.structure import BaseStructure
from loads import Sea_and_Inertia_Loads
from meshmagick.mesh import Mesh
import code_check as cc
import environmental_conditions as ec
import json
import os
from pathlib import Path


class Parameter_Space():
    def __init__(self):
        self._json_list = ['park', 'rna', 'tower', 'floater', 'analysis', 'design_basis']
        self._job_data = dict()

        # Collect template data
        for item in self._json_list:
            with open(Path(os.getcwd()).joinpath('templates').joinpath('{}_template.json'.format(item)), 'r') as f:
                self._job_data[item] = json.loads(f.read())

        # Save updated template to template dir
        for item in self._json_list:
            with open(Path(os.getcwd()).joinpath('templates').joinpath('{}_template.json'.format(item)), 'w') as f:
                f.write(json.dumps(self._job_data[item], indent=4, sort_keys=True))

    @property
    def nradial(self):
        return self._job_data['floater']['Number of radials']

    @nradial.setter
    def nradial(self, val):
        self._job_data['floater']['Number of radials'] = val

    @property
    def column_diameter(self):
        return self._job_data['floater']['Radial']['Column']['Diameter']

    @column_diameter.setter
    def column_diameter(self, val):
        self._job_data['floater']['Radial']['Column']['Diameter'] = val

    @property
    def ncol(self):
        return self._job_data['floater']['Number of columns']

    @ncol.setter
    def ncol(self, val):
        self._job_data['floater']['Number of columns'] = val

    @property
    def gap(self):
        return self._job_data['floater']['Gap factor']

    @gap.setter
    def gap(self, val):
        self._job_data['floater']['Gap factor'] = val

    @property
    def height(self):
        return self._job_data['floater']['Gap factor']

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
    def stiffener_height(self,val):
        self._job_data['floater']['Radial']['Flange']['Lower']['Stiffener']['Longitudinal']['Height'] = val



class Candidate():
    def __init__(self, analyses_root, park_label, wtg_label, parameter_space, case_label_type):

        self.settings = SettingsClass(analyses_root, park_label, wtg_label, parameter_space, case_label_type)
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

    def init_model(self, state=None):
        if state == 'New':
            print('\n--------------------------------------------------------------------------------------------')
            print('CREATE MODEL')
            print('--------------------------------------------------------------------------------------------')

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

    def intact_stability_ratio(self):
        print('\n--------------------------------------------------------------------------------------------')
        print('STABILITY ANALYSIS')
        print('--------------------------------------------------------------------------------------------')
        return stability.intact_stability(self.settings, self.hs_floater)

    def hydrodynamic_analysis(self):
        print('\n--------------------------------------------------------------------------------------------')
        print('NEMOH ANALYSIS')
        print('--------------------------------------------------------------------------------------------')
        self.settings.remove_old_db()
        self.queue = multiprocessing.Queue(-1)
        self.ql = QueueListener(self.queue, *logging.getLogger().handlers)
        self.ql.start()
        nf.run(self.settings._job_data['analysis'], self.queue)
        self.ql.stop()

    def structural_analysis(self):
        self.loads = Sea_and_Inertia_Loads(self.settings)
        self.tran_fun = calculations.TransferFunctions(self.settings, self.loads)
        self.panel_cc = cc.Panel(self.settings)

        self.imass, self.ipanel = self.tran_fun.get_flange_panel_index()
        self.f_part1 = self.tran_fun.assemble_forces(self.imass, self.ipanel, moment_ref_point=[0, 0, 0])

        # Consider first radial

        # Get section forces
        sp_x = self.settings.floater_data['Central column diameter'] / 2
        sp_z = self.settings.floater_data['Radial']['Heigth'] / 2 - self.hs_floater.hs_data['draught']
        section_point = [sp_x, 0, sp_z]  # Used for moment reference
        section_normal = [1, 0, 0]
        self.imass, self.ipanel = self.tran_fun.get_section_index(section_point, section_normal)
        self.f_sec1 = self.tran_fun.assemble_forces(self.imass, self.ipanel, moment_ref_point=section_point)

        yr = self.settings.park_data['Design Basis']['ULS']['Return period']
        area = self.settings.park_data['Design Basis']['Area']
        self.ltwc1 = ec.Long_Term_Wave_Conditions(area=area)
        self.cl = self.ltwc1.contour_line(yr)

        # Create list of short terms from contour line
        self.contourline = []
        for hs, tz in self.cl:
            self.contourline.append(ec.Short_Term_Wave_Conditions(hs=hs, tz=tz))

        dpu = np.zeros([self.loads.nbeta, 2],dtype=float)

        gamma_m = 1.15
        load_factor = 1.3
        sigma_y = 235000000 /gamma_m
        # Check lower, inner panel
        bc = 'pinned'

        for ibeta in range(self.loads.nbeta):
            f_sec = self.f_sec1['Dynamic']['Total'][:, ibeta, :]
            f_part = self.f_part1['Dynamic']['Total'][:, ibeta, :]
            sigma_x_p = self.panel_cc.axial_stress(f_sec, pos_y_side=True)
            sigma_x_n = self.panel_cc.axial_stress(f_sec, pos_y_side=False)
            p_lat = self.panel_cc.lateral_pressure(f_part)
            dpu[ibeta, 0] = self.panel_cc.minimize_panel_setion(sigma_y, bc, sigma_x_p, p_lat, self.contourline,
                                                                        freq=self.loads.w)
            dpu[ibeta, 1] = self.panel_cc.minimize_panel_setion(sigma_y, bc, sigma_x_n, p_lat, self.contourline,
                                                                        freq=self.loads.w)

        return {'Max UR':dpu.max(), 'panel_cc': self.panel_cc}

        # plt.plot(self.cl[:, 1], self.cl[:, 0], 'tab:orange')
        # plt.title('{} yr contourlines'.format(yr))
        # plt.ylabel('Hs')
        # plt.xlabel('Tz')
        # plt.show()

    def print_report(self):
        self.settings._report._doc.generate_pdf(clean_tex=False)
