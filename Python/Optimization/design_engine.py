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

class Parameter_Space():
    def __init__(self):
        pass

class Candidate():
    def __init__(self, analyses_root, park_label, wtg_label, case_label_type):

        self.settings = SettingsClass(analyses_root, park_label, wtg_label, case_label_type)
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
            pickle.dump(self.unit_model, open(self.settings.fio.data_io_dir.joinpath('unit_model.pkl'), 'wb'))
            pickle.dump(self.hs_floater, open(self.settings.fio.data_io_dir.joinpath('hs_floater.pkl'), 'wb'))

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
            self.unit_model = pickle.load(open(self.settings.fio.data_io_dir.joinpath('unit_model.pkl'), 'rb'))
            self.hs_floater = pickle.load(open(self.settings.fio.data_io_dir.joinpath('hs_floater.pkl'), 'rb'))
            for item in self.settings._json_list:
                with open(self.settings.fio.data_io_dir.joinpath('{}.json'.format(item)), 'r') as f:
                    self.settings._job_data[item] = json.loads(f.read())
        else:
            print('init_model state is either New or Old')
            exit()

    def run_intact_stability(self):
        print('\n--------------------------------------------------------------------------------------------')
        print('STABILITY ANALYSIS')
        print('--------------------------------------------------------------------------------------------')
        return stability.intact_stability(self.settings, self.hs_floater)

    def run_hydrodynamic_analysis(self):
        print('\n--------------------------------------------------------------------------------------------')
        print('NEMOH ANALYSIS')
        print('--------------------------------------------------------------------------------------------')
        self.settings.remove_old_db()
        self.queue = multiprocessing.Queue(-1)
        self.ql = QueueListener(self.queue, *logging.getLogger().handlers)
        self.ql.start()
        nf.run(self.settings._job_data['analysis'], self.queue)
        self.ql.stop()

    def run_postprocessing(self):
        self.loads = Sea_and_Inertia_Loads(self.settings)
        self.tran_fun = calculations.TransferFunctions(self.settings, self.loads)
        self.panel_cc = cc.Panel(self.settings)

        self.imass, self.ipanel = self.tran_fun.get_flange_panel_index()
        self.f_part1 = self.tran_fun.assemble_forces(self.imass, self.ipanel, moment_ref_point=[0, 0, 0])

        # Consider first radial

        # Get section forces
        sp_x = self.settings.floater_data['Central column diameter']/2
        sp_z = self.settings.floater_data['Radial height']/2 - self.hs_floater.hs_data['draught']
        section_point = [sp_x, 0, sp_z] # Used for moment reference
        section_normal = [1, 0, 0]
        self.imass, self.ipanel = self.tran_fun.get_section_index(section_point, section_normal)
        self.f_sec1 = self.tran_fun.assemble_forces(self.imass, self.ipanel, moment_ref_point=section_point)

        yr = 50
        self.ltwc1 = ec.Long_Term_Wave_Conditions(area=4)
        self.cl = self.ltwc1.contour_line(yr)

        # Create list of short terms from contour line
        self.stwc1_list = []
        for hs, tz in self.cl:
            self.stwc1_list.append(ec.Short_Term_Wave_Conditions(hs=hs, tz=tz))

        for stwcl in self.stwc1_list:
            print('\n\nUtilizations for Hs = {:5.2f} and Tz = {:5.2f}'.format(stwcl.hs, stwcl.tz))
            gamma_m = 1.15
            load_factor = 1.3
            sigma_y = 235000000 / 1.15
            bc = 'pinned'
            dpu = self.panel_cc.dynamic_panel_utilization(sigma_y, bc, self.f_sec1['Dynamic']['Total'],
                                                          self.f_part1['Dynamic']['Total'],
                                                          stwcl, freq=self.loads.w, pos_y_side=True)
            for i in range(self.loads._nbeta):
                print('\tWavedir {:5.1f} deg: {:6.2f}'.format(self.loads._beta[i] * 180 / np.pi, dpu[i] * load_factor))

            dpu = self.panel_cc.dynamic_panel_utilization(sigma_y, bc, self.f_sec1['Dynamic']['Total'],
                                                          self.f_part1['Dynamic']['Total'],
                                                          stwcl, freq=self.loads.w, pos_y_side=False)
            for i in range(self.loads._nbeta):
                print('\tWavedir {:5.1f} deg: {:6.2f}'.format(self.loads._beta[i] * 180 / np.pi, dpu[i] * load_factor))

        plt.plot(self.cl[:, 1], self.cl[:, 0], 'tab:orange')
        plt.title('{} yr contourlines'.format(yr))
        plt.ylabel('Hs')
        plt.xlabel('Tz')
        plt.show()

    def print_report(self):
        self.settings._report._doc.generate_pdf(clean_tex=False)
