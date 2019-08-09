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


class Candidate():
    def __init__(self, analyses_root, park_label, wtg_label):

        self.settings = SettingsClass(analyses_root, park_label, wtg_label)
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

        # mesh = Mesh(loads.pd.ppoints, loads.pd.ppanels[ipanel])
        # mesh.show()

        # Get section forces
        section_point = [3.5, 0, 0]
        section_normal = [1, 0, 0]
        self.imass, self.ipanel = self.tran_fun.get_section_index(section_point, section_normal)
        self.f_sec1 = self.tran_fun.assemble_forces(self.imass, self.ipanel, moment_ref_point=section_point)

        mesh = Mesh(self.loads.pd.ppoints, self.loads.pd.ppanels[self.ipanel])
        # mesh.show()

        self.sig_dyn_tot = self.panel_cc.axial_stress(self.f_sec1['Dynamic']['Total'])
        self.sig_stat_tot = self.panel_cc.axial_stress(self.f_sec1['Static']['Total'])

        p_dyn_lat = self.panel_cc.lateral_pressure(self.f_part1['Dynamic']['Total'])
        p_stat_lat = self.panel_cc.lateral_pressure(self.f_part1['Static']['Total'])

        self.hs = 12
        self.tp = 14
        # gamma = env.gamma(hs, tp)
        # sig_r = np.abs(sig_dyn_tot ** 2) * env.s_jonswap(hs=hs, wp=2 * np.pi / tp, w=loads.w, gamma=gamma)[:, np.newaxis]
        # dw = loads.w[1] - loads.w[0]
        # sig_r_m0 = sum(sig_r) * dw
        # tz = env.tp2tz(tp, gamma)
        # nz = 3 * 3600 / tz
        # sig_r_max = np.sqrt(sig_r_m0)*(np.sqrt(2*np.log(nz))+0.5772/np.sqrt(2*np.log(nz)))

        self.stwc1 = ec.Short_Term_Wave_Conditions(hs=self.hs, tp=self.tp)

        print('\nExpected largest maximum dynamic normal stress for Hs = {:4.1f} m and Tp = {:4.1f} s'.format(self.hs,
                                                                                                              self.tp))
        for i in range(self.loads._nbeta):
            print('\tWavedir {:5.1f} deg: {:6.1f} MPa'.format(self.loads._beta[i] * 180 / np.pi,
                                                              self.stwc1.expected_largest_maximum(
                                                                      self.sig_dyn_tot[:, i],
                                                                      self.loads.w) / 10 ** 6))
        print('\nStatic stress: {:1.1f} MPa'.format(self.sig_stat_tot / 10 ** 6))

        print(
                '\nExpected largest maximum dynamic lateral pressure for Hs = {:4.1f} m and Tp = {:4.1f} s'.format(
                    self.hs,
                    self.tp))
        for i in range(self.loads._nbeta):
            print('\tWavedir {:5.1f} deg: {:6.1f} kPa'.format(self.loads._beta[i] * 180 / np.pi,
                                                              self.stwc1.expected_largest_maximum(p_dyn_lat[:, i],
                                                                                                  self.loads.w) / 10 ** 3))
        print('\nStatic lateral force: {:1.1f} kPa'.format(p_stat_lat / 10 ** 3))

        print('\n\nUtilizations')
        gamma_m = 1.15
        load_factor = 1.3
        sigma_y = 235000000 / 1.15
        bc = 'pinned'
        dpu = self.panel_cc.dynamic_panel_utilization(sigma_y, bc, self.f_sec1['Dynamic']['Total'],
                                                      self.f_part1['Dynamic']['Total'],
                                                      self.hs, self.tp, freq=self.loads.w, pos_y_side=True)
        for i in range(self.loads._nbeta):
            print('\tWavedir {:5.1f} deg: {:6.2f}'.format(self.loads._beta[i] * 180 / np.pi, dpu[i] * load_factor))

        dpu = self.panel_cc.dynamic_panel_utilization(sigma_y, bc, self.f_sec1['Dynamic']['Total'],
                                                      self.f_part1['Dynamic']['Total'],
                                                      self.hs, self.tp, freq=self.loads.w, pos_y_side=False)
        for i in range(self.loads._nbeta):
            print('\tWavedir {:5.1f} deg: {:6.2f}'.format(self.loads._beta[i] * 180 / np.pi, dpu[i] * load_factor))

        ifreq = 0
        idir = 0
        irad = 4
        # nprob = len(NEMOH_DIR) + sum(NEMOH_DOF)
        # iprob = 2
        # nfreq  = len(w)
        # problem = (ifreq - 1) * nprob + iprob

        # print('Problem: {}'.format(problem))
        NEMOH_DOF = [1, 1, 1, 1, 1, 1]

        if False:
            # hdp.show_pressure(ifreq, idir, pressure_type='Froude-Krylof',axis=2)
            # hdp.show_pressure(ifreq, idir, pressure_type='Diffraction', axis=2)
            tran_fun.show_pressure(ifreq, irad, pressure_type='Radiation', axis=0)
            pass

        # print(hdp.ma[0,:,0,0])
        # print(hdp._fe_amp[0,:,dof])

        idof = np.array([i for i, x in enumerate(NEMOH_DOF) if x])

        if False:
            print('\nEigenvalue sollution WITH added mass')
            tb.eigenvalprint(loads.m + loads.ma[ifreq, :, :], loads.k)

            print('\n')

            print('\nRadiation damping:')
            tb.matprint(loads.c_hyd[ifreq])
            print('\nWater plane stiffness:')
            tb.matprint(loads.k)
            print('\nStatic mass [tonne]:')
            tb.matprint(loads.m / 1000)
            print('\nAdded mass [tonne]:')
            tb.matprint(loads.ma[ifreq] / 1000)
            print('\nExcitation force:')
            tb.matprint(np.abs(loads.fe[ifreq, idir, :]))
            f_fk = loads.p2f('Froude-Krylof', ifreq, idir)
            f_diff = loads.p2f('Diffraction', ifreq, idir)
            f_exc = f_fk + f_diff

            # print('\n')
            # tb.matprint(np.abs(nemoh.get_section_values(f_exc, loads.pd.ppanel_centers, [0, 0, 0], [1, 0, 0])))
        # np.set_printoptions(precision=3)

        idir = 0

        # --------------------------------------------------------------------------------------------------------------
        # PLOT RESULTS
        # --------------------------------------------------------------------------------------------------------------
        if False:
            h = tran_fun.get_rao(idir)
            fig, axs = plt.subplots(3, 2)
            w = 2 * np.pi / loads.w
            axs[0, 0].plot(w, abs(h[:, 2]), 'tab:orange')
            axs[0, 0].set_title('Rao heave')
            axs[0, 1].plot(w, abs(h[:, 4]), 'tab:orange')
            axs[0, 1].set_title('Rao pitch')
            axs[1, 0].plot(w, loads.ma[:, 2, 2] + loads.m[2, 2], 'tab:green')
            axs[1, 0].set_title('m+ma heave')
            axs[1, 1].plot(w, loads.ma[:, 4, 4] + loads.m[4, 4], 'tab:green')
            axs[1, 1].set_title('m+ma pitch')
            axs[2, 0].plot(w, abs(loads._fe[:, idir, 2]), 'tab:blue')
            axs[2, 0].set_title('fe heave')
            axs[2, 1].plot(w, abs(loads._fe[:, idir, 4]), 'tab:blue')
            axs[2, 1].set_title('fe pitch')
            plt.show()

    def print_report(self):
        self.settings._report._doc.generate_pdf(clean_tex=False)
