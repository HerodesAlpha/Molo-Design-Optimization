"""
Test bench script for complete analysis workflow.

This script performs a full analysis workflow including:
- Model creation
- Stability analysis
- NEMOH hydrodynamic analysis
- Postprocessing and structural analysis
"""

__author__ = "Eivind Sonju"
__copyright__ = "Copyright (C) 2017-2019 Verbun AS. All rights reserved."
__version__ = "2.0"

import logging
import multiprocessing
import pickle
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from logutils.queue import QueueListener

import core.code_check as cc
import core.environmental_conditions as ec
import core.meshmagick.hydrostatics as hs
import core.model_set_up as msu
import core.tool_box as tb
from core import nemoh, response, stability
from core.common import SettingsClass
from core.loads import Sea_and_Inertia_Loads
from core.meshmagick.mesh import Mesh
from nemoh_frontend import nemoh_frontend as nf
from pyNemoh_root.pyNemoh.structure import BaseStructure


if __name__ == '__main__':
    ANALYSES_ROOT = Path(r'C:\MOLO_Optimization')
    PARK_LABEL = 'site_01'
    WTG_LABEL = 'wtg_01'

    settings = SettingsClass(ANALYSES_ROOT, PARK_LABEL, WTG_LABEL)

    settings.create_model = True
    settings.calc_intact_stability = True
    settings.run_nemoh = False
    settings.postprocessing = True

    h5_bs = BaseStructure()

    settings.load_cases = {  # 121, np.pi / 15, np.pi
            "num_wave_frequencies": 41,
            "min_wave_frequencies": 2 * np.pi / 27,  # (rad/s)
            "max_wave_frequencies": 2 * np.pi / 4,
            "num_wave_directions": 2,
            "min_wave_directions": 0,  # deg
            "max_wave_directions": 90,
    }
    # TODO: Allow for none equidistant frequencies
    settings.case_label = 'molo_model'
    settings.set_file_structure_and_report()
    settings.simulation_dir = str(settings.fio.nemoh_root)
    settings.save_job_settings()

    settings.mesh_name = None
    settings.use_dipols = False
    settings.use_symmmetri = True

    settings.do_equilibriate = True

    settings.thin_panel_offset = 0.5

    if settings.create_model:
        print('\n--------------------------------------------------------------------------------------------')
        print('CREATE MODEL')
        print('--------------------------------------------------------------------------------------------')

        unit_model, hs_floater = msu.init_models(settings)
        with open(settings.fio.data_io_dir.joinpath('unit_model.pkl'), 'wb') as f:
            pickle.dump(unit_model, f)
        with open(settings.fio.data_io_dir.joinpath('hs_floater.pkl'), 'wb') as f:
            pickle.dump(hs_floater, f)

        settings.thin_panels = []

        m, k = tb.load_M_and_K(settings.fio.data_io_dir)

        print('\nEigenvalue solution WITHOUT added mass (given as lambda^0.5)')
        print('\nMass matrix:')
        tb.matprint(m)
        print('\nStiffness matrix:')
        tb.matprint(k)
        print('')
        tb.eigenvalprint(m, k)

    else:
        with open(settings.fio.data_io_dir.joinpath('unit_model.pkl'), 'rb') as f:
            unit_model = pickle.load(f)
        with open(settings.fio.data_io_dir.joinpath('hs_floater.pkl'), 'rb') as f:
            hs_floater = pickle.load(f)

    if settings.create_model and settings.calc_intact_stability:
        print('\n--------------------------------------------------------------------------------------------')
        print('STABILITY ANALYSIS')
        print('--------------------------------------------------------------------------------------------')
        area_ratio = stability.intact_stability(settings, hs_floater)

    if settings.create_model and settings.run_nemoh:
        print('\n--------------------------------------------------------------------------------------------')
        print('NEMOH ANALYSIS')
        print('--------------------------------------------------------------------------------------------')
        settings.remove_old_db()
        queue = multiprocessing.Queue(-1)
        ql = QueueListener(queue, *logging.getLogger().handlers)
        ql.start()
        nf.run(settings._job_data['analysis'], queue)
        ql.stop()

    if settings.postprocessing:

        loads = Sea_and_Inertia_Loads(settings)
        tran_fun = response.ResponseModel(settings, loads)
        panel_cc = cc.Panel(settings)

        imass, ipanel = tran_fun.get_flange_panel_index()
        f_part1 = tran_fun.assemble_forces(imass, ipanel, moment_ref_point=[0, 0, 0])

        # mesh = Mesh(loads.pd.ppoints, loads.pd.ppanels[ipanel])
        # mesh.show()

        # Get section forces
        section_point = [3.5, 0, 0]
        section_normal = [1, 0, 0]
        imass, ipanel = tran_fun.get_section_index(section_point, section_normal)
        f_sec1 = tran_fun.assemble_forces(imass, ipanel, moment_ref_point=section_point)

        mesh = Mesh(loads.pd.ppoints, loads.pd.ppanels[ipanel])
        # mesh.show()

        sig_dyn_tot = panel_cc.axial_stress(f_sec1['Dynamic']['Total'])
        sig_stat_tot = panel_cc.axial_stress(f_sec1['Static']['Total'])

        p_dyn_lat = panel_cc.lateral_pressure(f_part1['Dynamic']['Total'])
        p_stat_lat = panel_cc.lateral_pressure(f_part1['Static']['Total'])

        hs = 12
        tp = 14
        # gamma = env.gamma(hs, tp)
        # sig_r = np.abs(sig_dyn_tot ** 2) * env.s_jonswap(hs=hs, wp=2 * np.pi / tp, w=loads.w, gamma=gamma)[:, np.newaxis]
        # dw = loads.w[1] - loads.w[0]
        # sig_r_m0 = sum(sig_r) * dw
        # tz = env.tp2tz(tp, gamma)
        # nz = 3 * 3600 / tz
        # sig_r_max = np.sqrt(sig_r_m0)*(np.sqrt(2*np.log(nz))+0.5772/np.sqrt(2*np.log(nz)))

        stwc1 = ec.Short_Term_Wave_Conditions(hs=hs, tp=tp)

        print(f'\nExpected largest maximum dynamic normal stress for Hs = {hs:4.1f} m and Tp = {tp:4.1f} s')
        for i in range(loads._nbeta):
            print(f'\tWavedir {loads._beta[i] * 180 / np.pi:5.1f} deg: {stwc1.expected_largest_maximum(sig_dyn_tot[:, i], loads.w) / 10 ** 6:6.1f} MPa')
        print(f'\nStatic stress: {sig_stat_tot / 10 ** 6:1.1f} MPa')

        print(f'\nExpected largest maximum dynamic lateral pressure for Hs = {hs:4.1f} m and Tp = {tp:4.1f} s')
        for i in range(loads._nbeta):
            print(f'\tWavedir {loads._beta[i] * 180 / np.pi:5.1f} deg: {stwc1.expected_largest_maximum(p_dyn_lat[:, i], loads.w) / 10 ** 3:6.1f} kPa')
        print(f'\nStatic lateral force: {p_stat_lat / 10 ** 3:1.1f} kPa')

        print('\n\nUtilizations')
        gamma_m = 1.15
        load_factor = 1.3
        sigma_y = 235000000 / 1.15
        bc = 'pinned'
        dpu = panel_cc.dynamic_panel_utilization(sigma_y, bc, f_sec1['Dynamic']['Total'], f_part1['Dynamic']['Total'],
                                                 hs, tp, freq=loads.w, pos_y_side=True)
        for i in range(loads._nbeta):
            print(f'\tWavedir {loads._beta[i] * 180 / np.pi:5.1f} deg: {dpu[i] * load_factor:6.2f}')

        dpu = panel_cc.dynamic_panel_utilization(sigma_y, bc, f_sec1['Dynamic']['Total'], f_part1['Dynamic']['Total'],
                                                 hs, tp, freq=loads.w, pos_y_side=False)
        for i in range(loads._nbeta):
            print(f'\tWavedir {loads._beta[i] * 180 / np.pi:5.1f} deg: {dpu[i] * load_factor:6.2f}')

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
            print('\nEigenvalue solution WITH added mass')
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
            h = tran_fun.calc_rao_linear(idir)
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

    settings._report._doc.generate_pdf(clean_tex=False)
