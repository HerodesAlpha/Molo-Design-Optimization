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
import tool_box as tb
import calculations
import model_set_up as msu
import nemoh
import stability
from MOLO_Nemoh import nemoh_frontend as nf
from common import SettingsClass
from pyNemoh.structure import BaseStructure
from loads import Sea_and_Inertia_Loads

if __name__ == '__main__':
    ANALYSES_ROOT = Path(r'C:\analyses')
    PARK_LABEL = 'site_01'
    WTG_LABEL = 'wtg_01'

    settings = SettingsClass(ANALYSES_ROOT, PARK_LABEL, WTG_LABEL)

    settings.create_model = True
    settings.calc_intact_stability = True
    settings.run_nemoh = False
    settings.postprocessing = False

    h5_bs = BaseStructure()

    filling_ratio = [0.0, 0.0, 0.0]

    # settings.floater_data = {
    #         "Type"                    : "OY",
    #         "Central column diameter" : 7.0,
    #         "Central column thickness": 0.04,
    #         "Draught"                 : 0,
    #         "Gap factor"              : 0.8,
    #         "Lower flange thickness"  : 0.08,
    #         "Number of radial columns": 3,
    #         "Radial column diameter"  : 8.5,
    #         "Radial column thickness" : 0.04,
    #         "Radial height"           : 15,
    #         "Upper flange thickness"  : 0.08,
    #         "Ballast filling ratio"   : [0,
    #                                      filling_ratio
    #                                      ],
    #         "Thin panel offset"       : 0.2
    # }
    settings.load_cases = {  # 121, np.pi / 15, np.pi
            "num_wave_frequencies": 41,
            "min_wave_frequencies": 2 * np.pi / 27,  # (rad/s)
            "max_wave_frequencies": 2 * np.pi / 4,
            "num_wave_directions" : 2,
            "min_wave_directions" : 0,  # deg
            "max_wave_directions" : 90,
    }
    # TODO: Allow for none equidistant frequencies
    settings.case_label = 'stab1'
    settings.set_file_structure_and_report()
    settings.simulation_dir = str(settings.fio.nemoh_root)
    settings.save_job_settings()

    settings.mesh_name = None
    settings.use_dipols = False
    settings.use_symmmetri = True

    settings.do_equilibriate = True
    #fio = settings.fio

    settings.thin_panel_offset = 0.5

    if settings.create_model:
        print('\n--------------------------------------------------------------------------------------------')
        print('CREATE MODEL')
        print('--------------------------------------------------------------------------------------------')

        unit_model, hs_floater = msu.init_models(settings)
        pickle.dump(unit_model, open(settings.fio.data_io_dir.joinpath('unit_model.pkl'), 'wb'))
        pickle.dump(hs_floater, open(settings.fio.data_io_dir.joinpath('hs_floater.pkl'), 'wb'))

        settings.thin_panels = []

        m, k = tb.load_M_and_K(settings.fio.data_io_dir)

        print('\nEigenvalue sollution WITHOUT added mass (given as lambda^0.5)')
        print('\nMass matrix:')
        tb.matprint(m)
        print('\nStiffness matrix:')
        tb.matprint(k)
        print('')
        tb.eigenvalprint(m, k)

    else:
        unit_model = pickle.load(open(settings.fio.data_io_dir.joinpath('unit_model.pkl'), 'rb'))
        hs_floater = pickle.load(open(settings.fio.data_io_dir.joinpath('hs_floater.pkl'), 'rb'))

    if settings.create_model and settings.calc_intact_stability:
        print('\n--------------------------------------------------------------------------------------------')
        print('STABILITY ANALYSIS')
        print('--------------------------------------------------------------------------------------------')
        stability.intact_stability(settings, hs_floater)

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
        tran_fun = calculations.TransferFunctions(settings, loads)
        env = calculations.Environment(settings)
        f_sf1 = tran_fun.section_forces([1, 0, 0], [1, 0, 0])
        sig_dyn_tot = tran_fun.section_stress(f_sf1['Dynamic']['Total'])
        sig_stat_tot = tran_fun.section_stress(f_sf1['Static']['Total'])


        hs = 12
        tp = 14
        gamma = env.gamma(hs, tp)
        sig_r = np.abs(sig_dyn_tot ** 2) * env.s_jonswap(hs=hs, wp=2 * np.pi / tp, w=loads.w, gamma=gamma)[:, np.newaxis]
        dw = loads.w[1] - loads.w[0]
        sig_r_m0 = sum(sig_r) * dw
        tz = env.tp2tz(tp, gamma)
        nz = 3 * 3600 / tz
        sig_r_max = np.sqrt(sig_r_m0)*(np.sqrt(2*np.log(nz))+0.5772/np.sqrt(2*np.log(nz)))

        print('\nExpected largest maximum dynamic normal stress for Hs = {:4.1f} m and Tp = {:4.1f} s'.format(hs,tp))
        for i in range(loads._nbeta):
            print('\tWavedir {:5.1f} deg: {:6.1f} MPa'.format(loads._beta[i]*180/np.pi, sig_r_max[i] / 10 ** 6))
        print('\nStatic stress: {:1.1f} MPa'.format(sig_stat_tot / 10 ** 6))
        plt.plot(2 * np.pi / loads.w, sig_r)

        #plt.show()

        SELECT_DOF = 1
        SELECT_AXIS = 2
        SELECT_FREC = 1

        # print('\nSection forces')
        # for i in range(6):
        #     print('\tDOF{}: {:8.1f}'.format(i + 1, np.sqrt(
        #         abs(sum(hydro.spec_response(hdp.f_sec[:, i][::3], hdp.w, hs=10, wp=2 * np.pi / 14))))))
        # write_report(root, hdp, case_label)

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


        if True:
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

            print('\n')
            tb.matprint(np.abs(nemoh.get_section_values(f_exc, loads.pd.ppanel_centers, [0, 0, 0], [1, 0, 0])))
        # np.set_printoptions(precision=3)

        idir = 0

        # --------------------------------------------------------------------------------------------------------------
        # PLOT RESULTS
        # --------------------------------------------------------------------------------------------------------------
        if True:
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
        #
        # f = hdp.section_forces()
        #
        #
        # plt.plot(2 * np.pi / hdp.w, abs(f[:, 0, 4]))
        # plt.show()

        # np.set_printoptions(precision=3)

        # print(abs(sum(nemoh.p2f(hdp._p['Hydro static'], hdp.pd))) / 9.81)

        #
        f = tran_fun.section_forces([3.5, 0, 0], [1, 0, 0])

        f_dyn = f['Dynamic']
        w = 2 * np.pi / loads.w

        fig, axs = plt.subplots(2, 3)

        dir = 2
        for i, dof in enumerate([1, 5]):
            for key in f_dyn:
                if not key is 'Total':
                    axs[i, 0].plot(w, abs(f_dyn[key][:, dir, dof]), label=key)

            axs[i, 0].legend()

            y = abs(f_dyn['Total'][:, dir, dof])
            axs[i, 1].plot(w, y, label='Total')
            axs[i, 1].legend()

            axs[i, 2].plot(w, abs(tran_fun._rao[:, dir, dof]), label='RAO')
            axs[i, 2].axis([5, 15, 0, 0.01])
            axs[i, 2].legend()
            # axs[1, i].plot(w, abs(f_varying_buoyancy[:, dof]), 'tab:blue', label='Varying Buoyancy')
            # axs[1, i].plot(w, abs(f_inertia[:, dof]), 'tab:green', label='Inertia')
            # axs[1, i].set_title('Varying forces')
            # axs[1, i].legend()

        #plt.show()

        # Calc velocity and acc (not needed yet)
        # panel_vel = panel_pos * 1j * hdp.w(ifreq)
        # panel_acc = panel_pos * hdp.w(ifreq)**2

        # xyz = np.zeros([hdp._pd._npanel, 3])
        # xyz[:, 2] = 1
        # nemoh_mesh = Mesh(hdp._pd.ppoints, hdp._pd.ppanels)
        # h = force.show_force(nemoh_mesh, hdp._pd.ppanel_centers, np.imag(panel_pos * xyz))
        # h.show()

        # print('\nSum of all parts:\t{:5.2f} tonne'.format(sum([part.mass for part in part_list]) / 1000))

        # print('\nFz\t{: 7.2f} MN'.format(f_gravity[2] / 1000000))
        # print('Mx\t{: 7.2f} MNm'.format(f_gravity[3] / 1000000))
        # print('my\t{: 7.2f} MNm'.format(f_gravity[4] / 1000000))

        # print(hdp.ma_zero)
        # print(hdp.ma_inf)

        # this_mesh = Mesh(hdp.pd.ppoints, hdp.pd.ppanels)
        # this_mesh.show()
    settings._report._doc.generate_pdf(clean_tex=False)