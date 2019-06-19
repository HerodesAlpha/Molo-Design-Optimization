__author__ = "Eivind Sonju"
__copyright__ = "Copyright (C) 2017-2019 Verbun AS. All rights reserved."
__version__ = "2.0"

import logging
import multiprocessing
import pickle
from pathlib import Path
import scipy
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
import warnings
from meshmagick.mesh import Mesh

if __name__ == '__main__':
    ANALYSES_ROOT = Path(r'C:\analyses')
    PARK_LABEL = 'site_01'
    WTG_LABEL = 'wtg_01'

    settings = SettingsClass(ANALYSES_ROOT, PARK_LABEL, WTG_LABEL)

    settings.create_model = True
    settings.calc_gz = False
    settings.run_nemoh = False
    settings.postprocessing = True

    h5_bs = BaseStructure()

    filling_ratio = [0, 0, 0]

    settings.floater_data = {
            "Type"                    : "OY",
            "Central column diameter" : 9,
            "Central column thickness": 0.04,
            "Draught"                 : 0,
            "Gap factor"              : 0.8,
            "Lower flange thickness"  : 0.04,
            "Number of radial columns": 3,
            "Radial column diameter"  : 9,
            "Radial column thickness" : 0.04,
            "Radial height"           : 14,
            "Upper flange thickness"  : 0.04,
            "Ballast filling ratio"   : [0,
                                         filling_ratio
                                         ]
    }
    settings.load_cases = {  # 121, np.pi / 15, np.pi
            "num_wave_frequencies": 41,
            "min_wave_frequencies": 2 * np.pi / 25,  # (rad/s)
            "max_wave_frequencies": 2 * np.pi / 5,
            "num_wave_directions" : 1,
            "min_wave_directions" : 0,  # deg
            "max_wave_directions" : 90,
    }

    settings.case_label = 'floater_data'
    settings.set_file_structure()

    settings.mesh_name = None
    settings.use_dipols = False
    settings.use_symmmetri = True

    settings.simulation_dir = str(settings.fio.nemoh_root)

    settings.do_equilibriate = True
    fio = settings.fio

    if settings.create_model:

        unit_model, hs_floater = msu.init_models(settings)
        pickle.dump(unit_model, open(fio.data_io_dir.joinpath('unit_model.pkl'), 'wb'))
        pickle.dump(hs_floater, open(fio.data_io_dir.joinpath('hs_floater.pkl'), 'wb'))

        settings.thin_panels = []

        m, k = tb.load_M_and_K(settings.fio.data_io_dir)
        # m=m[0:5,0:5]
        # k=k[0:5,0:5]
        # Set values close to zero to zero
        for i in range(6):
            for j in range(6):
                if abs(m[i, j]) < 0:
                    m[i, j] = 0
                if abs(k[i, j]) < 0:
                    k[i, j] = 0
                #if (i == 0 and j == 0) or (i == 1 and j == 1) or (i == 5 and j == 5):
                #    k[i, j] = 1

        print('\nEigenvalue sollution WITHOUT added mass')
        print('\nMass matrix (heave, pitch and roll):')
        tb.matprint(m)
        print('\nStiffness matrix (heave, pitch and roll):')
        tb.matprint(k)
        print('')
        tb.eigenvalprint(m, k)

    else:
        unit_model = pickle.load(open(settings.fio.data_io_dir.joinpath('unit_model.pkl'), 'rb'))
        hs_floater = pickle.load(open(settings.fio.data_io_dir.joinpath('hs_floater.pkl'), 'rb'))

    if settings.create_model and settings.calc_gz:
        stability.gz_curve(fio, hs_floater)

    if settings.create_model and settings.run_nemoh:
        queue = multiprocessing.Queue(-1)
        ql = QueueListener(queue, *logging.getLogger().handlers)
        ql.start()
        nf.run(settings._job_data['analysis'], queue)
        ql.stop()
        # nemoh.runNemoh(fio, hydro_mesh_symmetri, mesh_file, NEMOH_DIR, RHO_SW, WATER_DEPTH, OMEGA_NEMOH_INP,
        #                NEMOH_DOF)

    if settings.postprocessing:

        hdp = calculations.TransferFunctions(settings)

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
        irad = 2
        # nprob = len(NEMOH_DIR) + sum(NEMOH_DOF)
        iprob = 2
        # nfreq  = len(w)
        # problem = (ifreq - 1) * nprob + iprob

        # print('Problem: {}'.format(problem))
        NEMOH_DOF = [1, 1, 1, 1, 1, 1]

        if False:
            # hdp.show_pressure(ifreq, idir, pressure_type='Froude-Krylof',axis=2)
            # hdp.show_pressure(ifreq, idir, pressure_type='Diffraction', axis=2)
            hdp.show_pressure(ifreq, irad, pressure_type='Radiation', axis=2)
            pass

        # print(hdp.ma[0,:,0,0])
        # print(hdp._fe_amp[0,:,dof])

        idof = np.array([i for i, x in enumerate(NEMOH_DOF) if x])

        print('Eigenvalue sollution WITH added mass')
        tb.eigenvalprint(hdp.m + hdp.ma[ifreq, :, :], hdp.k)
        #
        # warnings.filterwarnings("ignore", category=RuntimeWarning)
        # lamb = scipy.linalg.eigvalsh(hdp.k, hdp.m + hdp.ma[ifreq, :, :])
        # vT = 2 * np.pi / np.sqrt(lamb)
        # warnings.filterwarnings("default")

        # for x, T in enumerate(vT):
        #     print('Eigenval {}:\t{:5.1f} s'.format(x + 1, T))
        # available_dofs=np.array([2,4])
        # k=hdp.k[available_dofs,available_dofs]
        # m= hdp.m[available_dofs,available_dofs]+hdp.ma[ifreq,available_dofs,available_dofs]
        # print(np.linalg.eig(k,m))

        if True:
            print('\n')

            print('\nRadiation damping:')
            tb.matprint(hdp.c_hyd[ifreq])
            print('\nWater plane stiffness:')
            tb.matprint(hdp.k)
            print('\nStatic mass [tonne]:')
            tb.matprint(hdp.m / 1000)
            print('\nAdded mass [tonne]:')
            tb.matprint(hdp.ma[ifreq] / 1000)
            print('\nExcitation force:')
            tb.matprint(np.abs(hdp.fe[ifreq, idir, :]))
            # tmp = nemoh.get_section_forces(fio, problem, [-100,0,0], [1,0,0], sym=SYM)
            # print('\nSection force:\n{}'.format(np.abs(tmp)))

            f_fk = hdp.p2f('Froude-Krylof', ifreq, idir)
            f_diff = hdp.p2f('Diffraction', ifreq, idir)
            f_exc = f_fk + f_diff

            print('\n')
            tb.matprint(np.abs(nemoh.get_section_values(f_exc, hdp.pd.ppanel_centers, [0, 0, 0], [1, 0, 0])))
        # np.set_printoptions(precision=3)
        idof = 2
        idir = 0

        # rao = hdp.get_rao(idof, idir)
        h = hdp.get_h(idir)

        fig, axs = plt.subplots(3, 2)
        x = 2 * np.pi / hdp.w


        def y(idof):
            abs(h[:, idof])


        axs[0, 0].plot(x, abs(h[:, 2]), 'tab:orange')
        axs[0, 0].set_title('Rao heave')
        axs[0, 1].plot(x, abs(h[:, 4]), 'tab:orange')
        axs[0, 1].set_title('Rao pitch')
        axs[1, 0].plot(x, hdp.ma[:, 2, 2] + hdp.m[2, 2], 'tab:green')
        axs[1, 0].set_title('m+ma heave')
        axs[1, 1].plot(x, hdp.ma[:, 4, 4] + hdp.m[4, 4], 'tab:green')

        axs[1, 1].set_title('m+ma pitch')
        axs[2, 0].plot(x, abs(hdp._fe[:, idir, 2]), 'tab:blue')
        axs[2, 0].set_title('fe heave')
        axs[2, 1].plot(x, abs(hdp._fe[:, idir, 4]), 'tab:blue')
        axs[2, 1].set_title('fe pitch')

        # plt.plot(2 * np.pi / hdp.w, abs(h[:, idof]))
        # plt.show()
        #
        # plt.plot(2 * np.pi / hdp.w, abs(hdp.fe[:, idir, idof]))
        plt.show()

        # np.set_printoptions(precision=3)

        # print(abs(sum(nemoh.p2f(hdp._p['Hydro static'], hdp.pd))) / 9.81)

        f_static = hdp.p2f('Hydro_static')
        sum_f33_static = np.abs(nemoh.get_section_values(f_static, hdp.pd.ppanel_centers, [-1000, 0, 0], [1, 0, 0]))[2]

        print('\nTotal hydrostatic force:\t{:5.2f} tonne'.format(sum_f33_static / 9810))

        part_list = unit_model.get_parts()
        print('\nSum of all parts:\t{:5.2f} tonne'.format(sum([part.mass for part in part_list]) / 1000))

        f_grav = np.asarray([[0, 0, part.mass] for part in part_list]) * settings.grav
        mass_centers = np.asarray([-part.reduction_point for part in part_list])
        f_grav_s = np.abs(
                nemoh.get_section_values(f_grav, mass_centers, [0, 1, 0], [0, 1, 0], moment_ref='section'),
        )
        print('\nGravity force Fz\t{:5.2f} MN'.format(f_grav_s[2] / 1000000))
        print('Gravity force Mx\t{:5.2f} MNm'.format(f_grav_s[3] / 1000000))
        print('Gravity force My\t{:5.2f} MNm'.format(f_grav_s[4] / 1000000))

        # print(hdp.ma_zero)
        # print(hdp.ma_inf)

    # this_mesh = Mesh(hdp.pd.ppoints, hdp.pd.ppanels)
    # this_mesh.show()
