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

import calculations
import model_set_up as msu
import nemoh
import stability
from MOLO_Nemoh import nemoh_frontend as nf
from common import SettingsClass
from pyNemoh.structure import BaseStructure


if __name__ == '__main__':
    ANALYSES_ROOT = Path(r'C:\analyses')
    PARK_LABEL = 'site_01'
    WTG_LABEL = 'wtg_01'

    settings = SettingsClass(ANALYSES_ROOT, PARK_LABEL, WTG_LABEL)

    settings.create_model = True
    settings.calc_gz = False
    settings.run_nemoh = True
    settings.postprocessing = True

    h5_bs = BaseStructure()



    settings.floater_data = {
            "Type"                    : "OY",
            "Central column diameter" : 7,
            "Central column thickness": 0.04,
            "Draught"                 : 0,
            "Gap factor"              : 0.8,
            "Lower flange thickness"  : 0.04,
            "Number of radial columns": 3,
            "Radial column diameter"  : 9,
            "Radial column thickness" : 0.04,
            "Radial height"           : 18,
            "Upper flange thickness"  : 0.04
    }
    settings.load_cases = {  # 121, np.pi / 15, np.pi
            "num_wave_frequencies": 1,
            "min_wave_frequencies": 1,  # (rad/s)
            "max_wave_frequencies": 1,
            "num_wave_directions" : 1,
            "min_wave_directions" : 0,  # deg
            "max_wave_directions" : 90,
    }

    # settings.case_label = 'floater_data'
    # settings.case_label = None
    # settings.case_label = 'debug_plate_thin'
    settings.case_label = 'debug2_plate_thin'
    # settings.case_label = 'debug_plate_thick'
    settings.set_file_structure()

    # settings.mesh_name = 'debug_3c'
    settings.mesh_name = settings.case_label
    #settings.mesh_name = 'debug_plate_thin'
    settings.use_dipols = False


    settings.simulation_dir = str(settings.fio.nemoh_root)

    settings.do_equilibriate=False
    settings.draught=6.02

    fio = settings.fio

    if settings.create_model:

        # NEMOH_DOF = model_data["Analyses parameters"]["Degrees of Freedom"]
        # NEMOH_DIR = model_data["Analyses parameters"]["Number of wave directions, Min and Max (degrees)"]
        # RHO_SW = model_data["Analyses parameters"]["Density of sea water"]
        # WATER_DEPTH = model_data["Analyses parameters"]["Water depth"]
        # OMEGA_NEMOH_INP = model_data["Analyses parameters"]["Number of wave frequencies, Min, and Max (rad/s)"]
        # SYM = model_data["Analyses parameters"]["Use symmetri"]

        unit_model, hs_floater = msu.init_models(settings)
        pickle.dump(unit_model, open(fio.data_io_dir.joinpath('unit_model.pkl'), 'wb'))
        pickle.dump(hs_floater, open(fio.data_io_dir.joinpath('hs_floater.pkl'), 'wb'))

        # hydro_mesh_symmetri, mesh_file = nemoh.mesh(fio, hs_floater, SYM)
        #
        # hydro_mesh_symmetri.show()

        settings.thin_panels = []


    else:
        CASE_NUMBER = 2
        SYM = 1
        unit_model = pickle.load(open(fio.data_io_dir.joinpath('unit_model.pkl'), 'rb'))
        hs_floater = pickle.load(open(fio.data_io_dir.joinpath('hs_floater.pkl'), 'rb'))

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

        if 1:
            # hdp.show_pressure(ifreq, idir, pressure_type='Froude-Krylof',axis=2)
            # hdp.show_pressure(ifreq, idir, pressure_type='Diffraction', axis=1)
            hdp.show_pressure(ifreq, irad, pressure_type='Radiation', axis=2)

        # print(hdp.ma[0,:,0,0])
        # print(hdp._fe_amp[0,:,dof])

        idof = np.array([i for i, x in enumerate(NEMOH_DOF) if x])

        for i, x in enumerate(idof):
            w2 = hdp.k[x, x] / (hdp.m[x, x] + hdp.ma[ifreq, x, x])
            print('T{}:\t{:5.1f} s'.format(x + 1, 2 * np.pi / np.sqrt(w2)))

        np.set_printoptions(precision=3)
        idof = 2
        rao = hdp.get_rao(idof, idir)

        plt.plot(2 * np.pi / hdp.w, rao)
        # plt.show()

        np.set_printoptions(precision=3)

        if 1:
            print('\n')

            print('\nRadiation damping:\n{}'.format(hdp._c_hyd[ifreq]))
            print('\nWater plane stiffness:\n{}'.format(hdp.k))
            print('\nStatic mass:\n{}'.format(hdp.m))
            print('\nAdded mass:\n{}'.format(hdp._ma[ifreq]))
            print('\nExcitation force:\n{}'.format(np.abs(hdp._fe[ifreq, idir, :])))
            # tmp = nemoh.get_section_forces(fio, problem, [-100,0,0], [1,0,0], sym=SYM)
            # print('\nSection force:\n{}'.format(np.abs(tmp)))

            f_fk = nemoh.p2f(hdp._pressure['Froude-Krylof'][ifreq, idir, :], hdp.pd)
            f_diff = nemoh.p2f(hdp._pressure['Diffraction'][ifreq, idir, :], hdp.pd)
            f_exc = f_fk + f_diff

            print(np.abs(nemoh.get_section_values(f_exc, hdp.pd.ppanel_centers, [10, 0, 0], [1, 0, 0])))

            print(abs(hdp.p2f(ifreq, pressure_index=0, p='Hydro static')) / 9.81)
            # print(abs(sum(nemoh.p2f(hdp._p['Hydro static'], hdp.pd))) / 9.81)

        part_list = unit_model.get_parts_without_children()
        print(sum([part.mass for part in part_list]))
