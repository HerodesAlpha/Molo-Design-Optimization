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

    settings.create_model = False
    settings.calc_intact_stability = False
    settings.run_nemoh = False
    settings.postprocessing = True

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
    settings.case_label = 'molo_model'
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

        with open(settings.fio.data_io_dir.joinpath('unit_model.pkl'), 'rb') as f:
            unit_model = pickle.load(f)
            unit_model.move_reduction_point([0, 0, 10]) # Why?
            unit_model.print_vector_matrix_global()
            unit_model.move_reduction_point([0, 0, 10]) # Why?
            unit_model.print_vector_matrix_global()
