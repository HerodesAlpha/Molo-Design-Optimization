__author__ = "Eivind Sonju"
__copyright__ = "Copyright (C) 2017-2019 Verbun AS. All rights reserved."
__version__ = "2.0"

import json
from pathlib import Path
import pickle
import h5py
#import matplotlib.pyplot as plt
import numpy as np
from logutils.queue import QueueListener
import multiprocessing
import logging
import helper
import calculations
import model_set_up as msu
import nemoh
from common import FileIOClass
from common import SettingsClass

from MOLO_Nemoh import nemoh_frontend as nf
from MOLO_Nemoh import settings
from pyNemoh.structure import H5_STRUCTURE
from pyNemoh.postprocessor import read_results

import stability

if __name__ == '__main__':
    run_nemoh = False
    create_model = False
    calc_gz = False
    postprocessing = True

    ANALYSES_ROOT = Path(r'C:\analyses')
    PARK_LABEL='site_01'
    WTG_LABEL = 'wtg_01'

    settings = SettingsClass(ANALYSES_ROOT,PARK_LABEL,WTG_LABEL)
    settings.floater_data = {
                                "Type":                         "OY",
                                "Central column diameter":      10,
                                "Central column thickness":     0.04,
                                "Draught":                      0,
                                "Gap factor":                   0.8,
                                "Lower flange thickness":       0.04,
                                "Number of radial columns":     2,
                                "Radial column diameter":       11,
                                "Radial column thickness":      0.04,
                                "Radial height":                18,
                                "Upper flange thickness":       0.04
    }
    settings.case_label = 'floater_data'
    settings.set_file_structure()

    settings.simulation_dir = str(settings.fio.nemoh_root)

    fio = settings.fio


    if create_model:



        # NEMOH_DOF = model_data["Analyses parameters"]["Degrees of Freedom"]
        # NEMOH_DIR = model_data["Analyses parameters"]["Number of wave directions, Min and Max (degrees)"]
        # RHO_SW = model_data["Analyses parameters"]["Density of sea water"]
        # WATER_DEPTH = model_data["Analyses parameters"]["Water depth"]
        # OMEGA_NEMOH_INP = model_data["Analyses parameters"]["Number of wave frequencies, Min, and Max (rad/s)"]
        # SYM = model_data["Analyses parameters"]["Use symmetri"]

        unit_model, hs_floater = msu.launch_dipole(settings)
        pickle.dump(unit_model, open(fio.data_io_dir.joinpath('unit_model.pkl'), 'wb'))
        pickle.dump(hs_floater, open(fio.data_io_dir.joinpath('hs_floater.pkl'), 'wb'))

        # hydro_mesh_symmetri, mesh_file = nemoh.mesh(fio, hs_floater, SYM)
        #
        # hydro_mesh_symmetri.show()


    else:
        CASE_NUMBER = 2
        SYM = 1
        unit_model = pickle.load(open(fio.data_io_dir.joinpath('unit_model.pkl'), 'rb'))
        hs_floater = pickle.load(open(fio.data_io_dir.joinpath('hs_floater.pkl'), 'rb'))

    if create_model and calc_gz:
        stability.gz_curve(fio, hs_floater)

    if create_model and run_nemoh:
        queue = multiprocessing.Queue(-1)
        ql = QueueListener(queue, *logging.getLogger().handlers)
        ql.start()
        nf.run(settings._job_data['analysis'], queue)
        ql.stop()
        # nemoh.runNemoh(fio, hydro_mesh_symmetri, mesh_file, NEMOH_DIR, RHO_SW, WATER_DEPTH, OMEGA_NEMOH_INP,
        #                NEMOH_DOF)

    if postprocessing:
        structure = H5_STRUCTURE()
        # Read results, perform postprocessing and write pdf
        with h5py.File(fio.nemoh_root.joinpath('db.hdf5'), "r") as hdf5_db:
            for key in hdf5_db['results'].keys():
                print(key)  # Names of the groups in HDF5 file.
            #beta = hdf5_db['/results/case/beta'].value
            #print(beta)
            #results=read_results(hdf5_db)
            print(hdf5_db['results']['fk_pressure_raw'][0])

        hdp = calculations.TransferFunctions(fio)

