__author__ = "Eivind Sonju"
__copyright__ = "Copyright (C) 2017-2019 Verbun AS. All rights reserved."
__version__ = "2.0"

import json
import os
from pathlib import Path
import pickle

# import matplotlib.pyplot as plt
import numpy as np
from logutils.queue import QueueListener
import multiprocessing
import logging
import helper
import calculations
import model_set_up as msu
import nemoh

from MOLO_Nemoh import nemoh_frontend as nf

from pathlib import Path


class PhysicalQuantities():
    # This is the only place allowed to put these quantities
    def __init__(self):
        self._rho_sw = 1025
        self._grav = 9.81
        self._rho_st = 7850

    @property
    def rho_sw(self):
        return self._rho_sw

    @property
    def grav(self):
        return self._grav

    @property
    def rho_st(self):
        return self._rho_st



class FileIOClass(object):
    def __init__(self, root_dir, case_number=None):
        if case_number == None:
            i = 0
            while 1:
                i += 1
                if not root_dir.joinpath('case{:04d}'.format(i)).exists():
                    break
            self._case_label = 'case{:04d}'.format(i)
        else:
            self._case_label = 'case{:04d}'.format(case_number)

        self._case_dir = root_dir.joinpath(self._case_label)

        self._nemoh_root = self._case_dir.joinpath('nemoh')

        self._gmsh_root = self._case_dir.joinpath('gmsh')

        path_to_gmsh = r'C:\Users\eison\OneDrive - Verbun AS\Divisions\Offshore Wind\Library\Software\Bin\gmsh-4.2.2-Windows64\gmsh.exe'
        if Path(path_to_gmsh).exists():
            self._gmsh_exe = path_to_gmsh

        self._data_io_dir = self._case_dir.joinpath('data_io')

        self._nemoh_results_dir = self._nemoh_root.joinpath('results')
        self._nemoh_mesh_dir = self._nemoh_root.joinpath('mesh')

        self._stability_dir = self._case_dir.joinpath('stability')

        if case_number == None:
            self._nemoh_root.mkdir(parents=True, exist_ok=False)
            self._nemoh_results_dir.mkdir(parents=True, exist_ok=False)
            self._nemoh_mesh_dir.mkdir(parents=True, exist_ok=False)
            self._gmsh_root.mkdir(parents=True, exist_ok=False)
            self._data_io_dir.mkdir(parents=True, exist_ok=False)
            self._stability_dir.mkdir(parents=True, exist_ok=False)

        self._templates_dir = Path(os.getcwd()).joinpath('templates')
        if not self._templates_dir.exists():
            print('Templates folder missing')

    @property
    def case_dir(self):
        return self._case_dir

    @property
    def nemoh_dir(self):
        return self._nemoh_root

    @property
    def gmsh_dir(self):
        return self._gmsh_root

    @property
    def data_io_dir(self):
        return self._data_io_dir

    @property
    def templates_dir(self):
        return self._templates_dir

    @property
    def nemoh_results(self):
        return self._nemoh_results_dir

    @property
    def nemoh_root(self):
        return self._nemoh_root

    @property
    def nemoh_results_dir(self):
        return self._nemoh_results_dir

    @property
    def nemoh_mesh_dir(self):
        return self._nemoh_mesh_dir

    @property
    def stability_dir(self):
        return self._stability_dir

    @property
    def gmsh_exe(self):
        return self._gmsh_exe


class SettingsClass(PhysicalQuantities,object):
    def __init__(self, fio):
        super().__init__()
        self._json_list = ['park', 'rna','tower', 'floater', 'analysis']

        self.job_data = dict()
        self._fio = fio

        # Collect template data
        for item in self._json_list:
            with open(self._fio.templates_dir.joinpath('{}_template.json'.format(item)), 'r') as f:
                self.job_data[item] = json.loads(f.read())

        # Save updated template to template dir
        for item in self._json_list:
            with open(self._fio.templates_dir.joinpath('{}_template.json'.format(item)), 'w') as f:
                f.write(json.dumps(self.job_data[item], indent=4, sort_keys=True))

    def save_job_settings(self):
        # Save updated template to analysis directory
        for item in self._json_list:
            with open(self._fio.data_io_dir.joinpath('{}_template.json'.format(item)), 'w') as f:
                f.write(json.dumps(self.job_data[item], indent=4, sort_keys=True))

    @property
    def mesh_file(self):
        return self.job_data['analysis']['simulations']['sim01']['floating_bodies']['sim01.dat']['mesh_file']

    @mesh_file.setter
    def mesh_file(self, val):
        self.job_data['analysis']['simulations']['sim01']['floating_bodies']['sim01.dat']['mesh_file'] = val
        self.save_job_settings()

    @property
    def draught(self):
        return self.job_data['floater']['Draught']

    @draught.setter
    def draught(self, val):
        self.job_data['floater']['Draught'] = val
        self.save_job_settings()
