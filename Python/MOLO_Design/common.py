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
    def __init__(self, root_dir, park_label, wtg_label, case_label):
        case_parent_dir = root_dir.joinpath(park_label).joinpath(wtg_label)
        if case_label == None:
            i = 0
            while 1:
                i += 1
                if not case_parent_dir.joinpath('case{:04d}'.format(i)).exists():
                    break
            self._case_label = 'case{:04d}'.format(i)
        else:
            self._case_label = case_label

        self._case_dir = case_parent_dir.joinpath(self._case_label)

        self._nemoh_root = self._case_dir.joinpath('nemoh')
        self._gmsh_root = self._case_dir.joinpath('gmsh')
        self._data_io_dir = self._case_dir.joinpath('data_io')
        self._nemoh_results_dir = self._nemoh_root.joinpath('results')
        self._nemoh_mesh_dir = self._nemoh_root.joinpath('mesh')
        self._stability_dir = self._case_dir.joinpath('stability')

        self._templates_dir = Path(os.getcwd()).joinpath('templates')

        path_to_gmsh = r'C:\Users\eison\OneDrive - Verbun AS\Divisions\Offshore Wind\Library\Software\Bin\gmsh-4.2.2-Windows64\gmsh.exe'
        if Path(path_to_gmsh).exists():
            self._gmsh_exe = path_to_gmsh

        self.create_dir()

    def create_dir(self):
        exist_ok_bool = True
        self._nemoh_root.mkdir(parents=True, exist_ok=exist_ok_bool)
        self._nemoh_results_dir.mkdir(parents=True, exist_ok=exist_ok_bool)
        self._nemoh_mesh_dir.mkdir(parents=True, exist_ok=exist_ok_bool)
        self._gmsh_root.mkdir(parents=True, exist_ok=exist_ok_bool)
        self._data_io_dir.mkdir(parents=True, exist_ok=exist_ok_bool)
        self._stability_dir.mkdir(parents=True, exist_ok=exist_ok_bool)

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


class SettingsClass(PhysicalQuantities, object):
    def __init__(self, analyses_root, park_label, wtg_label):
        super().__init__()

        self._json_list = ['park', 'rna', 'tower', 'floater', 'analysis']

        self._job_data = dict()

        self._analyses_root = analyses_root
        self._park_label = park_label
        self._wtg_label = wtg_label
        self._case_label = None

        self._fio = None

        # Collect template data
        for item in self._json_list:
            with open(Path(os.getcwd()).joinpath('templates').joinpath('{}_template.json'.format(item)), 'r') as f:
                self._job_data[item] = json.loads(f.read())

        # Save updated template to template dir
        for item in self._json_list:
            with open(Path(os.getcwd()).joinpath('templates').joinpath('{}_template.json'.format(item)), 'w') as f:
                f.write(json.dumps(self._job_data[item], indent=4, sort_keys=True))

    def set_file_structure(self):
        self._fio = FileIOClass(self._analyses_root, self._park_label, self._wtg_label, self._case_label)
        print(self._job_data['analysis']['simulations'])
        self._job_data['analysis']['simulations']['sim01']['simulation_dir'] = str(self._fio.nemoh_root)
        self.save_job_settings()

    def save_job_settings(self):
        # Save updated settings to analysis directory
        for item in self._json_list:
            with open(self._fio.data_io_dir.joinpath('{}.json'.format(item)), 'w') as f:
                f.write(json.dumps(self._job_data[item], indent=4, sort_keys=True))

    def load_job_settings(self):
        # Load current settings from analysis directory
        for item in self._json_list:
            with open(self._fio.data_io_dir.joinpath('{}.json'.format(item)), 'r') as f:
                self._job_data[item] = json.loads(f.read())

    # @property
    # def simulation_dir(self):
    #     self.load_job_settings()
    #     return self._job_data['analysis']['simulations']['sim01']['simulation_dir']
    #
    # @simulation_dir.setter
    # def simulation_dir(self, val):
    #     self._job_data['analysis']['simulations']['sim01']['simulation_dir'] = val
    #     self.save_job_settings()

    @property
    def mesh_file(self):
        self.load_job_settings()
        return self._job_data['analysis']['simulations']['sim01']['floating_bodies']['sim01.dat']['mesh_file']

    @mesh_file.setter
    def mesh_file(self, val):
        self._job_data['analysis']['simulations']['sim01']['floating_bodies']['sim01.dat']['mesh_file'] = val
        self.save_job_settings()

    @property
    def draught(self):
        self.load_job_settings()
        return self._job_data['floater']['Draught']

    @draught.setter
    def draught(self, val):
        self._job_data['floater']['Draught'] = val
        self.save_job_settings()

    @property
    def fio(self):
        return self._fio

    @property
    def case_label(self):
        return self._case_label

    @case_label.setter
    def case_label(self, type=None):
        if type == None:
            self._case_label = None
        elif type == 'floater_data':
            nrc = self._job_data['floater']['Number of radial columns']
            rcd = self._job_data['floater']['Radial column diameter']
            gf = self._job_data['floater']['Gap factor']
            rh = self._job_data['floater']['Radial height']
            mt = self._job_data['floater']['Type']
            self._case_label = '{}{:0}C{:03.0f}-G{:02.0f}H{:03.0f}'.format(mt, nrc, rcd * 10, gf * 10, rh * 10)
        else:
            print('Cannot set case_label, \"{}\" is not a valid type\nUse \'floater_data\' or None (auto)'.format(type))
            exit()

        self.set_file_structure()

    @property
    def floater_data(self):
        self.load_job_settings()
        return self._job_data['floater']

    @floater_data.setter
    def floater_data(self, val):
        if type(val) is dict:
            for key in val:
                self._job_data['floater'][key] = val[key]
        else:
            print('Cannot set floater, {} is not dict'.format(val))
            exit()

    @property
    def load_cases(self):
        self.load_job_settings()
        return self._job_data['analysis']['simulations']['default']['load_cases']

    @load_cases.setter
    def load_cases(self, val):
        if type(val) is dict:
            for key in val:
                self._job_data['analysis']['simulations']['default']['load_cases'][key] = val[key]
        else:
            print('Cannot set analysis_data, {} is not dict'.format(val))
            exit()

    @property
    def job_data(self):
        return self._job_data
