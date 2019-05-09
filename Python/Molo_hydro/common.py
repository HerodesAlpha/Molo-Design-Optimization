from pathlib import Path


class PhysicalQuantities():
    # This is the only place allowed to put these quantities
    def __init__(self):
        self._rho_sw = 1025
        self._grav = 9.81


class FileIO(object):
    def __init__(self, root_dir, templates_dir, case_number=None):
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

        self._templates_dir = templates_dir
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
