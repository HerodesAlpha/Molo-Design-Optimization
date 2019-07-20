__author__ = "Eivind Sonju"
__copyright__ = "Copyright (C) 2017-2019 Verbun AS. All rights reserved."
__version__ = "2.0"

import json
import os

# import matplotlib.pyplot as plt

from pathlib import Path
import getpass
import sys
from pylatex import Document, PageStyle, Head, Foot, MiniPage, \
    StandAloneGraphic, MultiColumn, Tabu, LongTabu, LargeText, MediumText, \
    LineBreak, NewPage, Tabularx, TextColor, simple_page_number, Command
from pylatex.utils import bold, NoEscape


class PhysicalQuantities():
    # This is the only place allowed to put physical quantities
    def __init__(self):
        self._rho_sw = 1025
        self._gravity = -9.81  # In global coordinate system
        self._rho_st = 7850

    @property
    def rho_sw(self):
        return self._rho_sw

    @property
    def gravity(self):
        return self._gravity

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

        self._mesh_name = None

        self._case_dir = case_parent_dir.joinpath(self._case_label)

        self._nemoh_root = self._case_dir.joinpath('nemoh')
        self._gmsh_root = self._case_dir.joinpath('gmsh')
        self._data_io_dir = self._case_dir.joinpath('data_io')
        self._nemoh_results_dir = self._nemoh_root.joinpath('results')
        self._nemoh_mesh_dir = self._nemoh_root.joinpath('mesh')
        self._stability_dir = self._case_dir.joinpath('stability')

        self._templates_dir = Path(os.getcwd()).joinpath('templates')

        self._gmsh_exe = r'C:\Users\{}\OneDrive - Verbun AS\Divisions\Offshore Wind\Library\Software\Bin\gmsh-4.2.2-Windows64\gmsh.exe'.format(
                getpass.getuser())
        assert (Path(self._gmsh_exe).exists())

        self._freecad_path = r'C:\Program Files\FreeCAD 0.18\bin'
        assert (Path(self._freecad_path).exists())
        sys.path.append(self._freecad_path)

        # else:
        #     print('{} does not exits'.format(path_to_gmsh_exe))
        #     exit()

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

        self._json_list = ['park', 'rna', 'tower', 'floater', 'analysis', 'design_basis']

        self._mesh_name = None
        self._do_equilibrate = True

        self._job_data = dict()

        self._analyses_root = analyses_root
        self._park_label = park_label
        self._wtg_label = wtg_label
        self._case_label = None

        self._fio = None
        self._report = None

        # Collect template data
        for item in self._json_list:
            with open(Path(os.getcwd()).joinpath('templates').joinpath('{}_template.json'.format(item)), 'r') as f:
                self._job_data[item] = json.loads(f.read())

        # Save updated template to template dir
        for item in self._json_list:
            with open(Path(os.getcwd()).joinpath('templates').joinpath('{}_template.json'.format(item)), 'w') as f:
                f.write(json.dumps(self._job_data[item], indent=4, sort_keys=True))

        self._thin_panels = self._job_data['analysis']['simulations']['default']['calculation']['thin_panels']
        self._use_dipols = self._job_data['analysis']['simulations']['default']['calculation'][
            'use_dipoles_implementation']

        self._thin_panel_offset = self._job_data['floater']['Thin panel offset']
        self._flange_thickness = self._job_data['floater']['Radial']['Flange']['Lower']['Plate']['Thickness']
        self._create_model = False
        self._calc_gz = False
        self._run_nemoh = False
        self._postprocessing = False
        self.set_molo_label()

    def set_molo_label(self):
        nrc = self._job_data['floater']['Radial']['Number of columns']
        rcd = self._job_data['floater']['Radial']['Column']['Diameter']
        gf = self._job_data['floater']['Gap factor']
        rh = self._job_data['floater']['Radial']['Heigth']
        mt = self._job_data['floater']['Type']
        self._molo_label = '{}-{:0}C{:03.0f}-G{:02.0f}H{:03.0f}'.format(mt, nrc, rcd * 10, gf * 10, rh * 10)

    def set_file_structure(self):
        self._fio = FileIOClass(self._analyses_root, self._park_label, self._wtg_label, self._case_label)
        # print(self._job_data['analysis']['simulations'])
        self._job_data['analysis']['simulations']['sim01']['simulation_dir'] = str(self._fio.nemoh_root)
        self.save_job_settings()

        print(str(self._fio.case_dir))
        self.init_report()

        # self._report = Document(self._fio._case_dir.joinpath('report_{}'.format(self._molo_label)))
        # self._report.preamble.append(Command('title', '{}'.format(self._molo_label)))
        # self._report.preamble.append(Command('author', NoEscape(r'Eivind S{\o}nju')))
        # self._report.preamble.append(Command('date', NoEscape(r'\today')))
        # self._report.append(NoEscape(r'\maketitle'))

    def init_report(self):
        geometry_options = {
                "head"           : "40pt",
                "margin"         : "0.5in",
                "bottom"         : "0.6in",
                "includeheadfoot": True
        }
        self._report = Document(self._fio._case_dir.joinpath('report_{}'.format(self._molo_label)),
                                geometry_options=geometry_options)

        # Generating first page style
        first_page = PageStyle("firstpage")

        # Header image

        with first_page.create(Head("L")) as header_left:
            with header_left.create(MiniPage(width=NoEscape(r"0.49\textwidth"),
                                             pos='c')) as logo_wrapper:
                logo_file = str(Path(os.getcwd()).joinpath('templates').joinpath('logo.png'))
                #logo_file= '{' + logo_file + '}'
                logo_file =  logo_file.replace('\\','/')
                # print(str(logo_file))
                logo_wrapper.append(StandAloneGraphic(image_options="width=120px",
                                                      filename=logo_file))
                                                      # filename='logo.png'))

        # Add document title
        with first_page.create(Head("R")) as right_header:
            with right_header.create(MiniPage(width=NoEscape(r"0.49\textwidth"),
                                              pos='c', align='r')) as title_wrapper:
                title_wrapper.append(LargeText(bold('{}'.format(self._molo_label))))
                title_wrapper.append(LineBreak())
                title_wrapper.append(MediumText(bold(NoEscape(r'\today'))))

        # Add footer
        with first_page.create(Foot("C")) as footer:
            message = "Important message please read"
            with footer.create(Tabularx(
                    "X X X X",
                    width_argument=NoEscape(r"\textwidth"))) as footer_table:
                footer_table.add_row(
                        [MultiColumn(4, align='l', data=TextColor("blue", message))])
                footer_table.add_hline(color="blue")
                footer_table.add_empty_row()

                branch_address = MiniPage(
                        width=NoEscape(r"0.25\textwidth"),
                        pos='t')
                branch_address.append("Radyrveien 32")
                branch_address.append("\n")
                branch_address.append("1555 SON")

                document_details = MiniPage(width=NoEscape(r"0.25\textwidth"),
                                            pos='t', align='r')
                document_details.append("1000")
                document_details.append(LineBreak())
                document_details.append(simple_page_number())

                footer_table.add_row([branch_address, branch_address,
                                      branch_address, document_details])

        self._report.preamble.append(first_page)
        # End first page style

        # # Add customer information
        # with self._report.create(Tabu("X[l] X[r]")) as first_page_table:
        #     customer = MiniPage(width=NoEscape(r"0.49\textwidth"), pos='h')
        #     customer.append("Verna Volcano")
        #     customer.append("\n")
        #     customer.append("For some Person")
        #     customer.append("\n")
        #     customer.append("Address1")
        #     customer.append("\n")
        #     customer.append("Address2")
        #     customer.append("\n")
        #     customer.append("Address3")
        #
        #     # Add branch information
        #     branch = MiniPage(width=NoEscape(r"0.49\textwidth"), pos='t!',
        #                       align='r')
        #     branch.append("Branch no.")
        #     branch.append(LineBreak())
        #     branch.append(bold("1181..."))
        #     branch.append(LineBreak())
        #     branch.append(bold("TIB Cheque"))
        #
        #     first_page_table.add_row([customer, branch])
        #     first_page_table.add_empty_row()
        #

        self._report.add_color(name="lightgray", model="gray", description="0.80")
        #
        # # Add statement table
        # with self._report.create(LongTabu("X[l] X[2l] X[r] X[r] X[r]",
        #                                   row_height=1.5)) as data_table:
        #     data_table.add_row(["date",
        #                         "description",
        #                         "debits($)",
        #                         "credits($)",
        #                         "balance($)"],
        #                        mapper=bold,
        #                        color="lightgray")
        #     data_table.add_empty_row()
        #     data_table.add_hline()
        #     row = ["2016-JUN-01", "Test", "$100", "$1000", "-$900"]
        #     for i in range(30):
        #         if (i % 2) == 0:
        #             data_table.add_row(row, color="lightgray")
        #         else:
        #             data_table.add_row(row)

        self._report.preamble.append(Command('title', 'MOLO Conceptual Design Report'))
        self._report.change_document_style("firstpage")
        self._report.append(NewPage())


        # # Add cheque images
        # with self._report.create(LongTabu("X[c] X[c]")) as cheque_table:
        #     cheque_file = os.path.join(os.path.dirname(__file__),
        #                                'chequeexample.png')
        #     cheque = StandAloneGraphic(cheque_file, image_options="width=200px")
        #     for i in range(0, 20):
        #         cheque_table.add_row([cheque, cheque])

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

    def remove_old_db(self):
        db_file = self._fio.nemoh_root.joinpath('db.hdf5')
        if db_file.is_file():
            db_file.unlink()
            print('\ndb.hdf5 deleted from {}\n'.format(str(self._fio.nemoh_root)))

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
            self._case_label = None  # Auto numbering
        elif type == 'floater_data':
            self._case_label = self._molo_label
        else:
            self._case_label = type
        print('\nMODEL NAME: {}'.format(self._molo_label))

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

    @property
    def mesh_name(self):
        return self._mesh_name

    @mesh_name.setter
    def mesh_name(self, val):
        self._mesh_name = val

    @property
    def create_model(self):
        return self._create_model

    @create_model.setter
    def create_model(self, val):
        self._create_model = val

    @property
    def calc_intact_stability(self):
        return self._calc_gz

    @calc_intact_stability.setter
    def calc_intact_stability(self, val):
        self._calc_gz = val

    @property
    def run_nemoh(self):
        return self._run_nemoh

    @run_nemoh.setter
    def run_nemoh(self, val):
        self._run_nemoh = val

    @property
    def postprocessing(self):
        return self._postprocessing

    @postprocessing.setter
    def postprocessing(self, val):
        self._postprocessing = val

    @property
    def do_equilibrate(self):
        return self._do_equilibrate

    @do_equilibrate.setter
    def do_equilibrate(self, val):
        self._do_equilibrate = val

    @property
    def use_symmmetri(self):
        return self._use_symmmetri

    @use_symmmetri.setter
    def use_symmmetri(self, val):
        self._use_symmmetri = val

    @property
    def thin_panels(self):
        return self._thin_panels

    @thin_panels.setter
    def thin_panels(self, val):
        self._thin_panels = val
        self._job_data['analysis']['simulations']['default']['calculation']['thin_panels'] = self._thin_panels
        self.save_job_settings()

    @property
    def use_dipols(self):
        return self._use_dipols

    @use_dipols.setter
    def use_dipols(self, val):
        self._use_dipols = val
        self._job_data['analysis']['simulations']['default']['calculation'][
            'use_dipoles_implementation'] = self._use_dipols
        self.save_job_settings()

    @property
    def thin_panel_offset(self):
        return self._thin_panel_offset

    @thin_panel_offset.setter
    def thin_panel_offset(self, val):
        self._thin_panel_offset = val
        self._job_data['floater']['Thin panel offset'] = self._thin_panel_offset
        self.save_job_settings()

    @property
    def flange_thickness(self):
        return self._flange_thickness

    @flange_thickness.setter
    def flange_thickness(self, val):
        self._flange_thickness = val
        self._job_data['floater']['Lower flange thickness'] = self._flange_thickness
        self.save_job_settings()
