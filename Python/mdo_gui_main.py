import math
from tkinter import *
from tkinter import filedialog
from pathlib import Path
from optimization.design_engine import Candidate, Parameter_Space
import os
import pygubu
from plot_forces import get_axs
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import pickle
import core.tool_box as tb
import json
from tkinter import messagebox
import core.environmental_conditions as ec
import core.common as cc
import numpy as np
from pyNemoh_root.pyNemoh import utility
import h5py

import matplotlib.pyplot as plt

# import module_locator
nax=np.newaxis
# my_path = module_locator.module_path()
# print(my_path)

try:
    DATA_DIR = os.path.abspath(os.path.dirname(__file__))
except NameError:
    DATA_DIR = os.path.dirname(os.path.abspath(sys.argv[0]))

DEG2RAD = 4 * math.atan(1) * 2 / 360

nsp = tb.NumericStringParser()


class StdoutRedirector(object):
    def __init__(self, text_widget):
        self.text_space = text_widget

    def write(self, string):
        self.text_space.insert('end', string)
        self.text_space.see('end')

    def flush(self):
        pass


class Application():

    def __init__(self):
        self.about_dialog = None
        self.nautic_zones_dialog = None
        self.fig = None
        self.builder = b = pygubu.Builder()
        b.add_from_file(os.path.join(DATA_DIR, 'main.ui'))
        b.add_resource_path(os.path.join(DATA_DIR, 'imgs'))

        self.mainwindow = b.get_object('mainwindow')
        self.mainwindow.iconbitmap('.\imgs\molo_256.ico')

        self.mainwindow.protocol("WM_DELETE_WINDOW", self.quit)

        b.connect_callbacks(self)

        self.case_cfg = None
        self.global_cfg = None
        self.default_case_cfg = None
        self.local_case_cfg = None
        self.this_candidate = None
        self.data_io_dir = None
        self.pkl_path = None
        self.fio = None

        self.nautic_zones_gif = PhotoImage(file=r".\imgs\nautic_zones.gif")


        self.template_dir = Path(os.getcwd()).joinpath('templates')

        with open(self.template_dir.joinpath('wtg_template.json'), 'r') as f:
            self._wtg_template = json.loads(f.read())

        self.cb1=self.builder.get_object('wtg_label_input')
        a = list()
        for key in self._wtg_template:
            a.append(key)
        self.cb1['value']=a
        #print(self.cb1['values'])
        self.cb1.current(0)

        self.get_output('main_text')

        self.load_global_cfg()
        print('Global config loaded')

        self.set_fio()
        self.load_case_cfg()
        print('Case config loaded')

        self.wtg_label_variable=self.builder.tkvariables.__getitem__('wtg_label_variable')
        self.wtg_label_variable.trace('w', self.react_to_wtg_label_change)
        self.react_to_wtg_label_change()

        self.plot_pressure_type_variable=self.builder.tkvariables.__getitem__('plot_pressure_type_variable')
        self.plot_pressure_type_variable.trace('w', self.react_to_pressure_type_change)
        self.plot_pressure_index_label_variable = self.builder.tkvariables.__getitem__('plot_pressure_index_label_variable')
        #self.react_to_pressure_type_change()

        # Some tweaks
        self.cb = self.builder.get_object('plot_pressure_type_input')
        self.cb.current(0)

    def set_fio(self, case_label_type='molo_model'):
        analyses_root_input = Path(self.builder.get_object('analyses_root_input').get())
        park_label_input = self.builder.get_object('park_label_input').get().replace(' ','_')
        wtg_label_input = self.builder.get_object('wtg_label_input').get().replace(' ','_')
        nrc = int(self.get_numeric('ncol_input'))
        rh = self.get_numeric('height_input')
        rcd = self.get_numeric('column_diameter_input')
        gf = self.get_numeric('gap_input')

        self.fio = cc.FileIOClass(analyses_root_input, park_label_input, wtg_label_input, case_label_type, nrc,
                                  rcd, gf, rh, templates_dir=None)

    def get_output(self, text_box):
        self.text_box = self.builder.get_object(text_box)
        self.text_box.delete('1.0', 'end')
        sys.stdout = StdoutRedirector(self.text_box)

    def show_about_dialog(self):
        if self.about_dialog is None:
            dialog = self.builder.get_object('dlg_about', self.mainwindow)
            self.about_dialog = dialog

            def dialog_btnclose_clicked():
                dialog.close()

            btnclose = self.builder.get_object('about_btnclose')
            btnclose['command'] = dialog_btnclose_clicked
            dialog.run()
        else:
            self.about_dialog.show()

    def show_nautic_zones_dialog(self):
        if self.nautic_zones_dialog is None:
            dialog = self.builder.get_object('dlg_nautic_zones', self.mainwindow)
            canvas = self.builder.get_object('dlg_nautic_zones_canvas')
            canvas.create_image(0, 0, anchor=NW, image=self.nautic_zones_gif)

            self.nautic_zones_dialog = dialog

            def dialog_btnclose_clicked():
                dialog.close()

            btnclose = self.builder.get_object('nautic_zones_btnclose')
            btnclose['command'] = dialog_btnclose_clicked
            dialog.run()
        else:
            self.nautic_zones_dialog.show()

    def get_config(self, cfg):
        for id in cfg['user_input']['object']:
            cfg['user_input']['object'][id]['value'] = self.builder.get_object(id).get()

    def save_cfg(self):
        # Save current settings
        self.get_config(self.global_cfg)
        self.get_config(self.case_cfg)

        f_global = Path(DATA_DIR).joinpath('mdo_global_cfg.json')
        with open(f_global, 'w') as f:
            f.write(json.dumps(self.global_cfg, indent=4, sort_keys=True))

        if isinstance(self.data_io_dir, Path):
            f_case = self.data_io_dir.joinpath('mdo_case_cfg.json')
            with open(str(f_case), 'w') as f:
                f.write(json.dumps(self.case_cfg, indent=4, sort_keys=True))

        f_case_default = Path(DATA_DIR).joinpath('mdo_default_case_cfg.json')
        with open(f_case_default, 'w') as f:
            f.write(json.dumps(self.case_cfg, indent=4, sort_keys=True))

    def load_cfg(self):
        self.load_global_cfg()
        self.load_case_cfg()

    def load_global_cfg(self):
        # Read global config
        f_global = Path(DATA_DIR).joinpath('mdo_global_cfg.json')

        with open(f_global, 'r') as f:
            self.global_cfg = json.loads(f.read())

        objects = self.global_cfg['user_input']['object']
        for id in objects:
            value = objects[id]['value']
            try:
                self.builder.get_object(id).delete(0, END)
                self.builder.get_object(id).insert(0, value)
            except:
                print('{} not defined'.format(id))

    def load_case_cfg(self):

        # Read default case config
        f_case = Path(DATA_DIR).joinpath('mdo_default_case_cfg.json')
        if f_case.exists():
            with open(f_case, 'r') as f:
                self.default_case_cfg = json.loads(f.read())

        # Read local case config
        if isinstance(self.data_io_dir, Path):
            f_case = self.data_io_dir.joinpath('mdo_case_cfg.json')
            if f_case.exists():
                with open(str(f_case), 'r') as f:
                    self.local_case_cfg = json.loads(f.read())

        if self.local_case_cfg == None:

            self.case_cfg = self.default_case_cfg

        else:
            self.case_cfg = self.local_case_cfg

            # Loop trough items in default case config and add missing items
            for key in self.default_case_cfg['user_input']['object']:
                if key not in self.case_cfg['user_input']['object']:
                    self.case_cfg[key] = self.default_case_cfg['user_input']['object'][key]

        objects = self.case_cfg['user_input']['object']
        for id in objects:
            value = objects[id]['value']
            try:
                self.builder.get_object(id).delete(0, END)
                self.builder.get_object(id).insert(0, value)
            except:
                print('{} not defined'.format(id))

    def quit(self, event=None):
        if messagebox.askokcancel("Quit", "Do you want to quit?"):
            self.save_cfg()

            self.mainwindow.quit()

    def btnplot_clicked(self):
        if not self.fig == None:
            self.fig.clear()
        print('Plotting ...')
        self.fig, self.axs = get_axs()
        print('Fig retrieved')
        self.fig.show()
        master = self.builder.get_object('main_canvas')
        canvas = FigureCanvasTkAgg(self.fig, master=master)
        canvas.draw()

    def plot(self):

        btnplot = self.builder.get_object('plot_button')
        btnplot['command'] = self.btnplot_clicked()

    def run(self):
        self.mainwindow.mainloop()

    def get_numeric(self, id):
        return nsp.eval(self.builder.get_object(id).get())

    def react_to_pressure_type_change(self, *args):
        str=self.plot_pressure_type_variable.get().replace('_',' ')
        if str == 'Radiation':
            self.plot_pressure_index_label_variable='Degree of freedom'
        else:
            self.plot_pressure_index_label_variable='Wave direction'

    def react_to_wtg_label_change(self, *args):
        str=self.wtg_label_variable.get().replace('_',' ')
        #print(str)
        self.builder.get_object('model_wtg_label')['text']=str
        self.set_fio()
        #print('Case directory is {:s}'.format())
        print('Case parent directory is {}'.format(self.fio.case_parent_dir))

    def create_model(self):
        self.get_output('model_text')
        self.set_fio()

        p = Parameter_Space(templates_dir=self.template_dir)
        p.ncol = int(self.get_numeric('ncol_input'))
        print('Number of columns: {:6.1f}'.format(p.ncol))

        p.height = self.get_numeric('height_input')
        print('Column height: {:6.3f}'.format(p.height))

        p.radial_column_diameter = self.get_numeric('column_diameter_input')
        print('Column diameter: {:6.3f}'.format(p.radial_column_diameter))
        p.radial_column_thickness = self.get_numeric('column_thickness_input')

        p.gap = self.get_numeric('gap_input')
        print('Gap factor: {:6.3f}'.format(p.gap))



        p.upper_flange_thickness = self.get_numeric('upper_flange_thickness_input')
        p.lower_flange_thickness = self.get_numeric('lower_flange_thickness_input')


        p.lower_plate_width= p.radial_column_diameter + 2 * self.get_numeric('lower_flange_overwidth_input')
        p.lower_flange_overlength = self.get_numeric('lower_flange_overlength_input')



        p.filling_ratio = [0.1, 0.1, 0.1]
        self.this_candidate = Candidate(p,self.fio)
        self.this_candidate.settings.lower_face_corrected_z_pos = False
        self.this_candidate.settings.wtg_model = self.wtg_label_variable.get()
        print('\nTurbine type is the {:s}'.format(self.wtg_label_variable.get()))

        self.this_candidate.init_model(state='New')
        self.data_io_dir = self.this_candidate.settings.fio.data_io_dir
        self.pkl_path = self.data_io_dir.joinpath('this_candidate.pkl')
        with open(self.pkl_path, "wb") as f:
            pickle.dump(self.this_candidate, f)

    def calc_stability(self):

        self.get_output('stability_text')

        with open(self.pkl_path, "rb") as f:
            self.this_candidate = pickle.load(f)

        self.this_candidate.settings.create_stability_movie = self.builder.get_variable('stability_movie_chkbtn_var')
        if self.this_candidate.settings.create_stability_movie:
            print('A movie will be made', flush=True)

        r = self.this_candidate.intact_stability_ratio()
        if r >= 1.4:
            is_stable = True
            print('This candidate is stable with r = {:1.0f}%'.format(r * 100))
        else:
            is_stable = False
            print('This candidate is NOT stable with r = {:1.0f}%'.format(r * 100))

    def run_hydro(self):
        self.get_output('hydro_text')

        with open(self.pkl_path, "rb") as f:
            self.this_candidate = pickle.load(f)

        self.this_candidate.settings.load_cases = {
            "num_wave_frequencies": int(self.get_numeric('num_wave_frequencies_input')),
            "min_wave_frequencies": self.get_numeric('min_wave_frequencies_input'),
            "max_wave_frequencies": self.get_numeric('max_wave_frequencies_input'),
            "num_wave_directions": int(self.get_numeric('num_wave_directions_input')),
            "min_wave_directions": self.get_numeric('min_wave_directions_input'),
            "max_wave_directions": self.get_numeric('max_wave_directions_input'),
        }

        self.this_candidate.hydrodynamic_analysis()

    def open_workspace(self, title=None, dirName=None):
        options = {}
        options['initialdir'] = dirName
        options['title'] = title
        options['mustexist'] = False
        fileName = filedialog.askdirectory(**options)
        if fileName == "":
            return None
        else:
            return fileName

    def calc_response(self):
        self.get_output('response_text')

        area = int(self.builder.get_object('calc_response_area_input').get())
        ibeta = int(self.builder.get_object('calc_response_ibeta_input').get())
        yr = float(self.builder.get_object(
            'calc_response_yr_input').get())  # Return period in years, statistics conditioned for 3hr storms

        tc = self.this_candidate

        force_label = ['Fx [MN]', 'Fy [MN]', 'Fz [MN]', 'Mx [MNm]', 'My [MNm]', 'Mz [MNm]']

        # Create list of short terms from contour line
        # area = this_candidate.settings.park_data['Design Basis']['Area']

        tc.ltwc1 = ec.Long_Term_Wave_Conditions(area=area)
        tc.cl = tc.ltwc1.contour_line(27)
        tc.contourline = []
        # this_candidate.cl = np.asarray([[3.5, 13.5, 1],[9, 9.5, 5],[8, 11, 3.6],[7, 11, 2.6],[7, 11.5, 2.12]])

        for hs, tz, in tc.cl:
            # print(' {:6.1f} {:6.1f}'.format(hs, tz))
            tc.contourline.append(ec.Short_Term_Wave_Conditions(hs=hs, tz=tz, gamma=None))

        tc.settings.radiaton_damping_factor = 1

        bool_lin = self.builder.tkvariables.__getitem__('lin_visc_damp_var').get()
        tc.settings.do_linearize = bool_lin

        if tc.settings.do_linearize:
            print('\nLinearized viscous damping will be included')
        else:
            print('\nLinearized viscous damping will NOT be included')

        print('\nCase:\t{}'.format(self.this_candidate.settings.case_label))

        print('\nEigenvalue sollution WITH added mass')

        tb.eigenvalprint(tc.loads.m + tc.loads.ma[0, :, :], tc.loads.k)
        np.set_printoptions(precision=4)
        print(tc.loads.m)
        print(tc.loads.ma[0, :, :])
        print(tc.loads.k)


        d_col_central = tc.settings.floater_data['Central column diameter']
        d_col_radial = tc.settings.floater_data['Radial']['Column']['Diameter']
        n_col_radial = tc.settings.floater_data['Radial']['Number of columns']
        gap = tc.settings.floater_data['Gap factor']
        height = tc.settings.floater_data['Radial']['Heigth']
        draught = tc.hs_floater.hs_data['draught']

        # Get section forces
        sp_x = d_col_central * (0.5)  # + 0.8 + 1 + 0.8 + 1)
        sp_z = height / 2 - draught
        # sp_x = -100
        sp_z = 0
        sp1 = [sp_x, 0, sp_z]  # Used for moment reference
        sn1 = [1, 0, 0]

        radial_extreme = d_col_central * (0.5) + (1 + gap) * n_col_radial * d_col_radial

        print('\nAirgap point at {:1.2f}'.format(radial_extreme))

        c1 = np.asarray([[radial_extreme, 0, 0]])

        with h5py.File(tc.settings.fio.nemoh_root.joinpath('db.hdf5'), "r") as hdf5_db:
            environment = utility.read_environment(hdf5_db)

        k_wave = np.zeros(tc.loads.nw, dtype=float)
        for iw, val in enumerate(tc.loads.w):
            k_wave[iw] = utility.compute_wave_number(val, environment)
        w_bar = (radial_extreme - environment.x_eff) * np.cos(tc.loads.beta) + (0 - environment.y_eff) * np.sin(
            tc.loads.beta)
        eta = np.exp(utility.II * k_wave[:,nax] * w_bar[nax,:])

        print('\n {:^6s} {:^6s} {:^6s}  {:^6s}  {:^6s}  {:^6s}'.format('Hs', 'Tp', 'Gamma', force_label[2],
                                                                       force_label[4], 'AG'))
        for stwc in tc.contourline:
            tc.init_response(short_term_wave_condition=stwc)

            imass, ipanel, istrip = tc.response.get_section_index(sp1, sn1)
            f_sec1 = tc.response.assemble_forces(imass, ipanel, istrip, moment_ref_point=[0, 0, 0])
            del imass, ipanel

            fz = f_sec1['Dynamic']['SUM'][:, ibeta, 2] / 1000000
            my = f_sec1['Dynamic']['SUM'][:, ibeta, 4] / 1000000

            fz_elm = stwc.expected_largest_maximum(fz, tc.loads.w)
            my_elm = stwc.expected_largest_maximum(my, tc.loads.w)

            p1_rao = tc.response.point_rao(c1)[:, ibeta, 0, 2].flatten()
            ag_rao = eta[:,ibeta] + p1_rao
            ag1 = stwc.expected_largest_maximum(ag_rao, tc.loads.w)

            print(' {:6.1f} {:6.1f} {:6.1f}  {:6.2f}  {:6.1f}  {:6.1f} '.format(stwc.hs, stwc.tp, stwc.gamma, fz_elm,
                                                                                my_elm, ag1))

    def plot_pressure(self):
        self.get_output('plot_text')
        ifreq = int(self.builder.get_object('plot_pressure_ifreq_input').get())
        pressure_index = int(self.builder.get_object('plot_pressure_index_input').get())
        pressure_type = self.builder.get_object('plot_pressure_type_input').get()
        axis = int(self.builder.get_object('plot_pressure_axis_input').get())
        self.this_candidate.loads.show_pressure(ifreq, pressure_index, axis, pressure_type)

    def plot_rao(self):
        self.get_output('plot_text')
        ifreq = int(self.builder.get_object('plot_rao_ifreq_input').get())


        print('\nCase:\t{}'.format(self.this_candidate.settings.case_label))

        yr = self.this_candidate.settings.park_data['Design Basis']['ULS']['Return period']

        # Create list of short terms from contour line
        area = self.this_candidate.settings.park_data['Design Basis']['Area']
        self.this_candidate.ltwc1 = ec.Long_Term_Wave_Conditions(area=area)
        self.this_candidate.cl = self.this_candidate.ltwc1.contour_line(yr)
        self.this_candidate.contourline = []
        # this_candidate.cl = np.asarray([[3.5, 13.5, 1],[9, 9.5, 5],[8, 11, 3.6],[7, 11, 2.6],[7, 11.5, 2.12]])

        for hs, tz in self.this_candidate.cl:
            self.this_candidate.contourline.append(ec.Short_Term_Wave_Conditions(hs=hs, tz=tz))

        self.this_candidate.settings.radiaton_damping_factor = 1

        self.this_candidate.settings.do_linearize = True
        self.this_candidate.init_load()
        stwc = ec.Short_Term_Wave_Conditions(hs=9.5, tz=7.3)
        self.this_candidate.init_response(short_term_wave_condition=stwc)

        # Get section forces
        sp_x = self.this_candidate.settings.floater_data['Central column diameter'] * (0.5)  # + 0.8 + 1 + 0.8 + 1)
        sp_z = self.this_candidate.settings.floater_data['Radial']['Heigth'] / 2 - self.this_candidate.hs_floater.hs_data[
            'draught']
        #sp_x = -100
        #sp_z = 0
        section_point = [sp_x, 0, sp_z]  # Used for moment reference
        section_normal = [1, 0, 0]

        imass, ipanel, idamp = self.this_candidate.response.get_section_index(section_point, section_normal)
        self.this_candidate.f_sec1 = self.this_candidate.response.assemble_forces(imass, ipanel, idamp,
                                                                        moment_ref_point=[0, 0, 0])
        del imass, ipanel

        ibeta = 0
        idof = 2

        stat_force = self.this_candidate.f_sec1['Static']

        factor = 1 / 1000000

        print('\n--------------------------\n S T A T I C   F O R C E\n--------------------------')
        for key in stat_force:
            a = np.abs(stat_force[key][2]) * factor
            b = np.angle(stat_force[key][2])
            print('{:20} {:5.2f} {: 5.2f}'.format(key, a, b))

        w = self.this_candidate.loads.w
        fig, axs = plt.subplots(2, 2)
        dyn_force = self.this_candidate.f_sec1['Dynamic']

        x_tics = 1 / (2 * np.pi / w)
        x_label = 'Frequency [Hz]'

        # --------------
        # FREQUENCY
        # --------------
        ifreq_print = ifreq

        print('\n--------------------------\n D Y N A M I C   F O R C E\n--------------------------')
        print('{:16} {:5.3f}\n'.format(x_label, x_tics[ifreq_print]))
        all_keys = ['Froude-Krylof', 'Diffraction', 'Mass', 'Added mass', 'Radiation damping', 'Viscous damping',
                    'Buoyancy', 'SUM']
        plot_keys = list(all_keys[i] for i in [0, 1, 2, 3, 4, 5, 6, 7])
        for key in dyn_force:
            if key in plot_keys:
                if key == 'SUM':
                    linewidth = 2
                else:
                    linewidth = 1
                abs_val = np.abs(dyn_force[key][:, ibeta, idof]) * factor

                phase_val = np.angle(dyn_force[key][:, ibeta, idof])
                axs[0, 0].plot(x_tics, abs_val, label=key, linewidth=linewidth)
                axs[0, 1].plot(x_tics, phase_val, label=key)
                # print(abs_val[:])
                print('{:20} {:6.2f} {: 5.2f}'.format(key, abs_val[ifreq_print], phase_val[ifreq_print]))

        axs[0, 0].set_title('Amplitude')
        force_label = ['Fx [MN]', 'Fy [MN]', 'Fz [MN]', 'Mx [MNm]', 'My [MNm]', 'Mz [MNm]']
        axs[0, 0].set_ylabel(force_label[idof])
        # axs[0, 0].set_yscale('log')
        # axs[0, 0].set_ylim([0, 1])

        axs[0, 0].legend()
        axs[0, 0].grid()
        axs[0, 1].set_title('Phase')
        axs[0, 1].legend()
        axs[0, 1].grid()
        axs[1, 0].set_ylabel('RAO')
        axs[1, 0].set_xlabel(x_label)
        axs[1, 0].grid()
        axs[1, 1].set_xlabel(x_label)
        axs[1, 1].grid()
        print('\n--------------------------\n R A O\n--------------------------')
        d = {'Heave': 2, 'Pitch': 4}
        ax1 = axs[1, 0]
        ax2 = ax1.twinx()

        for r in [self.this_candidate.response.rao, self.this_candidate.response.rao_init]:

            for key in d:

                if key == 'Heave':
                    abs_val_rao = np.abs(r[:, ibeta, d[key]])
                    lns1 = ax1.plot(x_tics, abs_val_rao, label=key)
                    phase_val_rao = np.angle(r[:, ibeta, d[key]])
                    print('{:20} {:6.3f} {: 7.4f}'.format(key, abs_val_rao[ifreq_print],
                                                          phase_val_rao[ifreq_print]))


                else:
                    abs_val_rao = np.abs(r[:, ibeta, d[key]])
                    lns2 = ax2.plot(x_tics, abs_val_rao, '-r', label=key)

                    phase_val_rao = np.angle(r[:, ibeta, d[key]])
                    print('{:20} {:6.3f} {: 7.4f} ({: 5.1f} deg)'.format(key, abs_val_rao[ifreq_print],
                                                                         phase_val_rao[ifreq_print],
                                                                         abs_val_rao[ifreq_print] * 180 / np.pi))

                axs[1, 1].plot(x_tics, phase_val_rao, label=key)

        lns = lns1 + lns2
        labs = [l.get_label() for l in lns]
        axs[1, 0].legend(lns, labs, loc=0)
        axs[1, 1].legend()
        plt.suptitle(
            '{}\nSection point: [{sp[0]:1.2f}, {sp[1]:1.2f}, {sp[2]:1.2f}]\nSection normal: [{sn[0]:1.2f}, {sn[1]:1.2f}, {sn[2]:1.2f}]'.format(
                self.this_candidate.settings.case_label, sp=section_point, sn=section_normal))
        plt.show()

        # f_rad=np.sum(this_candidate.loads._force['Radiation'][ifreq_print,idof,:,2])
        # A=np.imag(f_rad)/w[ifreq_print]
        # B=-np.real(f_rad)
        # print('\n{:15} {:5.0f} {:5.0f}'.format('A and B', A, B))
        # z = A*np.exp(complex(0, 1)*B)
        # print('\n{:15} {:5.0f}'.format('Z', z))

        # a=this_candidate.response._point_mass[:, 2, 2]
        # b=this_candidate.response._point_mass_centers
        #
        # print()
        # print(sum(a)
        # print(this_candidate.loads._m[2,2])
        #
        # I=np.sum(a[:,np.newaxis]*b**2,axis=0)
        # print(I)
        # print(np.diag(this_candidate.loads._m[3:,3:]))
        np.set_printoptions(precision=3)
        #
        # print()
        # print(this_candidate.loads._m)
        # print()
        # print(this_candidate.loads._ma[ifreq_print,:,:])
        # print()
        # print(this_candidate.loads._k)
        # print()
        # print(np.diag(this_candidate.response._c_visc[ifreq_print,:]))

    def gui_init_load(self):
        if not self.this_candidate == None:
            if self.this_candidate.loads == None:
                self.this_candidate.init_load()
                print('Loads initialized')

        else:
            print('Please initialize Model and/or run Hydro')

    def gui_init_plot(self, event=None):
        self.get_output('plot_text')
        self.gui_init_load()

    def gui_init_response(self, event=None):
        self.get_output('response_text')
        self.gui_init_load()


if __name__ == '__main__':
    # root = tk.Tk()
    app = Application()
    app.run()
