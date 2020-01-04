# file: myapp.py
import sys
import os
import random
import math
from tkinter import *
from tkinter import filedialog
import tkinter as tk
from pathlib import Path, PurePath
import pickle
from design_engine import Candidate, Parameter_Space
import os
import io
from contextlib import redirect_stdout
import pygubu
from plot_forces import get_axs
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import pickle
import tool_box as tb
import json
from tkinter import messagebox
import numpy as np
import environmental_conditions as ec
import module_locator
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
        self.this_candidate=None
        self.data_io_dir=None
        self.pkl_path=None

        self.load_cfg()

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

    def save_cfg(self):
        # Save current settings
        for id in self.global_cfg['user_input']['object']:
            value = self.builder.get_object(id).get()
            self.global_cfg['user_input']['object'][id]['value'] = value

        f_global = Path(DATA_DIR).joinpath('mdo_global_cfg.json')
        with open(f_global, 'w') as f:
            f.write(json.dumps(self.global_cfg, indent=4, sort_keys=True))

        f_case = self.data_io_dir.joinpath('mdo_case_cfg.json')
        with open(f_case, 'w') as f:
            f.write(json.dumps(self.case_cfg, indent=4, sort_keys=True))

    def load_cfg(self):
        # Read global config
        f_global=Path(DATA_DIR).joinpath('mdo_global_cfg.json')

        with open(f_global, 'r') as f:
            self.global_cfg = json.loads(f.read())

        if  isinstance(self.data_io_dir, Path):
            f_case = self.data_io_dir.joinpath('mdo_case_cfg.json')
            if f_case.exists():
                with open(str(f_case), 'r') as f:
                    self.case_cfg = json.loads(f.read())
        else:
            f_case = Path(DATA_DIR).joinpath('mdo_default_case_cfg.json')
            if f_case.exists():
                with open(f_case, 'r') as f:
                    self.case_cfg = json.loads(f.read())

        objects=self.global_cfg['user_input']['object']
        objects.update(self.case_cfg['user_input']['object'])
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

    def create_model(self):
        self.get_output('model_text')

        analyses_root_input = Path(self.builder.get_object('analyses_root_input').get())
        print(analyses_root_input)

        park_label_input = self.builder.get_object('park_label_input').get()
        wtg_label_input = self.builder.get_object('wtg_label_input').get()
        template_dir = Path(os.getcwd()).parents[0].joinpath('Optimization').joinpath('templates')

        p = Parameter_Space(templates_dir=template_dir)
        p.ncol = int(self.get_numeric('ncol_input'))
        print('Number of columns: {:6.1f}'.format(p.ncol))

        p.height = self.get_numeric('height_input')
        print('Column height: {:6.3f}'.format(p.height))

        p.column_diameter = self.get_numeric('column_diameter_input')
        print('Column diameter: {:6.3f}'.format(p.column_diameter))

        p.gap = self.get_numeric('gap_input')
        print('Gap factor: {:6.3f}'.format(p.gap))

        p.filling_ratio = [0.1, 0.1, 0.1]
        self.this_candidate = Candidate(analyses_root_input, park_label_input, wtg_label_input, p, case_label_type='molo_model')
        self.this_candidate.settings.wtg_model = "Vestas 9.5"
        self.this_candidate.init_model(state='New')
        self.data_io_dir=self.this_candidate.settings.fio.data_io_dir
        self.pkl_path=self.data_io_dir.joinpath('this_candidate.pkl')
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
                "num_wave_directions" : int(self.get_numeric('num_wave_directions_input')),
                "min_wave_directions" : self.get_numeric('min_wave_directions_input'),
                "max_wave_directions" : self.get_numeric('max_wave_directions_input'),
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
        tc = self.this_candidate

        print('\nCase:\t{}'.format(self.this_candidate.settings.case_label))

        force_label = ['Fx [MN]', 'Fy [MN]', 'Fz [MN]', 'Mx [MNm]', 'My [MNm]', 'Mz [MNm]']

        yr = self.this_candidate.settings.park_data['Design Basis']['ULS']['Return period']

        # Create list of short terms from contour line
        #area = this_candidate.settings.park_data['Design Basis']['Area']
        area = 9
        tc.ltwc1 = ec.Long_Term_Wave_Conditions(area=area)
        tc.cl = tc.ltwc1.contour_line(27)
        tc.contourline = []
        #this_candidate.cl = np.asarray([[3.5, 13.5, 1],[9, 9.5, 5],[8, 11, 3.6],[7, 11, 2.6],[7, 11.5, 2.12]])

        for hs, tz, in tc.cl:
            #print(' {:6.1f} {:6.1f}'.format(hs, tz))
            tc.contourline.append(ec.Short_Term_Wave_Conditions(hs=hs, tz=tz, gamma=None))

        tc.settings.radiaton_damping_factor = 1



        tc.settings.do_linearize=True
        tc.init_load()
        ibeta = 0
        for stwc in tc.contourline:
            tc.init_response(short_term_wave_condition=stwc)

            # Get section forces
            sp_x = tc.settings.floater_data['Central column diameter'] * (0.5)  # + 0.8 + 1 + 0.8 + 1)
            sp_z = tc.settings.floater_data['Radial']['Heigth'] / 2 - tc.hs_floater.hs_data[
                'draught']
            # sp_x = -100
            sp_z = 0
            section_point = [sp_x, 0, sp_z]  # Used for moment reference
            section_normal = [1, 0, 0]

            imass, ipanel, istrip = tc.response.get_section_index(section_point, section_normal)
            tc.f_sec1 = tc.response.assemble_forces(imass, ipanel, istrip, moment_ref_point=[0, 0, 0])
            del imass, ipanel




            fz = tc.f_sec1['Dynamic']['SUM'][:, ibeta, 2] / 1000000
            my = tc.f_sec1['Dynamic']['SUM'][:, ibeta, 4] / 1000000



            fz_elm=stwc.expected_largest_maximum(fz, tc.loads.w)
            my_elm=stwc.expected_largest_maximum(my, tc.loads.w)


            print(' {:6.1f} {:6.1f} {:6.1f}  {:6.2f}  {:6.1f} '.format(stwc.hs, stwc.tp, stwc.gamma, fz_elm, my_elm))

if __name__ == '__main__':
    # root = tk.Tk()
    app = Application()
    app.run()
