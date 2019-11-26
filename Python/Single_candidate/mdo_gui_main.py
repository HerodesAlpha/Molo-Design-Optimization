# file: myapp.py
import sys
import os
import random
import math
from tkinter import *
from tkinter import filedialog
import tkinter as tk
from pathlib import Path
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

        self.load_cfg()
        self.this_candidate=None

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
        for id in self.cfg['user_input']['object']:
            value = self.builder.get_object(id).get()
            self.cfg['user_input']['object'][id]['value'] = value

        with open('mdo_cfg.json', 'w') as f:
            f.write(json.dumps(self.cfg, indent=4, sort_keys=True))

    def load_cfg(self):
        # Read config
        with open('mdo_cfg.json', 'r') as f:
            self.cfg = json.loads(f.read())

        objects=self.cfg['user_input']['object']
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

        analyses_root = Path(self.builder.get_object('analyses_root_input').get())
        print(analyses_root)

        park_label = self.builder.get_object('park_label_input').get()
        wtg_label = self.builder.get_object('wtg_label_input').get()
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

        p.filling_ratio = [0.1, 0.1]
        this_candidate = Candidate(analyses_root, park_label, wtg_label, p, case_label_type='molo_model')
        this_candidate.settings.wtg_model = "Vestas 9.5"
        this_candidate.init_model(state='New')
        with open("this_candidate.pkl", "wb") as f:
            pickle.dump(this_candidate, f)

    def calc_stability(self):

        self.get_output('stability_text')

        with open("this_candidate.pkl", "rb") as f:
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

        with open("this_candidate.pkl", "rb") as f:
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


if __name__ == '__main__':
    # root = tk.Tk()
    app = Application()
    app.run()
