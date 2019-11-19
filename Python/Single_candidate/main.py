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

try:
    DATA_DIR = os.path.abspath(os.path.dirname(__file__))
except NameError:
    DATA_DIR = os.path.dirname(os.path.abspath(sys.argv[0]))

DEG2RAD = 4 * math.atan(1) * 2 / 360

class StdoutRedirector(object):
    def __init__(self,text_widget):
        self.text_space = text_widget

    def write(self,string):
        self.text_space.insert('end', string)
        self.text_space.see('end')


class MyApplication:

    def __init__(self):
        self.about_dialog = None
        self.fig = None
        self.builder = b = pygubu.Builder()
        b.add_from_file(os.path.join(DATA_DIR, 'main.ui'))
        b.add_resource_path(os.path.join(DATA_DIR, 'imgs'))
        #self.canvas = b.get_object('main_canvas')

        self.mainwindow = b.get_object('mainwindow')
        self.mainwindow.iconbitmap('.\imgs\molo_256.ico')
        #print(dir(self.mainwindow))

     # Connect to Delete event
        self.mainwindow.protocol("WM_DELETE_WINDOW", self.quit)

        b.connect_callbacks(self)


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

    def quit(self, event=None):
        self.mainwindow.quit()

    def btnplot_clicked(self):
        if not self.fig==None:
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

    def create_model(self):
        self.text_box = self.builder.get_object('model_text')
        self.text_box.delete('1.0', '2.0')
        sys.stdout = StdoutRedirector(self.text_box)
        analyses_root=Path(self.builder.get_object('analyses_root_input').get())
        print(analyses_root)

        park_label =self.builder.get_object('park_label_input').get()
        wtg_label =self.builder.get_object('wtg_label_input').get()
        template_dir = Path(os.getcwd()).parents[0].joinpath('Optimization').joinpath('templates')

        p = Parameter_Space(templates_dir=template_dir)
        p.ncol = 2
        p.gap = 1.1
        p.height = 21
        p.column_diameter = 7.5
        p.filling_ratio = [0.1, 0.1]
        this_candidate = Candidate(analyses_root, park_label, wtg_label, p, case_label_type='molo_model')
        this_candidate.settings.wtg_model = "Vestas 9.5"
        this_candidate.init_model(state='New')
        with open("this_candidate.pkl", "wb") as f:
            pickle.dump(this_candidate, f)

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
    #root = tk.Tk()
    app = MyApplication()
    app.run()
