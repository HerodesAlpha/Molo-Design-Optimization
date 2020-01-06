# file: myapp.py
import sys
import os
import random
import math
from tkinter import *
from tkinter import filedialog
import tkinter as tk

import pygubu
from plot_forces import get_axs
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

try:
    DATA_DIR = os.path.abspath(os.path.dirname(__file__))
except NameError:
    DATA_DIR = os.path.dirname(os.path.abspath(sys.argv[0]))

DEG2RAD = 4 * math.atan(1) * 2 / 360


class MyApplication:

    def __init__(self):
        self.about_dialog = None
        self.fig = None
        self.builder = b = pygubu.Builder()
        b.add_from_file(os.path.join(DATA_DIR, 'main_test.ui'))
        b.add_resource_path(os.path.join(DATA_DIR, 'imgs'))
        #self.canvas = b.get_object('main_canvas')

        self.mainwindow = b.get_object('main_window')
        self.mainwindow.iconbitmap('.\imgs\molo_256.ico')
        #print(dir(self.mainwindow))

     # Connect to Delete event
        #self.mainwindow.protocol("WM_DELETE_WINDOW", self.quit)

        b.connect_callbacks(self)

        print('Plotting ...')

        get_axs(self)
        master = self.builder.get_object('main_canvas')

        canvas = FigureCanvasTkAgg(figure=self.fig, master=master)
        canvas.draw()



    def run(self):
        self.mainwindow.mainloop()

if __name__ == '__main__':
    #root = tk.Tk()
    app = MyApplication()
    app.run()
