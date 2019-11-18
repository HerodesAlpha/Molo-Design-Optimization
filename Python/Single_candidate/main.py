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
        self.ax = None
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


    def plot(self):

        def btnplot_clicked():
            if not self.ax==None:
                self.ax.clear()
            print('Plotting...')
            self.fig, self.axs = get_axs()
            master = self.builder.get_object('main_canvas')
            self.canvas = FigureCanvasTkAgg(self.fig, master=master)
            self.canvas.draw()

        btnplot = self.builder.get_object('plot_button')
        btnplot['command'] = btnplot_clicked





    def run(self):
        self.mainwindow.mainloop()


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
