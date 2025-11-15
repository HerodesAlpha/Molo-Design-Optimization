"""
Main test application for plotting forces using pygubu GUI builder.

This script creates a GUI application that displays force and RAO plots
using matplotlib embedded in a tkinter window.
"""

import os
import sys
import math
from pathlib import Path

import pygubu
import tkinter as tk
from tkinter import filedialog
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

from plot_forces import get_axs

try:
    DATA_DIR = os.path.abspath(os.path.dirname(__file__))
except NameError:
    DATA_DIR = os.path.dirname(os.path.abspath(sys.argv[0]))

DEG2RAD = math.pi / 180.0


class MyApplication:
    """
    Main application class for force plotting GUI.
    """
    
    def __init__(self):
        """Initialize the application and create GUI."""
        self.about_dialog = None
        self.fig = None
        self.builder = b = pygubu.Builder()
        b.add_from_file(os.path.join(DATA_DIR, 'main_test.ui'))
        b.add_resource_path(os.path.join(DATA_DIR, 'imgs'))
        # self.canvas = b.get_object('main_canvas')

        self.mainwindow = b.get_object('main_window')
        icon_path = os.path.join(DATA_DIR, 'imgs', 'molo_256.ico')
        if os.path.exists(icon_path):
            self.mainwindow.iconbitmap(icon_path)
        # print(dir(self.mainwindow))

        # Connect to Delete event
        # self.mainwindow.protocol("WM_DELETE_WINDOW", self.quit)

        b.connect_callbacks(self)

        print('Plotting ...')

        get_axs(self)
        master = self.builder.get_object('main_canvas')

        canvas = FigureCanvasTkAgg(figure=self.fig, master=master)
        canvas.draw()

    def run(self):
        """Start the main event loop."""
        self.mainwindow.mainloop()


if __name__ == '__main__':
    # root = tk.Tk()
    app = MyApplication()
    app.run()
