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
import numpy as np
from pyNemoh_root.pyNemoh import utility
import h5py

#import module_locator

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



        self.nautic_zones_gif= PhotoImage(file=r".\imgs\nautic_zones.gif")

        self.load_cfg()
        #self.create_model()
        cb = self.builder.get_object('plot_pressure_type_input')
        #cb = self.builder.get_object('plot_pressure_type_input')
        #cb.delete(0, END)
        #cb["values"] = '\"{:s}\"'.format('\" \"'.join(list(self.this_candidate.loads.pressure.keys())))

        #cb["values"] = list(self.this_candidate.loads.pressure.keys())
        #print(list(self.this_candidate.loads.pressure.keys()))
        cb.current(0)




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
            canvas.create_image(0,0, anchor = NW, image = self.nautic_zones_gif)

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
        # Read global config
        f_global = Path(DATA_DIR).joinpath('mdo_global_cfg.json')

        with open(f_global, 'r') as f:
            self.global_cfg = json.loads(f.read())

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

            self.case_cfg=self.default_case_cfg

        else:
            self.case_cfg=self.local_case_cfg

            # Loop trough items in default case config and add missing items
            for key in self.default_case_cfg['user_input']['object']:
                if key not in self.case_cfg['user_input']['object']:
                    self.case_cfg[key]=self.default_case_cfg['user_input']['object'][key]


        objects = self.global_cfg['user_input']['object']
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
        template_dir = Path(os.getcwd()).joinpath('templates')

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
        self.this_candidate = Candidate(analyses_root_input, park_label_input, wtg_label_input, p,
                                        case_label_type='molo_model')
        self.this_candidate.settings.wtg_model = "Vestas 9.5"
        self.this_candidate.settings.wtg_model = "Haliade X"
        self.this_candidate.settings.wtg_model = "Generic 8MW"


        self.this_candidate.init_model(state='Old')
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
        yr = float(self.builder.get_object('calc_response_yr_input').get())  # Return period in years, statistics conditioned for 3hr storms

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

        bool_lin=self.builder.tkvariables.__getitem__('lin_visc_damp_var').get()
        tc.settings.do_linearize = bool_lin

        if tc.settings.do_linearize:
            print('\nLinearized viscous damping will be included')
        else:
            print('\nLinearized viscous damping will NOT be included')

        print('\nCase:\t{}'.format(self.this_candidate.settings.case_label))

        print('\nEigenvalue sollution WITH added mass')

        tb.eigenvalprint(tc.loads.m + tc.loads.ma[0, :, :], tc.loads.k)

        d_col_central=tc.settings.floater_data['Central column diameter']
        d_col_radial=tc.settings.floater_data['Radial']['Column']['Diameter']
        n_col_radial=tc.settings.floater_data['Radial']['Number of columns']
        gap = tc.settings.floater_data['Gap factor']
        height=tc.settings.floater_data['Radial']['Heigth']
        draught=tc.hs_floater.hs_data['draught']

        # Get section forces
        sp_x =  d_col_central* (0.5)  # + 0.8 + 1 + 0.8 + 1)
        sp_z = height / 2 - draught
        # sp_x = -100
        sp_z = 0
        sp1 = [sp_x, 0, sp_z]  # Used for moment reference
        sn1 = [1, 0, 0]

        radial_extreme=d_col_central* (0.5)+(1+gap)*n_col_radial*d_col_radial

        print('\nAirgap point at {:1.2f}'.format(radial_extreme))

        c1 = np.asarray([[radial_extreme, 0, 0]])

        with h5py.File(tc.settings.fio.nemoh_root.joinpath('db.hdf5'), "r") as hdf5_db:
            environment = utility.read_environment(hdf5_db)

        k_wave=np.zeros(tc.loads.nw, dtype=float)
        for iw, val in enumerate(tc.loads.w):
            k_wave[iw] = utility.compute_wave_number(val, environment)
        w_bar = (radial_extreme - environment.x_eff) * np.cos(tc.loads.beta) + (0 - environment.y_eff) * np.sin(tc.loads.beta[ibeta])
        eta=np.exp(utility.II * k_wave * w_bar)








        print('\n {:^6s} {:^6s} {:^6s}  {:^6s}  {:^6s}  {:^6s}'.format('Hs', 'Tp', 'Gamma', force_label[2], force_label[4], 'AG'))
        for stwc in tc.contourline:
            tc.init_response(short_term_wave_condition=stwc)



            imass, ipanel, istrip = tc.response.get_section_index(sp1, sn1)
            f_sec1 = tc.response.assemble_forces(imass, ipanel, istrip, moment_ref_point=[0, 0, 0])
            del imass, ipanel

            fz = f_sec1['Dynamic']['SUM'][:, ibeta, 2] / 1000000
            my = f_sec1['Dynamic']['SUM'][:, ibeta, 4] / 1000000

            fz_elm = stwc.expected_largest_maximum(fz, tc.loads.w)
            my_elm = stwc.expected_largest_maximum(my, tc.loads.w)



            p1_rao=tc.response.point_rao(c1)[:, ibeta, 0,2].flatten()
            ag_rao = eta+p1_rao
            ag1 = stwc.expected_largest_maximum(ag_rao, tc.loads.w)

            print(' {:6.1f} {:6.1f} {:6.1f}  {:6.2f}  {:6.1f}  {:6.1f} '.format(stwc.hs, stwc.tp, stwc.gamma, fz_elm, my_elm, ag1))

    def plot_pressure(self):
        self.get_output('plot_text')
        ifreq = int(self.builder.get_object('plot_pressure_ifreq_input').get())
        pressure_index = int(self.builder.get_object('plot_pressure_index_input').get())
        pressure_type = self.builder.get_object('plot_pressure_type_input').get()
        axis = int(self.builder.get_object('plot_pressure_axis_input').get())
        self.this_candidate.loads.show_pressure(ifreq, pressure_index, axis, pressure_type)

    def gui_init_load(self):
        if not self.this_candidate == None:
            if self.this_candidate.loads == None:
                self.this_candidate.init_load()
                print('Loads initialized')

        else:
            print('self.this_candidate == None')

    def gui_init_plot(self, event=None):
        self.get_output('plot_text')
        self.gui_init_load()

    def gui_init_response(self, event=None):

        self.gui_init_load()










if __name__ == '__main__':
    # root = tk.Tk()
    app = Application()
    app.run()
