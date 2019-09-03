from pathlib import Path
import numpy as np
from design_engine import Candidate, Parameter_Space
import sys
import warnings
import tool_box as tb
import pandas
from contextlib import redirect_stdout
import os

import matplotlib.pyplot as plt
class BreakIt(Exception): pass


if __name__ == '__main__':
    if not sys.warnoptions:
        warnings.simplefilter("default")

    analyses_root = Path(r'C:\MOLO_Optimization')
    park_label = 'site_01'
    wtg_label = 'wtg_01'
    candidate_parent_dir = analyses_root.joinpath(park_label).joinpath(wtg_label)


    template_dir = Path(os.getcwd()).parents[0].joinpath('Optimization').joinpath('templates')

    p = Parameter_Space(templates_dir=template_dir)

    p.ncol = 3

    is_stable = False
    #    for d in np.linspace(8.8, 8.8, 1, dtype=float):
    irow = -1

    p.gap = 0.8
    p.height = 20
    p.column_diameter = 8.1
    p.filling_ratio = [0] * 3
    this_candidate = Candidate(analyses_root, park_label, wtg_label, p, case_label_type='molo_model')
    with open(this_candidate.settings.fio.case_dir.joinpath('stdout_redirect.txt'), 'w') as fout:

        if not this_candidate.has_model():
            state = 'New'

        else:
            state = "Old"

        print(state)
        with redirect_stdout(fout):
            this_candidate.init_model(state=state)
        print('\nCase:\t{}'.format(this_candidate.settings.case_label))
        this_candidate.init_load_response()

        # Consider first radial

        # Get section forces
        sp_x = this_candidate.settings.floater_data['Central column diameter'] * (0.5)  # + 0.8 + 1 + 0.8 + 1)
        sp_z = this_candidate.settings.floater_data['Radial']['Heigth'] / 2 - this_candidate.hs_floater.hs_data['draught']
        section_point = [sp_x, 0, sp_z]  # Used for moment reference
        section_normal = [1, 0, 0]

        imass, ipanel = this_candidate.response.get_section_index(section_point, section_normal)
        this_candidate.f_sec1 = this_candidate.response.assemble_forces(imass, ipanel, moment_ref_point=section_point)
        del imass, ipanel

        w = this_candidate.loads.w
        f_fk = this_candidate.f_sec1['Dynamic']['Froude-Krylof'][:,0,2]/9810
        f_rad = this_candidate.f_sec1['Dynamic']['Radiation'][:,0,2]/9810
        f_diff = this_candidate.f_sec1['Dynamic']['Diffraction'][:,0,2]/9810
        f_bouy = this_candidate.f_sec1['Dynamic']['Buoyancy'][:,0,2]/9810
        f_inert = this_candidate.f_sec1['Dynamic']['Inertia'][:,0,2]/9810

        fig, axs = plt.subplots(2)
        axs[0].plot(w, abs(f_fk), label='Froude-Krylof')
        axs[0].plot(w, abs(f_rad), label='Diffraction')
        axs[0].plot(w, abs(f_diff), label='Radiation')
        axs[0].plot(w, abs(f_bouy), label='Buoyancy')
        axs[0].plot(w, abs(f_inert), label='Inertia')
        axs[0].set_title('Magnitude')
        axs[0].legend()
        axs[1].plot(w, np.angle(f_fk), label='Froude-Krylof')
        axs[1].plot(w, np.angle(f_rad), label='Diffraction')
        axs[1].plot(w, np.angle(f_diff), label='Radiation')
        axs[1].plot(w, np.angle(f_bouy), label='Buoyancy')
        axs[1].plot(w, np.angle(f_inert), label='Inertia')
        axs[1].set_title('Phase')
        axs[1].legend()
        plt.show()
        #
        # fig, axs = plt.subplots(2)
        # axs[0, 0].plot(w, abs(f_fk), 'tab:blue', label='Froude-Krylof')
        # axs[0, 0].plot(w, abs(f_rad), 'tab:green', label='Diffraction')
        # axs[0, 0].plot(w, abs(f_diff), 'tab:orange', label='Radiation')
        # axs[0, i].plot(w, abs(f_bouy), 'tab:yellow', label='Buoyancy')
        # axs[0, i].plot(w, abs(f_inert), 'tab:red', label='Inertia')
        # axs[0, i].set_title('All forces')
        # axs[0, i].legend()


        # force_out = dict()
        # force_out['Static'] = dict()
        # force_out['Dynamic'] = dict()
        #
        # force_out['Static']['Total'] = f_gravity + f_bouyancy
        # force_out['Static']['Gravity'] = f_gravity
        # force_out['Static']['Buoyancy'] = f_bouyancy
        # force_out['Dynamic']['Total'] = f_rad + f_fk + f_diff + f_dz_s + f_inertia
        # force_out['Dynamic']['Radiation'] = f_rad
        # force_out['Dynamic']['Froude-Krylof'] = f_fk
        # force_out['Dynamic']['Diffraction'] = f_diff
        # force_out['Dynamic']['Buoyancy'] = f_dz_s
        # force_out['Dynamic']['Inertia'] = f_inertia
