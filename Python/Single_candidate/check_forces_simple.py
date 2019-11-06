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
    this_candidate = Candidate(analyses_root, park_label, wtg_label, p, case_label_type='test_simple')

    this_candidate.settings.mesh_name = this_candidate.settings.case_label
    this_candidate.settings.use_symmmetri = False
    this_candidate.settings.lower_face_corrected_z_pos = False
    with open(this_candidate.settings.fio.case_dir.joinpath('stdout_redirect.txt'), 'w') as fout:

        if not this_candidate.has_model():
            state = 'New'

        else:
            state = "Old"

        print(state)
        with redirect_stdout(fout):
            this_candidate.init_model(state=state)

        print('\nCase:\t{}'.format(this_candidate.settings.case_label))

        # Run Nemoh
        this_candidate.settings.load_cases = {  # 121, np.pi / 15, np.pi
                "num_wave_frequencies": 40,  # TODO: Implement adaptive frequency
                "min_wave_frequencies": 2 * np.pi / 30,  # (rad/s)
                "max_wave_frequencies": 2 * np.pi / 4,
                "num_wave_directions" : 3,
                "min_wave_directions" : 0,  # deg
                "max_wave_directions" : 90,
        }
        this_candidate.hydrodynamic_analysis()

        this_candidate.settings.critical_damping_ratio = 0.01
        this_candidate.init_load_response()

        # Consider first radial

        # Get section forces
        sp_x = this_candidate.settings.floater_data['Central column diameter'] * (0.5)  # + 0.8 + 1 + 0.8 + 1)
        sp_z = this_candidate.settings.floater_data['Radial']['Heigth'] / 2 - this_candidate.hs_floater.hs_data[
            'draught']
        sp_x = - 1000
        section_point = [sp_x, 0, sp_z]  # Used for moment reference
        section_normal = [1, 0, 0]

        imass, ipanel = this_candidate.response.get_section_index(section_point, section_normal)
        this_candidate.f_sec1 = this_candidate.response.assemble_forces(imass, ipanel, moment_ref_point=[0, 0, 0])
        del imass, ipanel

        ibeta = 0
        idof = 2

        stat_force = this_candidate.f_sec1['Static']
        stat_force['Gravity'] /= 7000

        factor = 1

        for key in stat_force:
            a = np.abs(stat_force[key][2]) * factor
            b = np.angle(stat_force[key][2])
            print('{}: {:1.2f} {:1.2f}'.format(key, a, b))

        w = this_candidate.loads.w
        fig, axs = plt.subplots(2,2)
        dyn_force = this_candidate.f_sec1['Dynamic']
        dyn_force['Inertia'] /=10 #
        dyn_force['Buoyancy'] *=10 #

        for key in dyn_force:
            if not key == 'Sum':
                axs[0,0].plot(w, np.abs(dyn_force[key][:, ibeta, idof]) * factor, label=key)
                axs[0,1].plot(w, np.angle(dyn_force[key][:, ibeta, idof]), label=key)

        axs[0,0].set_title('Magnitude')
        axs[0,0].legend()
        axs[0,1].set_title('Phase')
        axs[0,1].legend()



        axs[1, 0].plot(w, np.abs(this_candidate.response.rao[:, ibeta, idof]), label='RAO')
        axs[1, 1].plot(w, np.angle(this_candidate.response.rao[:, ibeta, idof]), label='RAO')

        plt.show()
        #
