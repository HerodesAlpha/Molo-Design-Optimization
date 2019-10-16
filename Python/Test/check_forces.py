from pathlib import Path
import numpy as np
from design_engine import Candidate, Parameter_Space
import sys
import warnings
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
    p.height = 21
    p.column_diameter = 7.5
    p.filling_ratio = [0] * 3
    this_candidate = Candidate(analyses_root, park_label, wtg_label, p, case_label_type='molo_model')
    this_candidate.settings.wtg_model = "Vestas 9.5"
    with open(this_candidate.settings.fio.case_dir.joinpath('stdout_redirect.txt'), 'w') as fout:
        #this_candidate.settings.wtg_model = "Vestas 9.5"

        if not this_candidate.has_model():
            state = 'New'
            print('No model found')

        else:
            state = "Old"
            print('Has old model')

        #state = "Old"

        with redirect_stdout(fout):
            this_candidate.init_model(state=state)
        print('\nCase:\t{}'.format(this_candidate.settings.case_label))

        this_candidate.settings.critical_damping_ratio = 0

        this_candidate.init_load_response()

        # Consider first radial

        # Get section forces
        sp_x = this_candidate.settings.floater_data['Central column diameter'] * (0.5)  # + 0.8 + 1 + 0.8 + 1)
        sp_z = this_candidate.settings.floater_data['Radial']['Heigth'] / 2 - this_candidate.hs_floater.hs_data[
            'draught']
        sp_x = -100
        sp_z = 0
        section_point = [sp_x, 0, sp_z]  # Used for moment reference
        section_normal = [1, 0, 0]

        imass, ipanel = this_candidate.response.get_section_index(section_point, section_normal)
        this_candidate.f_sec1 = this_candidate.response.assemble_forces(imass, ipanel, moment_ref_point=[0, 0, 0])
        del imass, ipanel

        ibeta = 0
        idof = 2

        stat_force = this_candidate.f_sec1['Static']

        factor = 1 / 1000000

        for key in stat_force:
            a = np.abs(stat_force[key][2]) * factor
            b = np.angle(stat_force[key][2])
            print('{}: {:1.2f} {:1.2f}'.format(key, a, b))

        w = this_candidate.loads.w
        fig, axs = plt.subplots(2, 2)
        dyn_force = this_candidate.f_sec1['Dynamic']

        x_tics = 1/(2 * np.pi / w)
        # x_label = 'Period [s]'
        x_label = 'Frequency [Hz]'

        for key in dyn_force:
            if not key == 'whatever':
                axs[0, 0].plot(x_tics, np.abs(dyn_force[key][:, ibeta, idof]) * factor, label=key)
                axs[0, 1].plot(x_tics, np.angle(dyn_force[key][:, ibeta, idof]), label=key)

        axs[0, 0].set_title('Amplitude')
        force_label = ['Fx [MN]', 'Fy [MN]', 'Fz [MN]', 'Mx [MNm]', 'My [MNm]', 'Mz [MNm]']
        axs[0, 0].set_ylabel(force_label[idof])

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

        d = {'Heave': 2, 'Pitch': 4}
        ax1 = axs[1, 0]
        ax2 = ax1.twinx()
        for key in d:

            if key == 'Heave':
                lns1 = ax1.plot(x_tics, np.abs(this_candidate.response.rao[:, ibeta, d[key]]), label=key)

            else:
                lns2 = ax2.plot(x_tics, np.abs(this_candidate.response.rao[:, ibeta, d[key]]), '-r', label=key)

            axs[1, 1].plot(x_tics, np.angle(this_candidate.response.rao[:, ibeta, d[key]]), label=key)

        lns = lns1 + lns2
        labs = [l.get_label() for l in lns]
        axs[1, 0].legend(lns, labs, loc=0)
        axs[1, 1].legend()
        plt.suptitle('{}\nSection point: [{sp[0]:1.2f}, {sp[1]:1.2f}, {sp[2]:1.2f}]\nSection normal: [{sn[0]:1.2f}, {sn[1]:1.2f}, {sn[2]:1.2f}]'.format(this_candidate.settings.case_label,sp=section_point,sn=section_normal))
        plt.show()
