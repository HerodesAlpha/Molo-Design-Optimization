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
        sp_z = this_candidate.settings.floater_data['Radial']['Heigth'] / 2 - this_candidate.hs_floater.hs_data[
            'draught']
        # sp_x = - 1000
        section_point = [sp_x, 0, sp_z]  # Used for moment reference
        section_normal = [1, 0, 0]

        imass, ipanel = this_candidate.response.get_section_index(section_point, section_normal)
        this_candidate.f_sec1 = this_candidate.response.assemble_forces(imass, ipanel, moment_ref_point=section_point)
        del imass, ipanel

        x = this_candidate.response._panel_pressure_centers
        y = np.squeeze(this_candidate.response._projected_panel_area)

        import csv

        with open("data.csv", 'w', newline='') as f:
            writer = csv.writer(f, dialect='excel')
            writer.writerows(np.hstack((x, y)))

        #
