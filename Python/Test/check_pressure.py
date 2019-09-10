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

        this_candidate.init_model(state=state)
        this_candidate.init_load_response()
        print('\nCase:\t{}'.format(this_candidate.settings.case_label))


        # this_candidate.loads.show_pressure(ifreq=20, pressure_index=2, pressure_type='Radiation', axis=2)

        this_candidate.loads.show_pressure(ifreq=20, pressure_index=0, pressure_type='Froude-Krylof', axis=2)

        # this_candidate.loads.show_pressure(ifreq=20, pressure_index=0, pressure_type='Diffraction', axis=2)
