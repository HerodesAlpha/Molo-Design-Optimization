"""
Build model script for creating and saving a design candidate.

This script creates a design candidate with specified parameters
and saves it as a pickle file for later use.
"""

import os
import pickle
from pathlib import Path

from optimization.design_engine import Candidate, Parameter_Space


if __name__ == '__main__':
    analyses_root = Path(r'C:\mdo_working_dir')
    park_label = 'site_01'
    wtg_label = 'test'
    candidate_parent_dir = analyses_root.joinpath(park_label).joinpath(wtg_label)

    p = Parameter_Space(templates_dir=None)

    p.ncol = 2

    is_stable = False
    #    for d in np.linspace(8.8, 8.8, 1, dtype=float):
    irow = -1

    p.gap = 1.7
    p.height = 21
    p.radial_column_diameter = 7.8
    p.filling_ratio = [0.1, 0.1]
    this_candidate = Candidate(analyses_root, park_label, wtg_label, p, case_label_type='molo_model')
    this_candidate.settings.wtg_model = "Vestas 9.5"
    # this_candidate.settings.wtg_model = "Haliade X"
    this_candidate.init_model(state='New')
    with open("this_candidate.pkl", "wb") as f:
        pickle.dump(this_candidate, f)
