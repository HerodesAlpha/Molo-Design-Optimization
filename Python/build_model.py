from pathlib import Path
import pickle
from core.design_engine import Candidate, Parameter_Space
import os


analyses_root = Path(r'C:\mdo_working_dir')
park_label = 'site_01'
wtg_label = 'test'
candidate_parent_dir = analyses_root.joinpath(park_label).joinpath(wtg_label)

template_dir = Path(os.getcwd()).parents[0].joinpath('optimization').joinpath('templates')

p = Parameter_Space(templates_dir=template_dir)

p.ncol = 2

is_stable = False
#    for d in np.linspace(8.8, 8.8, 1, dtype=float):
irow = -1

p.gap = 1.1
p.height = 21
p.column_diameter = 7.5
p.filling_ratio = [0.1, 0.1]
this_candidate = Candidate(analyses_root, park_label, wtg_label, p, case_label_type='molo_model')
this_candidate.settings.wtg_model = "Vestas 9.5"
# this_candidate.settings.wtg_model = "Haliade X"
this_candidate.init_model(state='New')
with open( "this_candidate.pkl", "wb" ) as f:
    pickle.dump( this_candidate,  f)