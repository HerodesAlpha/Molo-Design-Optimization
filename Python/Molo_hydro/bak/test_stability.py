import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import helper
import sea_loads
import model_set_up as ms
import nemoh
from common import FileIO
from stability import gz_curve

np.set_printoptions(precision=3)
step1 = False
step2 = False
step3 = True

ROOT = Path(r'C:\Users\EDUX\OneDrive - Verbun AS\Divisions\Offshore Wind\Projects\P2017.001\Work\MOLO_Analyses')
ROOT = Path(r'C:\Users\eison\OneDrive - Verbun AS\Divisions\Offshore Wind\Projects\P2017.001\Work\MOLO_Analyses_Stability')

TEMPLATES_DIR = Path(r'.\templates')

model_template = TEMPLATES_DIR.joinpath('model_template.json')

with open(model_template, 'r') as f:
    data = json.loads(f.read())

data["Analyses parameters"]["Degrees of Freedom"] = [0, 0, 1, 0, 1, 0]
data["Analyses parameters"]["Number of wave directions, Min and Max (degrees)"] = [1, 0, 0]
data["Analyses parameters"]["Density of sea water"] = 1025
data["Analyses parameters"]["Water depth"] = 100
data["Analyses parameters"]["Number of wave frequencies, Min, and Max (rad/s)"] = [81, np.pi / 20, np.pi]
data["Analyses parameters"]["Use symmetri"] = 1

with open(model_template, 'w') as f:
    f.write(json.dumps(data, indent=4, sort_keys=True))

NEMOH_DOF = data["Analyses parameters"]["Degrees of Freedom"]
NEMOH_DIR = data["Analyses parameters"]["Number of wave directions, Min and Max (degrees)"]
RHO_SW = data["Analyses parameters"]["Density of sea water"]
WATER_DEPTH = data["Analyses parameters"]["Water depth"]
OMEGA_NEMOH_INP = data["Analyses parameters"]["Number of wave frequencies, Min, and Max (rad/s)"]
SYM = data["Analyses parameters"]["Use symmetri"]

fio = FileIO(ROOT, TEMPLATES_DIR)
model_input = fio.data_io_dir.joinpath('model.json')
with open(model_input, 'w') as f:
    f.write(json.dumps(data, indent=4, sort_keys=True))

unit_model, hs_floater = ms.launch(fio, model_input)

gz_curve(hs_floater)


