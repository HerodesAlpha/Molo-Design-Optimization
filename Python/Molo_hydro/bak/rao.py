import json
from pathlib import Path

import numpy as np

import sea_loads
# import print_plot pp
import model_set_up as ms
import nemoh
from report import write_report

step1 = False
step2 = False
step3 = True

# Give the analysed DOFs
dof = [0, 0, 1, 0, 1, 0]
rho_sw = 1025
water_depth = 100
omega_calc = [41, 0.2, 20]

root = Path(r'C:\Users\eison\OneDrive - Verbun AS\Divisions\Offshore Wind\Projects\P2017.001\Work\MOLO_Analyses')

if step1:
    i = 0
    while 1:
        i += 1
        if not root.joinpath('case{:04d}'.format(i)).exists():
            break
    case_label = 'case{:04d}'.format(i)

    root = root.joinpath(case_label)

    nemoh_path = root.joinpath('nemoh')
    nemoh_path.mkdir(parents=True, exist_ok=False)

    model_template = Path(r'.\templates').joinpath('model_template.json')
    with open(model_template, 'r') as f:
        data = json.loads(f.read())

    model_input = root.joinpath('model.json')
    with open(model_input, 'w') as f:
        f.write(json.dumps(data, indent=4, sort_keys=True))

    unit_model, hs_floater = ms.launch(root, model_input)

    K33 = hs_floater.hs_data['stiffness_matrix'][0, 0]
    K44 = hs_floater.hs_data['stiffness_matrix'][1, 1]
    M33 = unit_model.inertias.mass_matrix_global[2][2]
    M44 = unit_model.inertias.mass_matrix_global[3][3]
    T33 = 2 * np.pi / (np.sqrt(K33 / M33))
    T44 = 2 * np.pi / (np.sqrt(K44 / M44))
    print('\n\nT33 = {:>4.1f} s'.format(T33))
    print('T44 = {:>4.1f} s'.format(T44))

    hydro_mesh_symmetri, mesh_file = nemoh.mesh(nemoh_path, hs_floater)

    if step2:
        nemoh.runNemoh(nemoh_path, hydro_mesh_symmetri, mesh_file, rho_sw, water_depth, omega_calc, dof)

else:
    case_select = 16
    case_label = 'case{:04d}'.format(case_select)
    root = root.joinpath(case_label)
    nemoh_path = root.joinpath('nemoh')

if step3:
    # Read results, perform postprocessing and write pdf
    w = nemoh.getOmega(nemoh_path)
    dir = nemoh.getDirections(nemoh_path)
    print(w)
    hdp = sea_loads.HydroCoefficients(dof, w, dir, root, nemoh_path)
    write_report(root, hdp, case_label)
