import json
from pathlib import Path

import matplotlib.pyplot as plt  # noqa
import numpy as np

import sea_loads
import model_set_up as ms
import nemoh
from common import FileIO

np.set_printoptions(precision=3)
step1 = True
step2 = True
step3 = True

# Give the analysed DOFs
dof = [0, 0, 1, 0, 1, 0]
rho_sw = 1025
water_depth = 100
omega_calc = [41, 0.2, 20]

root = Path(r'C:\Users\eison\OneDrive - Verbun AS\Divisions\Offshore Wind\Projects\P2017.001\Work\MOLO_Analyses')

if step1:
    fio = FileIO(root)

    model_template = fio.templates_dir.joinpath('model_template.json')

    with open(model_template, 'r') as f:
        data = json.loads(f.read())

    model_input = fio.data_io_dir.joinpath('model.json')
    with open(model_input, 'w') as f:
        f.write(json.dumps(data, indent=4, sort_keys=True))

    unit_model, hs_floater = ms.launch(fio, model_input)

    # yz_section = (unit_model.parts_list[1].gap + 1) * unit_model.parts_list[1].dia_rc / 2

    K33 = hs_floater.hs_data['stiffness_matrix'][0, 0]
    K44 = hs_floater.hs_data['stiffness_matrix'][1, 1]
    M33 = unit_model.inertias.mass_matrix_global[2][2]
    M44 = unit_model.inertias.mass_matrix_global[3][3]
    T33 = 2 * np.pi / (np.sqrt(K33 / M33))
    T44 = 2 * np.pi / (np.sqrt(K44 / M44))
    print('\n\nT33 = {:>4.1f} s'.format(T33))
    print('T44 = {:>4.1f} s'.format(T44))

    hydro_mesh_symmetri, mesh_file = nemoh.mesh(fio.nemoh_dir, hs_floater, sym=0)

    hydro_mesh_symmetri.show()

    if step2:
        nemoh.runNemoh(fio.nemoh_dir, hydro_mesh_symmetri, mesh_file, rho_sw, water_depth, omega_calc, dof)

else:

    fio = FileIO(root, case_number=16)

if step3:
    # Read results, perform postprocessing and write pdf
    w = nemoh.getOmega(fio.nemoh_results)
    dir = nemoh.getDirections(fio.nemoh_results)

    hdp = sea_loads.HydroCoefficients(dof, w, dir, fio)

    print('\nSection forces')
    for i in range(6):
        print('\tDOF{}: {:8.1f}'.format(i + 1, np.sqrt(
            abs(sum(sea_loads.spec_response(hdp.f_sec[:, i][::3], hdp.w, hs=10, wp=2 * np.pi / 14))))))
    # write_report(root, hdp, case_label)

    # plt.plot(2*np.pi/w,hdp._c_hyd[0,0,:,0])
    # plt.plot(2*np.pi/w,hdp._fe_amp[0,:,0])
    plt.plot(2 * np.pi / w, hdp.getRAO(5, 0))
    plt.show()

    dof = 0
    # print(hdp.ma[0,:,0,0])
    print(hdp._fe_amp[0, :, dof])
    print(hdp._c_hyd[0, dof, :, dof])
    print(hdp._ma[0, dof, :, dof])
    print(hdp._m[0, 0])
