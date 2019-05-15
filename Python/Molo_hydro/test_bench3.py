__author__ = "Eivind Sonju"
__copyright__ = "Copyright (C) 2017-2019 Verbun AS. All rights reserved."
__version__ = "2.0"

import json
from pathlib import Path
import pickle

import matplotlib.pyplot as plt
import numpy as np

import helper
import calculations
import model_set_up as ms
import nemoh
from common import FileIO


import stability

create_model = True
run_nemoh = True
calc_gz = True
postprocessing = True

ROOT = Path(r'C:\analyses')

TEMPLATES_DIR = Path(r'.\templates')

model_template_path = TEMPLATES_DIR.joinpath('model_template.json')

with open(model_template_path, 'r') as f:
    data = json.loads(f.read())

data["Analyses parameters"]["Degrees of Freedom"] = [0, 0, 1, 1, 1, 0]
data["Analyses parameters"]["Number of wave directions, Min and Max (degrees)"] = [3, 0, 90]
data["Analyses parameters"]["Density of sea water"] = 1025
data["Analyses parameters"]["Water depth"] = 100
data["Analyses parameters"]["Number of wave frequencies, Min, and Max (rad/s)"] = [121, np.pi / 15, np.pi]
data["Analyses parameters"]["Use symmetri"] = 1

with open(model_template_path, 'w') as f:
    f.write(json.dumps(data, indent=4, sort_keys=True))

NEMOH_DOF = data["Analyses parameters"]["Degrees of Freedom"]
NEMOH_DIR = data["Analyses parameters"]["Number of wave directions, Min and Max (degrees)"]
RHO_SW = data["Analyses parameters"]["Density of sea water"]
WATER_DEPTH = data["Analyses parameters"]["Water depth"]
OMEGA_NEMOH_INP = data["Analyses parameters"]["Number of wave frequencies, Min, and Max (rad/s)"]
SYM = data["Analyses parameters"]["Use symmetri"]

if create_model:
    CASE_NUMBER = None
    fio = FileIO(ROOT, TEMPLATES_DIR, CASE_NUMBER)
    model_input = fio.data_io_dir.joinpath('model.json')
    with open(model_input, 'w') as f:
        f.write(json.dumps(data, indent=4, sort_keys=True))

    unit_model, hs_floater = ms.launch(fio, model_input)
    pickle.dump(unit_model, open(fio.data_io_dir.joinpath('unit_model.pkl'), 'wb'))
    pickle.dump(hs_floater, open(fio.data_io_dir.joinpath('hs_floater.pkl'), 'wb'))

    hydro_mesh_symmetri, mesh_file = nemoh.mesh(fio, hs_floater, SYM)

    hydro_mesh_symmetri.show()


else:
    CASE_NUMBER = 2
    SYM = 1
    fio = FileIO(ROOT, TEMPLATES_DIR, CASE_NUMBER)
    unit_model = pickle.load(open(fio.data_io_dir.joinpath('unit_model.pkl'), 'rb'))
    hs_floater = pickle.load(open(fio.data_io_dir.joinpath('hs_floater.pkl'), 'rb'))

if create_model and calc_gz:
    stability.gz_curve(fio, hs_floater)

if create_model and run_nemoh:
    nemoh.runNemoh(fio, hydro_mesh_symmetri, mesh_file, NEMOH_DIR, RHO_SW, WATER_DEPTH, OMEGA_NEMOH_INP,
                   NEMOH_DOF)

if postprocessing:
    # Read results, perform postprocessing and write pdf
    w = nemoh.getOmega(fio.nemoh_results)
    dir = nemoh.getDirections(fio.nemoh_results)

    #hdp = sea_loads.HydroCoefficients(NEMOH_DOF, w, dir, fio, SYM)
    hdp = calculations.TransferFunctions(NEMOH_DOF, w, dir, fio, SYM)
    SELECT_DOF = 1
    SELECT_AXIS = 2
    SELECT_FREC = 1

    # print('\nSection forces')
    # for i in range(6):
    #     print('\tDOF{}: {:8.1f}'.format(i + 1, np.sqrt(
    #         abs(sum(hydro.spec_response(hdp.f_sec[:, i][::3], hdp.w, hs=10, wp=2 * np.pi / 14))))))
    # write_report(root, hdp, case_label)

    ifreq = 0
    idir = 0
    irad = 0
    nprob = len(NEMOH_DIR) + sum(NEMOH_DOF)
    iprob = 2
    nfreq = len(w)
    problem = (ifreq - 1) * nprob + iprob

    # print('Problem: {}'.format(problem))

    if 0:
        # hdp.show_pressure(ifreq, idir, pressure_type='Froude-Krylof',axis=2)
        # hdp.show_pressure(ifreq, idir, pressure_type='Diffraction',axis=2)
        hdp.show_pressure(ifreq, irad, pressure_type='Radiation', axis=2)

    # print(hdp.ma[0,:,0,0])
    # print(hdp._fe_amp[0,:,dof])

    idof = np.array([i for i, x in enumerate(NEMOH_DOF) if x])

    for i, x in enumerate(idof):
        w2 = hdp.k[x, x] / (hdp.m[x, x] + hdp.ma[ifreq, x, x])
        print('T{}:\t{:5.1f} s'.format(x + 1, 2 * np.pi / np.sqrt(w2)))

    np.set_printoptions(precision=3)
    idof = 2
    rao=hdp.getRAO(idof, idir)
    #print(rao)
    plt.plot(2 * np.pi / w, rao)
    plt.show()

    np.set_printoptions(precision=3)

    if 0:
        print('\n')




        print('\nRadiation damping:\n{}'.format(hdp._c_hyd[ifreq]))
        print('\nWater plane stiffness:\n{}'.format(hdp.k))
        print('\nStatic mass:\n{}'.format(hdp.m))
        print('\nAdded mass:\n{}'.format(hdp._ma[ifreq]))
        print('\nExcitation force:\n{}'.format(np.abs(hdp._fe[idir, ifreq, :])))
        # tmp = nemoh.get_section_forces(fio, problem, [-100,0,0], [1,0,0], sym=SYM)
        # print('\nSection force:\n{}'.format(np.abs(tmp)))

        f_fk = nemoh.p2f(hdp._p['Froude-Krylof'][ifreq, idir, :], hdp.pd)
        f_diff = nemoh.p2f(hdp._p['Diffraction'][ifreq, idir, :], hdp.pd)
        f_exc = f_fk + f_diff


        print(np.abs(nemoh.get_section_forces(f_exc, hdp.pd.ppanel_centers, [10, 0, 0], [1, 0, 0])))

        print(abs(hdp.p2f(ifreq, pressure_index=0, pressure_type='Hydro static')) / 9.81)
        #print(abs(sum(nemoh.p2f(hdp._p['Hydro static'], hdp.pd))) / 9.81)

    part_list = unit_model.get_parts()
    print(sum([part.mass for part in part_list]))
