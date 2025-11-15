"""
Single run script for automated design candidate analysis.

This script performs a complete analysis workflow including model creation,
stability check, hydrodynamic analysis, and structural optimization.
"""

import os
import sys
import warnings
from contextlib import redirect_stdout
from pathlib import Path

import numpy as np
import pandas

import core.tool_box as tb
from core.common import BreakIt
from optimization.design_engine import Candidate, Parameter_Space


if __name__ == '__main__':
    print_to_screen = True
    df = pandas.DataFrame(columns=('MODEL',
                                   'Height of lower flange stiffener',
                                   'Thickness of lower flange',
                                   'Thickness of lower flange stiffener',
                                   'SFX',
                                   'SFY',
                                   'SFZ',
                                   'SMX',
                                   'SMY',
                                   'SMZ',
                                   'PFX',
                                   'PFY',
                                   'PFZ',
                                   'PMX',
                                   'PMY',
                                   'PMZ'))

    if not sys.warnoptions:
        warnings.simplefilter("default")

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
    p.radial_column_diameter = 7.5
    p.filling_ratio = [0.1 , 0.1]
    this_candidate = Candidate(analyses_root, park_label, wtg_label, p, case_label_type='molo_model')

    this_candidate.settings.wtg_model = "Vestas 9.5"
    #this_candidate.settings.wtg_model = "Haliade X"


    with open(this_candidate.settings.fio.case_dir.joinpath('stdout_redirect.txt'), 'w') as fout:

        if not this_candidate.has_model():
            state = 'New'
        else:
            state = "Old"

        state = 'New'
        this_candidate.init_model(state=state)
        print(f'\nCase:\t{this_candidate.settings.case_label}')
        # this_candidate.init_model(state='New')
        this_candidate.settings.load_cases = {  # 121, np.pi / 15, np.pi
                "num_wave_frequencies": 80,  # TODO: Implement adaptive frequency
                "min_wave_frequencies": 2 * np.pi / 35,  # (rad/s)
                "max_wave_frequencies": 2 * np.pi / 4,
                "num_wave_directions" : 1,
                "min_wave_directions" : 0,  # deg
                "80" : 90,
        }

        if not this_candidate.has_stability_db() or 1:
            # if True:
            #try:
            r = this_candidate.intact_stability_ratio()
            if r >= 1.4:
                is_stable = True
                print('This candidate is stable with r = {:1.0f}%'.format(r * 100))
            else:
                is_stable = False
                print('This candidate is NOT stable with r = {:1.0f}%'.format(r * 100))
            # except:
            #     print('Stability check failed ...')
            #     is_stable = False
        else:
            print('This candidate has old stability database')

        #exit()
        if this_candidate.has_complete_hydrodynamic_db() and 0:
            print('This candidate has old hydro database. Do not perform hydrodynamic analysis')


        else:
            if this_candidate.is_stable() or 1:
                print(
                        'This candidate is stable, but has no hydro_database. Perform hydrodynamic analysis')
                this_candidate.hydrodynamic_analysis()

        if this_candidate.is_stable() or 1:
            if this_candidate.has_complete_hydrodynamic_db():
                # if False:
                print(
                        'This candidate is stable and has hydro database. Perform structural optimization')

                this_candidate.settings.critical_damping_ratio = 0.02
                this_candidate.init_load_response()

                irow += 1
                res = this_candidate.structural_analysis()

                if print_to_screen:
                    print(f"Max UR is : {res['Max UR']:1.2f}")
                    print(f"Height of lower flange stiffener : {res['panel_cc']._h_lfst:1.2f}")
                    print(f"Thickness of lower flange : {res['panel_cc']._t_lfst:1.2f}")
                    print(f"Thickness of lower flange stiffener : {res['panel_cc']._t_lf:1.2f}")

                    print(f"Section force: {np.array2string(res['section force'], precision=2)}")
                    print(f"Panel force: {np.array2string(res['panel force'], precision=2)}")
                    # this_candidate.loads.show_pressure(ifreq=20, pressure_index=1, pressure_type='Radiation', axis=2)

                    print('\nEigenvalue sollution WITH added mass')
                    tb.eigenvalprint(this_candidate.loads.m + this_candidate.loads.ma[0, :, :],
                                     this_candidate.loads.k)

                    tb.matprint(this_candidate.loads.m + this_candidate.loads.ma[0, :, :])
                    tb.matprint(this_candidate.loads.k)

                df.loc[irow] = [this_candidate.settings.case_label,
                                res['panel_cc']._h_lfst,
                                res['panel_cc']._t_lfst,
                                res['panel_cc']._t_lf,
                                res['section force'][0],
                                res['section force'][1],
                                res['section force'][2],
                                res['section force'][3],
                                res['section force'][4],
                                res['section force'][5],
                                res['panel force'][0],
                                res['panel force'][1],
                                res['panel force'][2],
                                res['panel force'][3],
                                res['panel force'][4],
                                res['panel force'][5],
                                ]


                print('optimization finished, results saved and report printed')

            else:
                print('This candidate is stable but has no hydro_database. No structural analysis')
        else:
            print('This candidate is not stable')

    this_candidate.print_report()
    df.to_csv(candidate_parent_dir.joinpath('output.csv'))
    df.to_excel(candidate_parent_dir.joinpath('output.xlsx'))
