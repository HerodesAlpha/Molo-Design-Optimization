"""
Auto sizer script for automated design optimization.

This script performs parameter sweep optimization across multiple
design variables (diameter, ballast, height, gap) and generates
results in CSV and Excel formats.
"""

import sys
import warnings
from contextlib import redirect_stdout
from pathlib import Path

import numpy as np
import pandas

import core.tool_box as tb
from optimization.design_engine import Candidate, Parameter_Space


class BreakIt(Exception):
    """Custom exception for breaking out of nested loops."""
    pass


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

    analyses_root = Path(r'C:\MOLO_Optimization')
    park_label = 'site_01'
    wtg_label = 'wtg_01'
    candidate_parent_dir = analyses_root.joinpath(park_label).joinpath(wtg_label)

    p = Parameter_Space()

    p.ncol = 3

    is_stable = False
    #    for d in np.linspace(8.8, 8.8, 1, dtype=float):
    irow = -1
    for d in np.linspace(8.1, 8.1, 1, dtype=float):
        # print('d = {:1.2f}'.format(d))
        # if d < 8.4: continue
        for b in np.linspace(0, 0.2, 3, dtype=float):
            # print('b = {:1.2f}'.format(b))
            # if b < 0.2: continue
            for h in np.linspace(15, 30, 4, dtype=float):
                # print('h = {:1.2f}'.format(h))
                # if h < 30: continue
                for g in np.linspace(0.6, 1, 3, dtype=float):
                    # if g < 0.8: continue

                    p.gap = g
                    p.height = h
                    p.radial_column_diameter = d
                    p.filling_ratio = [b] * 3
                    this_candidate = Candidate(analyses_root, park_label, wtg_label, p, case_label_type='molo_model')
                    with open(this_candidate.settings.fio.case_dir.joinpath('stdout_redirect.txt'), 'w') as fout:

                        if not this_candidate.has_model():
                            state = 'New'
                        else:
                            state = "Old"
                        with redirect_stdout(fout):
                            this_candidate.init_model(state=state)
                        print(f'\nCase:\t{this_candidate.settings.case_label}')
                        # this_candidate.init_model(state='New')
                        this_candidate.settings.load_cases = {  # 121, np.pi / 15, np.pi
                                "num_wave_frequencies": 40,  # TODO: Implement adaptive frequency
                                "min_wave_frequencies": 2 * np.pi / 30,  # (rad/s)
                                "max_wave_frequencies": 2 * np.pi / 4,
                                "num_wave_directions" : 3,
                                "min_wave_directions" : 0,  # deg
                                "max_wave_directions" : 90,
                        }

                        if not this_candidate.has_stability_db():
                            # if True:
                            try:
                                with redirect_stdout(fout):
                                    r = this_candidate.intact_stability_ratio()
                                if r >= 1.4:
                                    is_stable = True
                                    print(f'This candidate is stable with r = {r * 100:1.0f}%')
                                else:
                                    is_stable = False
                                    print('This candidate is not stable')
                            except Exception as e:
                                print(f'Stability check failed: {e}')
                                is_stable = False
                        else:
                            print('This candidate has old stability database')

                        if this_candidate.has_complete_hydrodynamic_db():
                            print('This candidate has old hydro database. Do not perform hydrodynamic analysis')


                        else:
                            if this_candidate.is_stable():
                                print(
                                        'This candidate is stable, but has no hydro_database. Perform hydrodynamic analysis')
                                with redirect_stdout(fout):
                                    # pass
                                    this_candidate.hydrodynamic_analysis()

                        if this_candidate.is_stable():
                            if this_candidate.has_complete_hydrodynamic_db():
                                # if False:
                                print(
                                        'This candidate is stable and has hydro database. Perform structural optimization')

                                this_candidate.settings.critical_damping_ratio = 0.05
                                this_candidate.init_load_response(sea_spectrum=None)

                                irow += 1
                                with redirect_stdout(fout):
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

                                this_candidate.print_report()
                                print('optimization finished, results saved and report printed')

                            else:
                                print('This candidate is stable but has no hydro_database. No structural analysis')
                        else:
                            print('This candidate is not stable')

    df.to_csv(candidate_parent_dir.joinpath('output.csv'))
    df.to_excel(candidate_parent_dir.joinpath('output.xlsx'))
