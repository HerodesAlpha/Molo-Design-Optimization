from pathlib import Path
import numpy as np
from design_engine import Candidate, Parameter_Space
import sys
import warnings
import tool_box as tb


if __name__ == '__main__':

    if not sys.warnoptions:
        warnings.simplefilter("default")

    analyses_root = Path(r'C:\MOLO_Optimization')
    park_label = 'site_01'
    wtg_label = 'wtg_01'
    p = Parameter_Space()
    p.height = 20
    p.ncol = 3
    p.gap = 0.8
#    for d in np.linspace(8.8, 8.8, 1, dtype=float):
    for d in np.linspace(8.1, 8.1, 1, dtype=float):
        p.column_diameter = d
        p.filling_ratio = [0.1]*3
        this_candidate = Candidate(analyses_root, park_label, wtg_label, p, case_label_type='molo_model')
        if not this_candidate.has_old_db():
            state = 'New'
        else:
            state = "Old"
        this_candidate.init_model(state=state)
        #this_candidate.init_model(state='New')
        this_candidate.settings.load_cases = {  # 121, np.pi / 15, np.pi
                "num_wave_frequencies": 60, # TODO: Implement adaptive frequency
                "min_wave_frequencies": 2 * np.pi / 35,  # (rad/s)
                "max_wave_frequencies": 2 * np.pi / 2,
                "num_wave_directions" : 1,
                "min_wave_directions" : 0,  # deg
                "max_wave_directions" : 0,
        }
        if not this_candidate.has_old_db():
            if this_candidate.intact_stability_ratio() >= 1.4:
                this_candidate.hydrodynamic_analysis()
            else:
                print('Intact stability ratio < 1.4')
                exit()

        res = this_candidate.structural_analysis()
        print('\nCase:\t{}'.format(this_candidate.settings.case_label))
        print('Max UR is : {:1.2f}'.format(res['Max UR']))
        print('Height of lower flange stiffener : {:1.2f}'.format(res['panel_cc']._h_lfst))
        print('Thickness of lower flange : {:1.2f}'.format(res['panel_cc']._t_lfst))
        print('Thickness of lower flange stiffener : {:1.2f}'.format(res['panel_cc']._t_lf))

        print('Setion force: {}'.format(np.array2string(res['section force'],precision=2)))

        print('Panel force: {}'.format(np.array2string(res['panel force'],precision=2)))
        #this_candidate.loads.show_pressure(ifreq=20, pressure_index=1, pressure_type='Radiation', axis=2)


        print('\nEigenvalue sollution WITH added mass')
        tb.eigenvalprint(this_candidate.loads.m + this_candidate.loads.ma[0, :, :], this_candidate.loads.k)

        tb.matprint(this_candidate.loads.m + this_candidate.loads.ma[0, :, :])
        tb.matprint(this_candidate.loads.k)



        #this_candidate.print_report()
