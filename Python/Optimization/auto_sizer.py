from pathlib import Path
import numpy as np
from design_engine import Candidate, Parameter_Space
import sys
import warnings

def create_candidate(analyses_root, park_label, wtg_label, p, state=None):
    if state == None:
        state = 'New'
    c = Candidate(analyses_root, park_label, wtg_label, p, case_label_type='molo_model')
    c.init_model(state=state)
    return c


if __name__ == '__main__':

    if not sys.warnoptions:
        warnings.simplefilter("default")

    analyses_root = Path(r'C:\MOLO_Optimization')
    park_label = 'site_01'
    wtg_label = 'wtg_01'
    p = Parameter_Space()
    p.height = 15
    p.ncol = 3
    p.gap = 0.8
    for d in np.arange(8.6, 8.8, 0.2, dtype=float):
        p.column_diameter = d
        state = 'Old'

        c1 = create_candidate(analyses_root, park_label, wtg_label, p, state)
        c1.settings.load_cases = {  # 121, np.pi / 15, np.pi
                "num_wave_frequencies": 41,
                "min_wave_frequencies": 2 * np.pi / 27,  # (rad/s)
                "max_wave_frequencies": 2 * np.pi / 4,
                "num_wave_directions" : 2,
                "min_wave_directions" : 0,  # deg
                "max_wave_directions" : 90,
        }
        if state == 'Old':
            res = c1.structural_analysis()
            # while ur > 1:
                # p.stiffener_height *= 1.1
                # c1 = create_candidate(analyses_root, park_label, wtg_label, p, state)
                # ur = c1.structural_analysis()['Max UR']
                # print('Max UR: {:1.2f}'.format(ur))
            print('\nMax UR is : {:1.2f}'.format(res['Max UR']))
            print('Height of lower flange stiffener : {:1.2f}'.format(res['panel_cc']._h_lfst))
            print('Thickness of lower flange : {:1.2f}'.format(res['panel_cc']._t_lfst))
            print('Thickness of lower flange stiffener : {:1.2f}'.format(res['panel_cc']._t_lf))


        else:
            if c1.intact_stability_ratio() >= 1.4:
                c1.hydrodynamic_analysis()
                results = c1.structural_analysis()
                print('Max UR: {:1.2f}'.format(results['Max UR']))
            else:
                print('Intact stability ratio < 1.4')

