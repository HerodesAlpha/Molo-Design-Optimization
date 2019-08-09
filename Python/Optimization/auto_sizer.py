from pathlib import Path
import numpy as np
from design_engine import Candidate




def create_candidate(analyses_root, park_label, wtg_label,state=None):
    if state == None:
        state = 'New'
    c = Candidate(analyses_root, park_label, wtg_label)
    c.settings.case_label = 'molo_model'
    c.init_model(state='New')
    return c

if __name__ == '__main__':

    analyses_root = Path(r'C:\MOLO_Optimization')
    park_label = 'site_01'
    wtg_label = 'wtg_01'
    c1 = create_candidate(analyses_root, park_label, wtg_label,'Old')
    c1.settings.load_cases = {  # 121, np.pi / 15, np.pi
            "num_wave_frequencies": 41,
            "min_wave_frequencies": 2 * np.pi / 27,  # (rad/s)
            "max_wave_frequencies": 2 * np.pi / 4,
            "num_wave_directions" : 2,
            "min_wave_directions" : 0,  # deg
            "max_wave_directions" : 90,
    }
    #c1.run_intact_stability()
    #c1.run_hydrodynamic_analysis()
    c1.run_postprocessing()