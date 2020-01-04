import pickle
import numpy as np

if __name__ == '__main__':
    with open( "this_candidate.pkl", "rb" ) as f:
        this_candidate=pickle.load(f)

    this_candidate.settings.load_cases = {  # 121, np.pi / 15, np.pi
            "num_wave_frequencies": 120,  # TODO: Implement adaptive frequency
            "min_wave_frequencies": 2 * np.pi / 35,  # (rad/s)
            "max_wave_frequencies": 2 * np.pi / 1,
            "num_wave_directions" : 1,
            "min_wave_directions" : 0,  # deg
            "max_wave_directions" : 90,
    }


    this_candidate.hydrodynamic_analysis()
