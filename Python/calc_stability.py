"""
Calculate stability script for evaluating design candidate stability.

This script loads a design candidate and calculates its intact stability ratio.
"""

import pickle
from multiprocessing import freeze_support


if __name__ == '__main__':
    # freeze_support()

    print('This is main')
    with open("this_candidate.pkl", "rb") as f:
        this_candidate = pickle.load(f)

    this_candidate.settings.create_stability_movie = True

    r = this_candidate.intact_stability_ratio()
    if r >= 1.4:
        is_stable = True
        print(f'This candidate is stable with r = {r * 100:1.0f}%')
    else:
        is_stable = False
        print(f'This candidate is NOT stable with r = {r * 100:1.0f}%')
