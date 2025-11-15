"""
Check pressure script for visualizing pressure distributions.

This script loads a design candidate and visualizes pressure distributions
on the mesh for different pressure types and frequencies.
"""

import pickle


if __name__ == '__main__':
    with open("this_candidate.pkl", "rb") as f:
        this_candidate = pickle.load(f)

    # this_candidate.loads.show_pressure(ifreq=20, pressure_index=2, pressure_type='Radiation', axis=2)

    this_candidate.loads.show_pressure(ifreq=15, pressure_index=0, pressure_type='Froude-Krylof', axis=2)

    # this_candidate.loads.show_pressure(ifreq=20, pressure_index=0, pressure_type='Diffraction', axis=2)
