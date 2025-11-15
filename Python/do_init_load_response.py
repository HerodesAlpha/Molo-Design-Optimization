"""
Initialize load response script for computing loads and responses.

This script loads a design candidate, initializes load and response calculations,
and saves the result for later use.
"""

import pickle


if __name__ == '__main__':
    with open("this_candidate.pkl", "rb") as f:
        this_candidate = pickle.load(f)

    this_candidate.init_load_response()

    with open("this_candidate_with_load_response.pkl", "wb") as f:
        pickle.dump(this_candidate, f)
