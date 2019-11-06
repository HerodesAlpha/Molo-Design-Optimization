import pickle
with open( "this_candidate.pkl", "rb" ) as f:
    this_candidate=pickle.load(f)

this_candidate.init_load_response()

with open( "this_candidate_with_load_response.pkl", "wb" ) as f:
    pickle.dump( this_candidate,  f)