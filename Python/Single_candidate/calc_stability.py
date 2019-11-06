import pickle

with open( "this_candidate.pkl", "rb" ) as f:
    this_candidate=pickle.load(f)
r = this_candidate.intact_stability_ratio()
if r >= 1.4:
    is_stable = True
    print('This candidate is stable with r = {:1.0f}%'.format(r * 100))
else:
    is_stable = False
    print('This candidate is NOT stable with r = {:1.0f}%'.format(r * 100))
