import pickle
import numpy as np
import matplotlib.pyplot as plt
import environmental_conditions as ec

with open("this_candidate.pkl", "rb") as f:
    this_candidate = pickle.load(f)
print('\nCase:\t{}'.format(this_candidate.settings.case_label))


yr = this_candidate.settings.park_data['Design Basis']['ULS']['Return period']

# Create list of short terms from contour line
area = this_candidate.settings.park_data['Design Basis']['Area']
this_candidate.ltwc1 = ec.Long_Term_Wave_Conditions(area=area)
this_candidate.cl = this_candidate.ltwc1.contour_line(yr)
this_candidate.contourline = []
#this_candidate.cl = np.asarray([[3.5, 13.5, 1],[9, 9.5, 5],[8, 11, 3.6],[7, 11, 2.6],[7, 11.5, 2.12]])

for hs, tz in this_candidate.cl:
    this_candidate.contourline.append(ec.Short_Term_Wave_Conditions(hs=hs, tz=tz))

this_candidate.settings.radiaton_damping_factor = 1



this_candidate.settings.do_linearize=True
this_candidate.init_load()
sea_spectrum=ec.Short_Term_Wave_Conditions(hs=6, tz=7).s_jonswap(this_candidate.loads.w)
this_candidate.init_response(sea_spectrum=sea_spectrum)

# Get section forces
sp_x = this_candidate.settings.floater_data['Central column diameter'] * (0.5)  # + 0.8 + 1 + 0.8 + 1)
sp_z = this_candidate.settings.floater_data['Radial']['Heigth'] / 2 - this_candidate.hs_floater.hs_data[
    'draught']
# sp_x = -100
sp_z = 0
section_point = [sp_x, 0, sp_z]  # Used for moment reference
section_normal = [1, 0, 0]

imass, ipanel = this_candidate.response.get_section_index(section_point, section_normal)
this_candidate.f_sec1 = this_candidate.response.assemble_forces(imass, ipanel, moment_ref_point=[0, 0, 0])
del imass, ipanel







def f_elm_contour(f, this_candidate):
    return np.array([stwcl.expected_largest_maximum(f, this_candidate.loads.w) for stwcl in this_candidate.contourline])

force_label = ['Fx [MN]', 'Fy [MN]', 'Fz [MN]', 'Mx [MNm]', 'My [MNm]', 'Mz [MNm]']
dyn_force = this_candidate.f_sec1['Dynamic']

ibeta=0

fz = f_elm_contour(dyn_force['SUM'][:, ibeta, 2] / 1000000, this_candidate)
my = f_elm_contour(dyn_force['SUM'][:, ibeta, 4] / 1000000, this_candidate)
for i, cl in enumerate(this_candidate.contourline):

    print(' {:6.1f} {:6.1f} {:6.1f}  {:6.2f}  {:6.1f} '.format(cl.hs, cl.tp, cl.gamma, fz[i], my[i]))



