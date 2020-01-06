import pickle
import numpy as np
import matplotlib.pyplot as plt
import core.environmental_conditions as ec

with open("this_candidate.pkl", "rb") as f:
    this_candidate = pickle.load(f)
print('\nCase:\t{}'.format(this_candidate.settings.case_label))

force_label = ['Fx [MN]', 'Fy [MN]', 'Fz [MN]', 'Mx [MNm]', 'My [MNm]', 'Mz [MNm]']

yr = this_candidate.settings.park_data['Design Basis']['ULS']['Return period']

# Create list of short terms from contour line
#area = this_candidate.settings.park_data['Design Basis']['Area']
area = 9
this_candidate.ltwc1 = ec.Long_Term_Wave_Conditions(area=area)
this_candidate.cl = this_candidate.ltwc1.contour_line(27)
this_candidate.contourline = []
#this_candidate.cl = np.asarray([[3.5, 13.5, 1],[9, 9.5, 5],[8, 11, 3.6],[7, 11, 2.6],[7, 11.5, 2.12]])

for hs, tz, in this_candidate.cl:
    print(' {:6.1f} {:6.1f}'.format(hs, tz))
    this_candidate.contourline.append(ec.Short_Term_Wave_Conditions(hs=hs, tz=tz, gamma=None))

this_candidate.settings.radiaton_damping_factor = 1



this_candidate.settings.do_linearize=True
this_candidate.init_load()
ibeta = 0
for stwc in this_candidate.contourline:
    this_candidate.init_response(short_term_wave_condition=stwc)

    # Get section forces
    sp_x = this_candidate.settings.floater_data['Central column diameter'] * (0.5)  # + 0.8 + 1 + 0.8 + 1)
    sp_z = this_candidate.settings.floater_data['Radial']['Heigth'] / 2 - this_candidate.hs_floater.hs_data[
        'draught']
    # sp_x = -100
    sp_z = 0
    section_point = [sp_x, 0, sp_z]  # Used for moment reference
    section_normal = [1, 0, 0]

    imass, ipanel, istrip = this_candidate.response.get_section_index(section_point, section_normal)
    this_candidate.f_sec1 = this_candidate.response.assemble_forces(imass, ipanel, istrip, moment_ref_point=[0, 0, 0])
    del imass, ipanel




    fz = this_candidate.f_sec1['Dynamic']['SUM'][:, ibeta, 2] / 1000000
    my = this_candidate.f_sec1['Dynamic']['SUM'][:, ibeta, 4] / 1000000



    fz_elm=stwc.expected_largest_maximum(fz, this_candidate.loads.w)
    my_elm=stwc.expected_largest_maximum(my, this_candidate.loads.w)


    print(' {:6.1f} {:6.1f} {:6.1f}  {:6.2f}  {:6.1f} '.format(stwc.hs, stwc.tp, stwc.gamma, fz_elm, my_elm))



