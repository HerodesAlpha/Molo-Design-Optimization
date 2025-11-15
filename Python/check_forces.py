"""
Check forces script for analyzing and visualizing section forces and RAOs.

This script loads a design candidate, sets up environmental conditions,
and plots dynamic forces, RAOs, and static forces for analysis.
"""

import pickle

import matplotlib.pyplot as plt
import numpy as np

import core.environmental_conditions as ec


if __name__ == '__main__':
    with open("this_candidate.pkl", "rb") as f:
        this_candidate = pickle.load(f)
    print(f'\nCase:\t{this_candidate.settings.case_label}')

    yr = this_candidate.settings.park_data['Design Basis']['ULS']['Return period']

    # Create list of short terms from contour line
    area = this_candidate.settings.park_data['Design Basis']['Area']
    this_candidate.ltwc1 = ec.Long_Term_Wave_Conditions(area=area)
    this_candidate.cl = this_candidate.ltwc1.contour_line(yr)
    this_candidate.contourline = []
    # this_candidate.cl = np.asarray([[3.5, 13.5, 1],[9, 9.5, 5],[8, 11, 3.6],[7, 11, 2.6],[7, 11.5, 2.12]])

    for hs, tz in this_candidate.cl:
        this_candidate.contourline.append(ec.Short_Term_Wave_Conditions(hs=hs, tz=tz))

    this_candidate.settings.radiation_damping_factor = 1

    this_candidate.settings.do_linearize = True
    this_candidate.init_load()
    stwc = ec.Short_Term_Wave_Conditions(hs=9.5, tz=7.3)
    this_candidate.init_response(short_term_wave_condition=stwc)

    # Get section forces
    sp_x = this_candidate.settings.floater_data['Central column diameter'] * (0.5)  # + 0.8 + 1 + 0.8 + 1)
    sp_z = this_candidate.settings.floater_data['Radial']['Heigth'] / 2 - this_candidate.hs_floater.hs_data[
        'draught']
    sp_x = -100
    sp_z = 0
    section_point = [sp_x, 0, sp_z]  # Used for moment reference
    section_normal = [1, 0, 0]

    imass, ipanel, idamp = this_candidate.response.get_section_index(section_point, section_normal)
    this_candidate.f_sec1 = this_candidate.response.assemble_forces(imass, ipanel, idamp, moment_ref_point=[0, 0, 0])
    del imass, ipanel

    ibeta = 0
    idof = 2

    stat_force = this_candidate.f_sec1['Static']

    factor = 1 / 1000000

    print('\n--------------------------\n S T A T I C   F O R C E\n--------------------------')
    for key in stat_force:
        a = np.abs(stat_force[key][2]) * factor
        b = np.angle(stat_force[key][2])
        print(f'{key:20} {a:5.2f} {b: 5.2f}')

    w = this_candidate.loads.w
    fig, axs = plt.subplots(2, 2)
    dyn_force = this_candidate.f_sec1['Dynamic']

    x_tics = 1 / (2 * np.pi / w)
    x_label = 'Frequency [Hz]'

    # --------------
    # FREQUENCY
    # --------------
    ifreq_print = 30

    print('\n--------------------------\n D Y N A M I C   F O R C E\n--------------------------')
    print(f'{x_label:16} {x_tics[ifreq_print]:5.3f}\n')
    all_keys = ['Froude-Krylof', 'Diffraction', 'Mass', 'Added mass', 'Radiation damping', 'Viscous damping', 'Buoyancy', 'SUM']
    plot_keys = [all_keys[i] for i in [0, 1, 2, 3, 4, 5, 6, 7]]
    for key in dyn_force:
        if key in plot_keys:
            if key == 'SUM':
                linewidth = 2
            else:
                linewidth = 1
            abs_val = np.abs(dyn_force[key][:, ibeta, idof]) * factor

            phase_val = np.angle(dyn_force[key][:, ibeta, idof])
            axs[0, 0].plot(x_tics, abs_val, label=key, linewidth=linewidth)
            axs[0, 1].plot(x_tics, phase_val, label=key)
            # print(abs_val[:])
            print(f'{key:20} {abs_val[ifreq_print]:6.2f} {phase_val[ifreq_print]: 5.2f}')

    axs[0, 0].set_title('Amplitude')
    force_label = ['Fx [MN]', 'Fy [MN]', 'Fz [MN]', 'Mx [MNm]', 'My [MNm]', 'Mz [MNm]']
    axs[0, 0].set_ylabel(force_label[idof])
    # axs[0, 0].set_yscale('log')
    # axs[0, 0].set_ylim([0, 1])

    axs[0, 0].legend()
    axs[0, 0].grid()
    axs[0, 1].set_title('Phase')
    axs[0, 1].legend()
    axs[0, 1].grid()
    axs[1, 0].set_ylabel('RAO')
    axs[1, 0].set_xlabel(x_label)
    axs[1, 0].grid()
    axs[1, 1].set_xlabel(x_label)
    axs[1, 1].grid()
    print('\n--------------------------\n R A O\n--------------------------')
    d = {'Heave': 2, 'Pitch': 4}
    ax1 = axs[1, 0]
    ax2 = ax1.twinx()

    for r in [this_candidate.response.rao, this_candidate.response.rao_init]:

        for key in d:

            if key == 'Heave':
                abs_val_rao = np.abs(r[:, ibeta, d[key]])
                lns1 = ax1.plot(x_tics, abs_val_rao, label=key)
                phase_val_rao = np.angle(r[:, ibeta, d[key]])
                print(f'{key:20} {abs_val_rao[ifreq_print]:6.3f} {phase_val_rao[ifreq_print]: 7.4f}')

            else:
                abs_val_rao = np.abs(r[:, ibeta, d[key]])
                lns2 = ax2.plot(x_tics, abs_val_rao, '-r', label=key)

                phase_val_rao = np.angle(r[:, ibeta, d[key]])
                print(f'{key:20} {abs_val_rao[ifreq_print]:6.3f} {phase_val_rao[ifreq_print]: 7.4f} ({abs_val_rao[ifreq_print] * 180 / np.pi: 5.1f} deg)')

            axs[1, 1].plot(x_tics, phase_val_rao, label=key)

    lns = lns1 + lns2
    labs = [l.get_label() for l in lns]
    axs[1, 0].legend(lns, labs, loc=0)
    axs[1, 1].legend()
    plt.suptitle(f'{this_candidate.settings.case_label}\nSection point: [{section_point[0]:1.2f}, {section_point[1]:1.2f}, {section_point[2]:1.2f}]\nSection normal: [{section_normal[0]:1.2f}, {section_normal[1]:1.2f}, {section_normal[2]:1.2f}]')
    plt.show()

    # f_rad=np.sum(this_candidate.loads._force['Radiation'][ifreq_print,idof,:,2])
    # A=np.imag(f_rad)/w[ifreq_print]
    # B=-np.real(f_rad)
    # print('\n{:15} {:5.0f} {:5.0f}'.format('A and B', A, B))
    # z = A*np.exp(complex(0, 1)*B)
    # print('\n{:15} {:5.0f}'.format('Z', z))

    # a=this_candidate.response._point_mass[:, 2, 2]
    # b=this_candidate.response._point_mass_centers
    #
    # print()
    # print(sum(a)
    # print(this_candidate.loads._m[2,2])
    #
    # I=np.sum(a[:,np.newaxis]*b**2,axis=0)
    # print(I)
    # print(np.diag(this_candidate.loads._m[3:,3:]))
    np.set_printoptions(precision=3)
    #
    # print()
    # print(this_candidate.loads._m)
    # print()
    # print(this_candidate.loads._ma[ifreq_print,:,:])
    # print()
    # print(this_candidate.loads._k)
    # print()
    # print(np.diag(this_candidate.response._c_visc[ifreq_print,:]))
