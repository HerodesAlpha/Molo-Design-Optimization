import pickle
import numpy as np
import matplotlib.pyplot as plt
import core.environmental_conditions as ec

def get_axs(object):
    with open("this_candidate.pkl", "rb") as f:
        this_candidate = pickle.load(f)

    this_candidate.settings.do_linearize=True
    this_candidate.init_load()
    stwc=ec.Short_Term_Wave_Conditions(hs=9.5, tz=7.3)
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
    this_candidate.f_sec1 = this_candidate.response.assemble_forces(imass, ipanel, idamp, moment_ref_point=[0,0,0])
    del imass, ipanel

    ibeta = 0
    idof = 2

    stat_force = this_candidate.f_sec1['Static']

    factor = 1 / 1000000

    w = this_candidate.loads.w
    object.fig, object.axs = plt.subplots(2, 2)
    dyn_force = this_candidate.f_sec1['Dynamic']
    object_axs = object.axs
    x_tics = 1/(2 * np.pi / w)
    x_label = 'Frequency [Hz]'

    ifreq_print = 30

    all_keys=['Froude-Krylof','Diffraction','Mass', 'Added mass','Radiation damping','Viscous damping','Buoyancy','SUM']
    plot_keys = list( all_keys[i] for i in [0,1,2,3,4,5,6,7] )
    for key in dyn_force:
        if key in plot_keys:
            if key == 'SUM':
                linewidth=2
            else:
                linewidth = 1
            abs_val=np.abs(dyn_force[key][:, ibeta, idof]) * factor

            phase_val = np.angle(dyn_force[key][:, ibeta, idof])
            object_axs[0, 0].plot(x_tics,abs_val, label=key, linewidth=linewidth)
            object_axs[0, 1].plot(x_tics,phase_val , label=key)

    object_axs[0, 0].set_title('Amplitude')
    force_label = ['Fx [MN]', 'Fy [MN]', 'Fz [MN]', 'Mx [MNm]', 'My [MNm]', 'Mz [MNm]']
    object_axs[0, 0].set_ylabel(force_label[idof])

    object_axs[0, 0].legend()
    object_axs[0, 0].grid()
    object_axs[0, 1].set_title('Phase')
    object_axs[0, 1].legend()
    object_axs[0, 1].grid()
    object_axs[1, 0].set_ylabel('RAO')
    object_axs[1, 0].set_xlabel(x_label)
    object_axs[1, 0].grid()
    object_axs[1, 1].set_xlabel(x_label)
    object_axs[1, 1].grid()

    d = {'Heave': 2, 'Pitch': 4}
    ax1 = object_axs[1, 0]
    ax2 = ax1.twinx()
    lns1 = None
    lns2 = None
    for r in [this_candidate.response.rao, this_candidate.response.rao_init]:

        for key in d:

            if key == 'Heave':
                abs_val_rao = np.abs(r[:, ibeta, d[key]])
                lns1 = ax1.plot(x_tics, abs_val_rao, label=key)
                phase_val_rao = np.angle(r[:, ibeta, d[key]])


            else:
                abs_val_rao = np.abs(r[:, ibeta, d[key]])
                lns2 = ax2.plot(x_tics, abs_val_rao, '-r', label=key)

                phase_val_rao = np.angle(r[:, ibeta, d[key]])


            object_axs[1, 1].plot(x_tics, phase_val_rao, label=key)



    lns = lns1 + lns2
    labs = [l.get_label() for l in lns]
    object_axs[1, 0].legend(lns, labs, loc=0)
    object_axs[1, 1].legend()


