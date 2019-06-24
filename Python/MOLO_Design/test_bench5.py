__author__ = "Eivind Sonju"
__copyright__ = "Copyright (C) 2017-2019 Verbun AS. All rights reserved."
__version__ = "2.0"

import logging
import multiprocessing
import pickle
from pathlib import Path
import scipy
import matplotlib.pyplot as plt
import numpy as np
from logutils.queue import QueueListener
import tool_box as tb
import calculations
import model_set_up as msu
import nemoh
import stability
from MOLO_Nemoh import nemoh_frontend as nf
from common import SettingsClass
from pyNemoh.structure import BaseStructure
import warnings
from meshmagick.mesh import Mesh
import force

if __name__ == '__main__':
    ANALYSES_ROOT = Path(r'C:\analyses')
    PARK_LABEL = 'site_01'
    WTG_LABEL = 'wtg_01'

    settings = SettingsClass(ANALYSES_ROOT, PARK_LABEL, WTG_LABEL)

    settings.create_model = True
    settings.calc_gz = True
    settings.run_nemoh = False
    settings.postprocessing = True

    h5_bs = BaseStructure()

    filling_ratio = [0.0, 0.0, 0.05]

    settings.floater_data = {
            "Type"                    : "OY",
            "Central column diameter" : 7.5,
            "Central column thickness": 0.04,
            "Draught"                 : 0,
            "Gap factor"              : 0.8,
            "Lower flange thickness"  : 0.04,
            "Number of radial columns": 3,
            "Radial column diameter"  : 7.5,
            "Radial column thickness" : 0.04,
            "Radial height"           : 15,
            "Upper flange thickness"  : 0.04,
            "Ballast filling ratio"   : [0,
                                         filling_ratio
                                         ],
            "Thin panel offset"       : 0.2
    }
    settings.load_cases = {  # 121, np.pi / 15, np.pi
            "num_wave_frequencies": 3,
            "min_wave_frequencies": 2 * np.pi / 30,  # (rad/s)
            "max_wave_frequencies": 2 * np.pi / 4,
            "num_wave_directions" : 2,
            "min_wave_directions" : 0,  # deg
            "max_wave_directions" : 90,
    }

    settings.case_label = 'floater_data'
    settings.set_file_structure()
    settings.simulation_dir = str(settings.fio.nemoh_root)
    settings.save_job_settings()

    settings.mesh_name = None
    settings.use_dipols = False
    settings.use_symmmetri = True

    settings.do_equilibriate = True
    fio = settings.fio

    settings.thin_panel_offset = 1.0

    if settings.create_model:
        print('\n--------------------------------------------------------------------------------------------')
        print('CREATE MODEL')
        print('--------------------------------------------------------------------------------------------')

        unit_model, hs_floater = msu.init_models(settings)
        pickle.dump(unit_model, open(fio.data_io_dir.joinpath('unit_model.pkl'), 'wb'))
        pickle.dump(hs_floater, open(fio.data_io_dir.joinpath('hs_floater.pkl'), 'wb'))

        settings.thin_panels = []

        m, k = tb.load_M_and_K(settings.fio.data_io_dir)
        # m=m[0:5,0:5]
        # k=k[0:5,0:5]
        # Set values close to zero to zero
        # for i in range(6):
        #     for j in range(6):
        #         if abs(m[i, j]) < 0:
        #             m[i, j] = 0
        #         if abs(k[i, j]) < 0:
        #             k[i, j] = 0
        # if (i == 0 and j == 0) or (i == 1 and j == 1) or (i == 5 and j == 5):
        #    k[i, j] = 1

        print('\nEigenvalue sollution WITHOUT added mass (given as lambda^0.5)')
        print('\nMass matrix:')
        tb.matprint(m)
        print('\nStiffness matrix:')
        tb.matprint(k)
        print('')
        tb.eigenvalprint(m, k)

    else:
        unit_model = pickle.load(open(settings.fio.data_io_dir.joinpath('unit_model.pkl'), 'rb'))
        hs_floater = pickle.load(open(settings.fio.data_io_dir.joinpath('hs_floater.pkl'), 'rb'))

    if settings.create_model and settings.calc_gz:
        print('\n--------------------------------------------------------------------------------------------')
        print('STABILITY ANALYSIS')
        print('--------------------------------------------------------------------------------------------')
        stability.gz_curve(settings, hs_floater)

    if settings.create_model and settings.run_nemoh:
        print('\n--------------------------------------------------------------------------------------------')
        print('NEMOH ANALYSIS')
        print('--------------------------------------------------------------------------------------------')
        settings.remove_old_db()
        queue = multiprocessing.Queue(-1)
        ql = QueueListener(queue, *logging.getLogger().handlers)
        ql.start()
        nf.run(settings._job_data['analysis'], queue)
        ql.stop()
        # nemoh.runNemoh(fio, hydro_mesh_symmetri, mesh_file, NEMOH_DIR, RHO_SW, WATER_DEPTH, OMEGA_NEMOH_INP,
        #                NEMOH_DOF)

    if settings.postprocessing:

        hdp = calculations.TransferFunctions(settings)

        SELECT_DOF = 1
        SELECT_AXIS = 2
        SELECT_FREC = 1

        # print('\nSection forces')
        # for i in range(6):
        #     print('\tDOF{}: {:8.1f}'.format(i + 1, np.sqrt(
        #         abs(sum(hydro.spec_response(hdp.f_sec[:, i][::3], hdp.w, hs=10, wp=2 * np.pi / 14))))))
        # write_report(root, hdp, case_label)

        ifreq = 0
        idir = 0
        irad = 4
        # nprob = len(NEMOH_DIR) + sum(NEMOH_DOF)
        iprob = 2
        # nfreq  = len(w)
        # problem = (ifreq - 1) * nprob + iprob

        # print('Problem: {}'.format(problem))
        NEMOH_DOF = [1, 1, 1, 1, 1, 1]

        if True:
            # hdp.show_pressure(ifreq, idir, pressure_type='Froude-Krylof',axis=2)
            # hdp.show_pressure(ifreq, idir, pressure_type='Diffraction', axis=2)
            # hdp.show_pressure(ifreq, irad, pressure_type='Radiation', axis=0)
            pass

        # print(hdp.ma[0,:,0,0])
        # print(hdp._fe_amp[0,:,dof])

        idof = np.array([i for i, x in enumerate(NEMOH_DOF) if x])

        print('Eigenvalue sollution WITH added mass')
        tb.eigenvalprint(hdp.m + hdp.ma[ifreq, :, :], hdp.k)

        if True:
            print('\n')

            print('\nRadiation damping:')
            tb.matprint(hdp.c_hyd[ifreq])
            print('\nWater plane stiffness:')
            tb.matprint(hdp.k)
            print('\nStatic mass [tonne]:')
            tb.matprint(hdp.m / 1000)
            print('\nAdded mass [tonne]:')
            tb.matprint(hdp.ma[ifreq] / 1000)
            print('\nExcitation force:')
            tb.matprint(np.abs(hdp.fe[ifreq, idir, :]))
            f_fk = hdp.p2f('Froude-Krylof', ifreq, idir)
            f_diff = hdp.p2f('Diffraction', ifreq, idir)
            f_exc = f_fk + f_diff

            print('\n')
            tb.matprint(np.abs(nemoh.get_section_values(f_exc, hdp.pd.ppanel_centers, [0, 0, 0], [1, 0, 0])))
        # np.set_printoptions(precision=3)

        idir = 0

        # --------------------------------------------------------------------------------------------------------------
        # PLOT RESULTS
        # --------------------------------------------------------------------------------------------------------------
        if True:
            h = hdp.get_rao(idir)
            fig, axs = plt.subplots(3, 2)
            w = 2 * np.pi / hdp.w
            axs[0, 0].plot(w, abs(h[:, 2]), 'tab:orange')
            axs[0, 0].set_title('Rao heave')
            axs[0, 1].plot(w, abs(h[:, 4]), 'tab:orange')
            axs[0, 1].set_title('Rao pitch')
            axs[1, 0].plot(w, hdp.ma[:, 2, 2] + hdp.m[2, 2], 'tab:green')
            axs[1, 0].set_title('m+ma heave')
            axs[1, 1].plot(w, hdp.ma[:, 4, 4] + hdp.m[4, 4], 'tab:green')
            axs[1, 1].set_title('m+ma pitch')
            axs[2, 0].plot(w, abs(hdp._fe[:, idir, 2]), 'tab:blue')
            axs[2, 0].set_title('fe heave')
            axs[2, 1].plot(w, abs(hdp._fe[:, idir, 4]), 'tab:blue')
            axs[2, 1].set_title('fe pitch')
            plt.show()
        #
        # plt.plot(2 * np.pi / hdp.w, abs(hdp.fe[:, idir, idof]))
        # plt.show()

        # np.set_printoptions(precision=3)

        # print(abs(sum(nemoh.p2f(hdp._p['Hydro static'], hdp.pd))) / 9.81)

        print('\n--------------------------------------------------------------------------------------------')
        print('SECTION FORCES')
        print('--------------------------------------------------------------------------------------------')

        section_point = [0, 0, 0]
        section_normal = [1, 0, 0]

        # Prepare RAO's for motion dependent response variables
        rao = hdp.get_rao(idir)

        # Collect all forces acting on the section
        if True:

            # Hydro static / Buoyancy
            f_bouyancy = nemoh.get_section_values(np.real(hdp.p2f('Hydro_static')), hdp.pd.ppanel_centers, section_point,
                                                  section_normal)

            print('\nBuoyancy force')
            tb.matprint(f_bouyancy)

            # Gravity
            part_list = unit_model.get_parts()
            point_mass_gravity_force = np.asarray([[0, 0, part.mass] for part in part_list]) * settings.grav
            point_mass_centers = np.asarray([-part.reduction_point for part in part_list])
            f_gravity = nemoh.get_section_values(point_mass_gravity_force, point_mass_centers, section_point,
                                                 section_normal)

            print('\nGravity force')
            tb.matprint(f_gravity)

            # Froude-Krylof and diffraction
            f_fk = np.zeros([hdp.nw, 6], dtype=complex)
            f_diff = np.zeros([hdp.nw, 6], dtype=complex)

            for ifreq in range(hdp.nw):
                f_fk[ifreq, :] = nemoh.get_section_values(hdp.p2f('Froude-Krylof', ifreq, idir), hdp.pd.ppanel_centers,
                                                          section_point,
                                                          section_normal)
                f_diff[ifreq, :] = nemoh.get_section_values(hdp.p2f('Diffraction', ifreq, idir), hdp.pd.ppanel_centers,
                                                            section_point,
                                                            section_normal)

            # Calculate radiation force transferfunctions R = H * eta
            f_rad = np.zeros([hdp.nw, 6], dtype=complex)

            # get_rao create complex motion at origin per freq in all dofs for given wave dir
            # p2f takes pressure and create global x,y,z force at center of each panel
            # rao_at_panel transform motion at origin to motion and panel_centers

            for ifreq in range(hdp.nw):
                for irad in range(6):
                    # Here the RAO for each DOF is multiplied with each RAO dependent panel force (x,y,z)
                    this_f_rad = hdp.p2f('Radiation', ifreq, irad) * rao[ifreq, irad]  # TODO: Check if correct
                    f_rad[ifreq, :] += nemoh.get_section_values(this_f_rad, hdp.pd.ppanel_centers,
                                                                section_point,
                                                                section_normal)

            f_varying_buoyancy = np.zeros([hdp.nw, 6], dtype=complex)
            f_inertia = np.zeros([hdp.nw, 6], dtype=complex)

            part_list = unit_model.get_parts()
            part_mass = np.asarray([part.mass for part in part_list])
            part_meanpos = np.asarray([-part.reduction_point for part in part_list])

            for ifreq in range(hdp.nw):

                # Rotate the panels according to RAO
                rao_rot_mat = tb.rotation_matrix(rao[ifreq, 3:6])  # rao_rot_mat is complex

                # Gen dynamic position of panels and calc hydro static pressure
                panel_pos = np.transpose(np.dot(rao_rot_mat, hdp.pd.ppanel_centers.T))
                panel_pos += rao[ifreq, 0:3]
                p_dz = (hdp._rho_sw * hdp._grav) * panel_pos[:, 2]
                f_dz = hdp.p2f(p_dz)
                f_varying_buoyancy[ifreq, :] = nemoh.get_section_values(f_dz,
                                                                        hdp.pd.ppanel_centers, section_point,
                                                                        section_normal)



                # Get dynamic acceleration of part masses and calc inertia force

                part_dynpos = np.transpose(np.dot(rao_rot_mat, part_meanpos.T))
                part_dynpos += rao[ifreq, 0:3]
                part_dynacc = part_dynpos * hdp.w[ifreq] ** 2
                part_inertia_force = part_dynacc * part_mass[:, np.newaxis]
                f_inertia[ifreq, :] = nemoh.get_section_values(part_inertia_force,
                                                               part_meanpos, section_point,
                                                               section_normal)






            f_tot_dyn = f_fk + f_diff + f_rad + f_varying_buoyancy + f_inertia

            w = 2 * np.pi / hdp.w

            fig, axs = plt.subplots(2, 2)
            heave = 2
            pitch = 4
            for i,dof in enumerate([heave,pitch]):
                axs[0, i].plot(w, abs(f_fk[:, dof]), 'tab:blue', label='Froude-Krylof')
                axs[0, i].plot(w, abs(f_diff[:, dof]), 'tab:green', label='Diffraction')
                axs[0, i].plot(w, abs(f_rad[:, dof]), 'tab:orange', label='Radiation')
                axs[0, i].set_title('Potential forces')
                axs[0,i].legend()

                #axs[1, i].plot(w, abs(f_bouyancy[:, dof]), 'tab:blue', label='Buoyancy')
                #axs[1, i].plot(w, abs(f_gravity[:, dof]), 'tab:green', label='Gravity')
                #axs[1, i].set_title('Mean forces')
                #axs[1,i].legend()


                axs[1, i].plot(w, abs(f_varying_buoyancy[:, dof]), 'tab:blue', label='Varying Buoyancy')
                axs[1, i].plot(w, abs(f_inertia[:, dof]), 'tab:green', label='Inertia')
                axs[1, i].set_title('Varying forces')
                axs[1,i].legend()


            #plt.plot(2 * np.pi / hdp.w, abs(f_fk[:, 2]), label='Froude-Krylof')
            #plt.plot(2 * np.pi / hdp.w, abs(f_diff[:, 2]), label='Diffraction')
            #plt.plot(2 * np.pi / hdp.w, abs(f_rad[:, 2]), label='Radiation')
            #plt.plot(2 * np.pi / hdp.w, abs(f_hydro_static_rao[:, 2]), label='Hydro pressure')
            #plt.plot(2 * np.pi / hdp.w, abs(f_inertia[:, 2]), label='Inertia')
            #plt.plot(2 * np.pi / hdp.w, abs(f_tot_dyn[:, 2]), label='Total')


            plt.show()



        f_varying_buoyancy = np.zeros([hdp.nw, 6], dtype=complex)







        # Calc velocity and acc (not needed yet)
        # panel_vel = panel_pos * 1j * hdp.w(ifreq)
        # panel_acc = panel_pos * hdp.w(ifreq)**2

        xyz = np.zeros([hdp._pd._npanel, 3])
        xyz[:, 2] = 1
        nemoh_mesh = Mesh(hdp._pd.ppoints, hdp._pd.ppanels)
        h = force.show_force(nemoh_mesh, hdp._pd.ppanel_centers, np.imag(panel_pos * xyz))
        #h.show()



        # print('\nSum of all parts:\t{:5.2f} tonne'.format(sum([part.mass for part in part_list]) / 1000))

        # print('\nFz\t{: 7.2f} MN'.format(f_gravity[2] / 1000000))
        # print('Mx\t{: 7.2f} MNm'.format(f_gravity[3] / 1000000))
        # print('my\t{: 7.2f} MNm'.format(f_gravity[4] / 1000000))

        # print(hdp.ma_zero)
        # print(hdp.ma_inf)

        # this_mesh = Mesh(hdp.pd.ppoints, hdp.pd.ppanels)
        # this_mesh.show()
