__author__ = "Eivind Sonju"
__copyright__ = "Copyright (C) 2017-2019 Verbun AS. All rights reserved."
__version__ = "2.0"

# Use meshmagic eq. solver until somthing faster can be implemented
import core.meshmagick.hydrostatics as hs
import numpy as np
from core.meshmagick.MMviewer import MMViewer
import vtk
import imageio.v2 as imageio
import os
import matplotlib.pyplot as plt
from pylatex import Section, Figure, NoEscape, NewPage
import h5py
import multiprocessing
from multiprocessing import Process
from joblib import Parallel, delayed
import multiprocessing as mp
from multiprocessing import freeze_support
from itertools import repeat
from copy import deepcopy

# print("Number of cpu's : ", multiprocessing.cpu_count())

# Stiffness is linear for heel < 10 deg, param. excit. not likely

def calc_residual(hs_floater, thetay):
    thetax = 0
    rot_matrix = hs_floater.mesh.rotate([thetax, thetay, 0.])
    hs_floater._gravity_center = np.dot(rot_matrix, hs_floater._gravity_center)
    hs_floater._rotation = np.dot(rot_matrix, hs_floater._rotation)
    hs_floater.set_displacement(hs_floater.mass)
    hs_floater._reinit_clipper()
    return hs_floater.residual


def righting_moment_parallell(settings, hs_floater): # Out of comission

    # Not working well!!!
    hs_floater.verbose_off()
    dthetay = (np.pi / 180)
    righting_moments=[]
    heel_angles=[]
    int_max_deg=90
    if 1:

        # Step 1: Init multiprocessing.Pool()
        pool = mp.Pool(mp.cpu_count())

        # Step 2: `pool.apply` the `howmany_within_range()`
        floater_list = []
        heel_angles=[]
        for i in range(int_max_deg):
            floater_list.append(deepcopy(hs_floater))
            heel_angles.append(i*dthetay)
        residual= pool.starmap(calc_residual, zip(floater_list,heel_angles))
        del floater_list


        # Step 3: Don't forget to close
        #pool.join()
        pool.close()

    # n_cpu = multiprocessing.cpu_count()

    # curve = Parallel(n_jobs=n_cpu)(delayed(calc)(hs_floater, thetay) for thetay in range(90))

    for i in range(int_max_deg):
        print('{:7.1f} {val[0]:7.2f} {val[1]:7.2f} {val[2]:7.2f}'.format(heel_angles[i] * 180 / np.pi,
                                                                     val=-residual[i] / 1000000),
          flush=True)

    with open(settings.fio.stability_dir.joinpath('gz.txt'), 'w+') as f_gz:
        f_gz.write('{:7s} {:7s} {:7s} {:7s}\n'.format('theta', 'fz', 'my', 'mz'))

    righting_moments=[-row[2] for row in residual]

    return np.asarray([heel_angles, righting_moments])

def righting_moment_curve(settings, hs_floater):
    parallell = False
    if parallell:
        return righting_moment_parallell(settings, hs_floater)
    else:
        return righting_moment_movie(settings, hs_floater)


def righting_moment_movie(settings, hs_floater):
    # np.linalg.solve
    hs_floater.verbose_off()
    # Init
    thetax = thetay = 0
    dthetay = (np.pi / 180)
    dthetax = 0

    ifile = 0

    x = 0

    heel_angles = []
    righting_moments = []

    with imageio.get_writer(settings.fio.stability_dir.joinpath('stability.mp4'), mode='I') as writer:
        with open(settings.fio.stability_dir.joinpath('gz.txt'), 'w+') as f_gz:
            f_gz.write('{:7s} {:7s} {:7s} {:7s}\n'.format('theta', 'fz', 'my', 'mz'))
            while not x < 0:
                rot_matrix = hs_floater.mesh.rotate([thetax, dthetay, 0.])
                hs_floater._gravity_center = np.dot(rot_matrix, hs_floater._gravity_center)
                hs_floater._rotation = np.dot(rot_matrix, hs_floater._rotation)
                hs_floater.set_displacement(hs_floater.mass)
                hs_floater._reinit_clipper()
                # hs_floater._update_hydrostatic_properties()

                thetay += dthetay
                # print(thetay)
                x = -hs_floater.residual[2]
                f_gz.write('{:7.1f} {val[0]:7.2f} {val[1]:7.2f} {val[2]:7.2f}\n'.format(thetay * 180 / np.pi,
                                                                                        val=-hs_floater.residual / 1000000))
                heel_angles.append(thetay)
                righting_moments.append(-hs_floater.residual[2])

                if settings.create_stability_movie:
                    vtk_polydata = hs_floater.mesh._vtk_polydata()
                    hs_floater.viewer = MMViewer(use_interactor=False)
                    hs_floater.viewer.add_polydata(vtk_polydata)
                    hs_floater.viewer.plane_on()
                    corner_annotation = vtk.vtkCornerAnnotation()
                    corner_annotation.SetLinearFontScaleFactor(2)
                    corner_annotation.SetNonlinearFontScaleFactor(1)
                    corner_annotation.SetMaximumFontSize(20)
                    corner_annotation.SetText(3, '{:5.1f} deg {:6.1f} MNm'.format(thetay * 180 / np.pi, x / 1000000))
                    corner_annotation.GetTextProperty().SetColor(0., 0., 0.)
                    hs_floater.viewer.renderer.AddViewProp(corner_annotation)
                    hs_floater.viewer.render_window.SetOffScreenRendering(1)
                    # hs_floater.viewer.ShowWindowOff()
                    hs_floater.viewer.show_no_interactive()
                    ifile += 1

                    # fio.stability_dir.joinpath('stability.mp4')

                    filename = str(settings.fio.stability_dir.joinpath('gz_{:05d}.gif'.format(ifile)))
                    hs_floater.viewer.save_png(filename)
                    image = imageio.imread(filename)
                    writer.append_data(image)
                    os.remove(filename)

                    # hs_floater.render_window_interactor.GetRenderWindow().Finalize()
                    # hs_floater.render_window_interactor.TerminateApp()
                    hs_floater.viewer.finalize()

                print('{:7.1f} {val[0]:7.2f} {val[1]:7.2f} {val[2]:7.2f}'.format(thetay * 180 / np.pi,
                                                                                 val=-hs_floater.residual / 1000000),
                      flush=True)

    return np.asarray([heel_angles, righting_moments])


def wind_heeling_moment_curve(settings, heel_angles):
    m_xy_unit = 0.905  # Must be multiplied with hub height, rotor diameter, and wind speed squared. Calibrated against 116.5m hub height
    # 167 m rotor diameter and 70.7 m/s wind speed at hub height
    z = settings._job_data['wtg'][settings.wtg_model]['Hub height']
    rd = settings._job_data['wtg'][settings.wtg_model]['Rotor diameter']
    u = settings._job_data['design_basis']["Wind"]['ESS']['u']
    H = settings._job_data['design_basis']["Wind"]['ESS']['Reference height']
    z0 = settings._job_data['design_basis']["Wind"]['Surface friction coefficient']
    uh = u * (1 + np.log(z / H) / np.log(H / z0))

    m = m_xy_unit * z * rd * uh ** 2
    heeling_moments = np.zeros(heel_angles.shape[0])
    for i, angle in enumerate(heel_angles):
        heeling_moments[i] = m * np.cos(angle)

    return np.asarray([heel_angles, heeling_moments])


def intact_stability(settings, hs_floater):
    rmc = righting_moment_curve(settings, hs_floater)
    whm = wind_heeling_moment_curve(settings, rmc[0, :])

    # Find second intercept
    a = [i > j for i, j in zip(rmc[1, :], whm[1, :])]
    if not any(a):
        r = 0
        intercept = 0
    else:
        intercept = [i for i, x in enumerate(a) if x][-1]  # Index of last righting moment greater than heeling moment
        r = sum(rmc[1, :intercept]) / sum(whm[1, :intercept])
    if r < 1.4:
        print('Requirements for intact stability is NOT fulfilled')
    else:
        print('Requirements for intact stability is fulfilled')
    print('\tArea ratio is {: 7.1f}% (requirement is 140%)'.format(r * 100))

    width = r'1\textwidth'

    if 0:
        settings._report._doc.append(NewPage())
        with settings._report._doc.create(Section('Stability')) as stability_section:
            textstr = 'Area ratio is {:1.0f}%\nU_10min = {:1.1f} m/s\nz = {:1.0f} m'.format(r * 100,
                                                                                            settings._job_data[
                                                                                                'design_basis'][
                                                                                                "Wind"]['ESS']['u'],
                                                                                            settings._job_data[
                                                                                                'design_basis'][
                                                                                                "Wind"]['ESS'][
                                                                                                'Reference height'])
            stability_section.append(textstr)
            with stability_section.create(Figure(position='htbp')) as plot:
                fig = plt.figure(1, figsize=(8, 5))
                ax = fig.add_subplot(111)
                ax.plot(rmc[0, :] * 180 / np.pi, rmc[1, :] / 1000000, label='Righting moment')
                ax.plot(whm[0, :] * 180 / np.pi, whm[1, :] / 1000000, label='Heeling moment')
                ib = whm[0, intercept] * 180 / np.pi
                ax.annotate('Second intercept',
                            xy=(ib, whm[1, intercept]), xycoords='data',
                            xytext=(0.8, 0.5), textcoords='axes fraction',
                            arrowprops=dict(arrowstyle="->"))
                ax.set_xlabel('Angle of inclination [degrees]')
                ax.set_ylabel('Moment [MNm]')
                ax.legend()
                plt.grid(visible=True, which='major', color='#666666', linestyle='-')
                # Show the minor grid lines with very faint and almost transparent grey lines
                plt.minorticks_on()
                plt.grid(visible=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
                plot.add_plot(width=NoEscape(width))
                plot.add_caption('Intact Stability')
                plt.close()

    if 1:
        fig = plt.figure(1, figsize=(8, 5))
        ax = fig.add_subplot(111)
        ax.plot(rmc[0, :] * 180 / np.pi, rmc[1, :] / 1000000, label='Righting moment')
        ax.plot(whm[0, :] * 180 / np.pi, whm[1, :] / 1000000, label='Heeling moment')
        ib = whm[0, intercept] * 180 / np.pi
        ax.annotate('Second intercept',
                    xy=(ib, whm[1, intercept]), xycoords='data',
                    xytext=(0.8, 0.5), textcoords='axes fraction',
                    arrowprops=dict(arrowstyle="->"))
        ax.set_xlabel('Angle of inclination [degrees]')
        ax.set_ylabel('Moment [MNm]')
        ax.legend()
        plt.grid(visible=True, which='major', color='#666666', linestyle='-')
        # Show the minor grid lines with very faint and almost transparent grey lines
        plt.minorticks_on()
        plt.grid(visible=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
        # plot.add_plot(width=NoEscape(width))
        # plot.add_caption('Intact Stability')
        plt.show()


    # with h5py.File(settings.fio.stability_dir.joinpath('stability.hdf5'), "a") as hdf5_stability_db:
    #    hdf5_stability_db.create_dataset('intact_stability_area_ratio', data=r)

    return r
