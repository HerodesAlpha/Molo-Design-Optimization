# Use meshmagic eq. solver until somthing faster can be implemented
import meshmagick.hydrostatics as hs
import numpy as np
from meshmagick.MMviewer import MMViewer
import vtk
import imageio
import os

# TODO: Evaluate Mathieu instability

def gz_curve(fio,hs_floater):
    #np.linalg.solve
    hs_floater.verbose_off()
    # Init
    thetax = thetay = 0
    dthetay = (np.pi / 180)/4
    dthetax = 0

    ifile=0

    x = 0

    with imageio.get_writer(fio.stability_dir.joinpath('stability.mp4'), mode='I') as writer:
        while not x < 0:



            rot_matrix = hs_floater.mesh.rotate([thetax, dthetay, 0.])
            hs_floater._gravity_center = np.dot(rot_matrix, hs_floater._gravity_center)
            hs_floater._rotation = np.dot(rot_matrix, hs_floater._rotation)
            hs_floater.set_displacement(hs_floater.mass)
            hs_floater._reinit_clipper()
            #hs_floater._update_hydrostatic_properties()

            thetay += dthetay
            x = -hs_floater.residual[2]

            if 1:
                vtk_polydata = hs_floater.mesh._vtk_polydata()
                hs_floater.viewer = MMViewer()
                hs_floater.viewer.add_polydata(vtk_polydata)
                hs_floater.viewer.plane_on()
                corner_annotation = vtk.vtkCornerAnnotation()
                corner_annotation.SetLinearFontScaleFactor(2)
                corner_annotation.SetNonlinearFontScaleFactor(1)
                corner_annotation.SetMaximumFontSize(20)
                corner_annotation.SetText(3, '{:5.1f} deg {:6.1f} MNm'.format(thetay*180/np.pi, x/1000000))
                corner_annotation.GetTextProperty().SetColor(0., 0., 0.)
                hs_floater.viewer.renderer.AddViewProp(corner_annotation)
                hs_floater.viewer.render_window.SetOffScreenRendering(1)
                #hs_floater.viewer.ShowWindowOff()
                hs_floater.viewer.show()
                ifile += 1
                filename = 'gz_{:05d}.gif'.format(ifile)
                hs_floater.viewer.save_png(filename)
                image = imageio.imread(filename)
                writer.append_data(image)
                os.remove(filename)

                #hs_floater.render_window_interactor.GetRenderWindow().Finalize()
                #hs_floater.render_window_interactor.TerminateApp()
                hs_floater.viewer.finalize()




            print('{:7.1f} {val[0]:7.2f} {val[1]:7.2f} {val[2]:7.2f}'.format(thetay*180/np.pi, val=-hs_floater.residual/1000000))












        #
        #
        # dM += hs_floater.S55*np.pi/1800 # Add moment giving one deg rotation based on current water plane stiffness
        # dF1 = hs.Force(point=(0, 0, -0.5), value=(-dM, 0, 0), name="F1")
        # dF2 = hs.Force(point=(0, 0, +0.5), value=(dM, 0, 0), name="F2")
        # hs_floater.add_force(dF1)
        # hs_floater.add_force(dF2)
        # hs_floater.equilibrate()
        # v = np.dot(hs_floater._rotation,[0,0,1])
        # angle=np.arctan2(v[2], v[0]) * 180 / np.pi
        # print('{:7.2f} {:7.2f}'.format(angle, dM/1000000))