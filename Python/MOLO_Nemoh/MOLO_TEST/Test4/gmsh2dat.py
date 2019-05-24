from meshmagick import mmio
from meshmagick.mesh import Mesh
import numpy as np


class PanelData(object):
    def __init__(self, vertices, faces):
        self._ppoints = vertices
        self._ppanels = faces
        self._npoints = len(vertices)
        self._npanel = len(faces)
        _a1 = np.linalg.norm(np.cross(self._ppoints[self._ppanels[:, 1]] - self._ppoints[self._ppanels[:, 0]],
                                      self._ppoints[self._ppanels[:, 2]] - self._ppoints[self._ppanels[:, 0]]),
                             axis=1) * 0.5
        _a2 = np.linalg.norm(np.cross(self._ppoints[self._ppanels[:, 3]] - self._ppoints[self._ppanels[:, 0]],
                                      self._ppoints[self._ppanels[:, 2]] - self._ppoints[self._ppanels[:, 0]]),
                             axis=1) * 0.5
        self._ppanel_areas = _a1 + _a2

        self._c1 = np.sum(self._ppoints[self._ppanels[:, :3]], axis=1) / 3.
        self._c2 = (np.sum(self._ppoints[self._ppanels[:, 2:4]], axis=1) + self._ppoints[self._ppanels[:, 0]]) / 3.

        self._ppanel_centers = (np.array(([_a1, ] * 3)).T * self._c1 + np.array(([_a2, ] * 3)).T * self._c2)
        self._ppanel_centers /= np.array(([self._ppanel_areas, ] * 3)).T

        self._ppanel_cross = np.cross(self._ppoints[self._ppanels[:, 2]] - self._ppoints[self._ppanels[:, 0]],
                                      self._ppoints[self._ppanels[:, 3]] - self._ppoints[self._ppanels[:, 1]])
        self._norm = np.linalg.norm(self._ppanel_cross, axis=1)
        self._ppanel_normals = np.true_divide(self._ppanel_cross, self._norm[:, np.newaxis])

    @property
    def ppoints(self):
        return self._ppoints

    @property
    def ppanels(self):
        return self._ppanels

    @property
    def ppanel_normals(self):
        return self._ppanel_normals

    @property
    def ppanel_areas(self):
        return self._ppanel_areas

    @property
    def ppanel_centers(self):
        return self._ppanel_centers

    @property
    def npoints(self):
        return self._ppoints.shape[0]

    @property
    def npanels(self):
        return self._ppanels.shape[0]


np.set_printoptions(precision=3)
vtol = 0.1
drc = 7
dcc = 6
gaf = 0.8
ncol = 3
hgt = 7
xO = 0
yO = 0
zO = -7
ncol = 2

da = (1 + gaf) * drc
dtheta = 2 * np.pi / 3
theta = [i * dtheta for i in range(3)]

vertices, faces = mmio.load_MSH('MOLO_2c.msh')
pd = PanelData(vertices, faces)

is_flange_element = []
dipol = []
#
for i in range(pd.npanels):
    xp = pd.ppanel_centers[i][0]
    yp = pd.ppanel_centers[i][1]

    if abs(pd.ppanel_normals[i][2]) == 1:  # Flange element
        print('Process flange element {}'.format(i))
        is_flange_element.append(i)
        # Check normal and flip if positive up
        if pd.ppanel_normals[i][2] == 1:
            print('   Flip normal')

            faces[i] = faces[i][::-1]
        # Next, find dipols, i.e. elements not inside the cylinders
        # Assume dipol if not found inside cylinders
        found_inside = False
        for irad in range(3):
            dxc = np.cos(theta[irad]) * da
            dyc = np.sin(theta[irad]) * da
            for icol in range(ncol):
                xc = dxc * (icol + 1)
                yc = dyc * (icol + 1)
                if np.sqrt((xp - xc) ** 2 + (yp - yc) ** 2) < drc / 2:
                    found_inside = True
                    break
        if np.sqrt((xp) ** 2 + (yp) ** 2) < dcc / 2:
            found_inside = True

        if not found_inside:
            print('   Is dipol')
            # print(i)
            dipol.append(i)
    else:  # Cylinder element
        # Find the cylinder x,y to which the element belongs
        print('Process cylinder element {}'.format(i))
        not_found=True
        for irad in range(3):
            dxc = np.cos(theta[irad]) * da
            dyc = np.sin(theta[irad]) * da
            for icol in range(ncol):
                xc = dxc * (icol + 1)
                yc = dyc * (icol + 1)
                if np.sqrt((xp - xc) ** 2 + (yp - yc) ** 2) < drc / 2 * (1 + vtol):
                    print('   Found on radial {}, column {}'.format(irad+1,icol+1))
                    not_found = False
                    break
        if not_found:
            if np.sqrt((xp) ** 2 + (yp) ** 2) < dcc / 2 * (1 + vtol):
                print('   Found on center column'.format(i))
                not_found = False
                xc = yc = 0
        if not_found:
            print('   Not found on any column'.format(i))

        # Check if normal is outwards from cylinder center
        vec= pd.ppanel_centers[i] - np.asarray([xc,yc,pd.ppanel_centers[i][2]])
        dot = np.dot(vec, pd.ppanel_normals[i])
        #print(dot)
        if dot < 0:  # dot product i positive for coordinates on the positive side of the plane
            print('   Flip normal')
            faces[i] = faces[i][::-1]

print(dipol)

# print(pd.ppanel_centers[i])

# vec = pd.ppanel_centers - section_point
# dot = np.dot(vec, section_normal)  # dot product i positive for coordinates on the positive side of the plane
# ind = dot >= 0
# gsf = gpf[ind, :]  # Global Section Forces
# gsm = np.cross(ppanel_centers[ind], gsf)
# return sum(np.concatenate((gsf, gsm), axis=1))


start_mesh = Mesh(vertices, faces)
dipol_mesh = Mesh(vertices, faces[dipol])
start_mesh.show()
mmio.write_MAR('MOLO_2c.dat', start_mesh.vertices, start_mesh.faces)
