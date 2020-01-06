from core.meshmagick import mmio
from core.meshmagick.mesh import Mesh
import numpy as np
import json
import core.tool_box as tb

class BreakIt(Exception): pass




def printv(string):
    if 1:
        print(string)


np.set_printoptions(precision=3)



vtol = 0.01
drc = 7
dcc = 6
gaf = 0.8
hgt = 7
xO = 0
yO = 0
zO = -7
ncol = 2

da = (1 + gaf) * drc
dtheta = 2 * np.pi / 3
theta = [i * dtheta for i in range(3)]
limit = drc / 2 * (1 + vtol)
vertices, faces = mmio.load_MSH('MOLO_2c_thin_small.msh')
#start_mesh = Mesh(vertices, faces)
#start_mesh.show()
pd = tb.PanelData(vertices, faces)

is_flange_element = []
dipol = []
flipped = []
#
for i in range(pd.npanels):
    xp = pd.ppanel_centers[i][0]
    yp = pd.ppanel_centers[i][1]
    xc = yc = 0

    if abs(pd.ppanel_normals[i][2]) == 1:  # Flange element
        printv('Process flange element {}'.format(i))
        is_flange_element.append(i)
        # Check normal and flip if positive up
        if pd.ppanel_normals[i][2] == 1:
            flipped.append(i)
            printv('   Flip normal')

            faces[i] = faces[i][::-1]
        # Next, find dipols, i.e. elements not inside the cylinders
        # Assume dipol if not found inside cylinders
        found_inside = False
        try:
            for irad in range(3):
                dxc = np.cos(theta[irad]) * da
                dyc = np.sin(theta[irad]) * da
                for icol in range(ncol):
                    xc = dxc * (icol + 1)
                    yc = dyc * (icol + 1)
                    if np.sqrt((xp - xc) ** 2 + (yp - yc) ** 2) < drc / 2:
                        found_inside = True
                        raise BreakIt
        except BreakIt:
            pass
        if np.sqrt((xp) ** 2 + (yp) ** 2) < dcc / 2:
            found_inside = True

        if not found_inside:
            printv('   Is dipol')
            # print(i)
            dipol.append(i)
    else:  # Cylinder element
        # Find the cylinder x,y to which the element belongs
        printv('Process cylinder element {}'.format(i))
        not_found = True
        try:
            for irad in range(3):
                dxc = np.cos(theta[irad]) * da
                dyc = np.sin(theta[irad]) * da
                for icol in range(ncol):
                    xc = dxc * (icol + 1)
                    yc = dyc * (icol + 1)
                    dx = np.abs(xp - xc)
                    dy = np.abs(yp - yc)
                    if dx < limit and dy < limit:
                        printv('   Found on radial {}, column {}'.format(irad + 1, icol + 1))
                        not_found = False
                        raise BreakIt
        except BreakIt:
            pass
        if not_found:
            xc = yc = 0 # Check center column
            if np.sqrt((xp) ** 2 + (yp) ** 2) < dcc / 2 * (1 + vtol):
                printv('   Found on center column'.format(i))
                not_found = False
        if not_found:
            printv('   Not found on any column'.format(i))
            exit()

            # Check if normal is outwards from cylinder center
        printv('   xp = {: 6.2f}, yp = {: 6.2f}'.format(xp, yp))
        printv('   xc = {: 6.2f}, yc = {: 6.2f}'.format(xc, yc))
        vec = np.asarray([xp - xc, yp - yc])
        pn = pd.ppanel_normals[i][0:2]
        if np.dot(vec, pn) < 0:  # dot product i positive for coordinates on the positive side of the plane
            printv('   Flip normal')
            faces[i] = faces[i][::-1]
            flipped.append(i)

print(faces[dipol])
model_template_path='test4.json'
with open(model_template_path, 'r') as f:
    data = json.loads(f.read())

# TODO: Update analysis.json in analysis with thin plates
data["simulations"]["default"]["calculation"]["thin_panels"] = dipol
#data["simulations"]["default"]["calculation"]["thin_panels"] = -1
with open(model_template_path, 'w') as f:
    f.write(json.dumps(data, indent=4, sort_keys=True))

full_mesh = Mesh(vertices, faces)
#dipol_mesh = Mesh(vertices, faces[dipol])
#flipped_mesh = Mesh(vertices, faces[flipped])
#full_mesh.show()
mmio.write_MAR('MOLO_2c.dat', full_mesh.vertices, full_mesh.faces)
