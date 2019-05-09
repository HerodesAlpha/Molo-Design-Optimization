# from meshmagick.mesh import Mesh

import meshmagick.hydrostatics as hs
import meshmagick.mesh_clipper as mc
import meshmagick.mmio as mmio
import model_set_up
import tool_box as tb
# from meshmagick.mesh import Mesh
from meshmagick.mesh import *

# unit_model, floater_model, wtg_model = model_set_up.read_excel_model('MOLO.xlsx', floater_id=1, tower_id=1, rna_id=1)
unit_model, floater_model, wtg_model = model_set_up.read_json(r'templates\model_template.json')
tb.write_gmsh(root, floater_model)

vertices, faces = mmio.load_MSH('MOLO_{}c.msh'.format(floater_model.nc))
MOLO = Mesh(vertices, faces)
MOLO.merge_duplicates()
MOLO.heal_normals()
MOLO.heal_mesh()
MOLO.rotate_z(-np.pi / 2)  # IMPORTANT
print('Before equilibrium calc')
unit_model.print_vector_matrix_global()
hs_MOLO = hs.Hydrostatics(MOLO, verbose=True)
hs_MOLO.gravity = 9.81
hs_MOLO.rho_water = 1025.
hs_MOLO.mass = unit_model.mass / 1000  # Give mass in tons
print('Mass given to hydro is {:5.2f} t'.format(hs_MOLO.mass))
# M=0.1*unit_model.mass
# F1=hs.Force(point=(0,0,7+5),value=(M,0,0),name="F1")
# F2=hs.Force(point=(0,0,7-5),value=(-M,0,0),name="F2")
# hs_MOLO.add_force(F1)
# hs_MOLO.add_force(F2)
hs_MOLO.gravity_center = -unit_model.inertias.reduction_point
hs_MOLO.equilibrate()
# Update model with calculated draft
z_d = hs_MOLO.hs_data['draught']
unit_model.inertias.reduction_point = [0, 0, unit_model.inertias.reduction_point[2] + z_d]
print('\n\nAfter equilibrium calc giving {:5.2f} m draught'.format(z_d))
unit_model.print_vector_matrix_global()

K33 = hs_MOLO.hs_data['stiffness_matrix'][0, 0]
K44 = hs_MOLO.hs_data['stiffness_matrix'][1, 1]
M33 = unit_model.inertias.mass_matrix_global[2][2]
M44 = unit_model.inertias.mass_matrix_global[3][3]
T33 = 1 / (np.sqrt(K33 / M33) / (2 * np.pi))
T44 = 1 / (np.sqrt(K44 / M44) / (2 * np.pi))
print('\n\nT33 = {:>4.1f} s'.format(T33))
print('T44 = {:>4.1f} s'.format(T44))

# with open('hydrostatic_report.txt', 'w+') as f:
#    f.write(hs_MOLO.get_hydrostatic_report())
print(hs_MOLO.get_hydrostatic_report())

tb.save_M_and_K(M=unit_model.inertias.mass_matrix_global, MMK=hs_MOLO.hs_data['stiffness_matrix'])

MOLO_waterline = mc.MeshClipper(hs_MOLO.mesh)

hydro_mesh = Mesh(MOLO_waterline.lower_mesh.vertices, MOLO_waterline.lower_mesh.faces)
hydro_mesh.merge_duplicates()
hydro_mesh.heal_normals()
hydro_mesh.heal_mesh()
hydro_mesh.show()

MOLO_symmetri = hydro_mesh.de_symmetrize_xz()
MOLO_symmetri.merge_duplicates()
MOLO_symmetri.heal_normals()
MOLO_symmetri.heal_mesh()
# MOLO_symmetri.show()

print('Number of faces: {}'.format(MOLO_symmetri.nb_faces))
print('Number of vertices: {}'.format(MOLO_symmetri.nb_vertices))

mmio.write_MAR('MOLO_sym.dat', MOLO_symmetri.vertices, MOLO_symmetri.faces)
# mmio.write_MAR('MOLO.dat', hydro_mesh.vertices, hydro_mesh.faces)
