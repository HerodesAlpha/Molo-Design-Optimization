__author__ = "Eivind Sonju"
__copyright__ = "Copyright (C) 2017-2019 Verbun AS. All rights reserved."
__version__ = "2.0"

import json
import os

import numpy as np
import pandas as pd

import meshmagick.hydrostatics as hs
import meshmagick.mmio as mmio
import model
import tool_box as tb
from meshmagick.mesh import Mesh


def rem_list(df):
    # Dirty :-(
    for item in df:
        try:
            df[item] = df[item][0]
        except:
            pass
    return df


def read_excel_model(excel_file, floater_id, tower_id, rna_id):
    fdst = 'tmp.xlsx'
    # shutil.copy(excel_file,fdst)
    cmd = 'xcopy /r /y {} {}'.format(excel_file, fdst)
    try:
        os.system(cmd)
    except OSError:
        # File most likely does not exist
        print("File most likely does not exist")
        exit()

    df_floater = pd.read_excel(fdst, sheet_name='Floater', index=False)
    df_rna = pd.read_excel(fdst, sheet_name='RNA', index=False)
    df_tower = pd.read_excel(fdst, sheet_name='Tower', index=False)
    df_b2c = pd.read_excel(fdst, sheet_name='Ballast_2C', index=False)
    df_b3c = pd.read_excel(fdst, sheet_name='Ballast_3C', index=False)
    df_b4c = pd.read_excel(fdst, sheet_name='Ballast_4C', index=False)

    a = df_floater[df_floater.ID == floater_id].to_dict(orient='list')

    dict_tower_data = rem_list(df_tower[df_tower.ID == tower_id].to_dict(orient='list'))
    dict_floater_data = rem_list(df_floater[df_floater.ID == floater_id].to_dict(orient='list'))
    dict_rna_data = rem_list(rem_list(df_rna[df_rna.ID == rna_id].to_dict(orient='list')))

    twr_data = model.TowerDataClass(dict_tower_data)
    floater_data = model.FloaterDataClass(dict_floater_data)
    rna_data = model.RNADataClass(dict_rna_data)

    # Set global z
    interface_point = floater_data.t_lf + floater_data.hgt + floater_data.t_uf
    twr_data.p = [0, 0, -(interface_point + twr_data.h / 2)]
    rna_data.p = [0, 0, -(interface_point + twr_data.h)]
    # Write json
    # with open('model_template.json', 'w') as f:
    #     f.write(json.dumps({'twr': dict_tower_data, 'rna': dict_rna_data, 'floater': dict_floater_data,'test':[[1,0,1],[1,0,1]]}, indent=4,
    #                        sort_keys=True))

    rho_st = 7850
    #
    # Build floater
    floater_model = model.FloaterClass(floater_data, rho_st)  # CoGz will be set for BOS at z=0

    # Build tower and RNA
    wtg_model = model.WtgClass(twr_data, rna_data, rho_st)  # CoGz will be controlled by p

    # Assemble parts into complete unit
    unit_model = model.UnitClass([wtg_model, floater_model])
    return unit_model, floater_model, wtg_model


def model01():
    rna_mass = 500000
    turbine_diameter = 167
    hub_height_above_interface_point = turbine_diameter / 2 + 30

    fdi = dict()
    fdi['type'] = 'MOLO'
    fdi['nc'] = 4  # 2, 3 or 4
    fdi['gap'] = 0.8
    fdi['dia_rc'] = 7
    fdi['thi_rc'] = 0.05
    fdi['dia_hc'] = 6
    fdi['thi_hc'] = 0.04
    fdi['hgt'] = 14
    fdi['t_lf'] = 0.04
    fdi['t_uf'] = 0.04
    fdi['ballast'] = [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]
    fdi['rho_bal'] = 1025
    floater_data = model.FloaterDataClass(fdi)

    interface_point = floater_data.t_lf + floater_data.hgt + floater_data.t_uf

    tdi = dict()
    tdi['type'] = 'Tower'
    tdi['h'] = (hub_height_above_interface_point / 2 + 35)
    tdi['t'] = 0.05
    tdi['d'] = 6
    tdi['p'] = [0, 0, - (interface_point + hub_height_above_interface_point / 2)]

    rdi = dict()
    rdi['type'] = 'Vestas 9.5MW'
    rdi['d_rot'] = turbine_diameter
    rdi['m'] = rna_mass
    rdi['ixx'] = 5 ** 2 * rna_mass
    rdi['iyy'] = 10 ** 2 * rna_mass
    rdi['izz'] = 10 ** 2 * rna_mass
    rdi['iyz'] = 0
    rdi['ixz'] = 0
    rdi['ixy'] = 0
    rdi['p'] = [0, 0, - (interface_point + hub_height_above_interface_point)]

    rho_st = 7850
    #
    # Build floater
    floater_data = model.FloaterDataClass(fdi)
    floater_model = model.FloaterClass(floater_data, rho_st)  # CoGz will be set for BOS at z=0

    # Build tower and RNA
    twr_data = model.TowerDataClass(tdi)
    rna_data = model.RNADataClass(rdi)
    wtg_model = model.WtgClass(twr_data, rna_data, rho_st)  # CoGz will be controlled by p

    # Assemble parts into complete unit
    unit_model = model.UnitClass([wtg_model, floater_model])
    return unit_model, floater_model, wtg_model


def create_mass_models(settings):
    twr_data_found = False
    rna_data_found = False
    floater_data_found = False
    with open(json_file, 'r') as f:
        data = json.loads(f.read())
        for key in data:
            if key == 'Tower':
                if twr_data_found == True:
                    print('Tower data given more than once!')
                    exit()
                else:
                    dict_tower_data = data[key]
                    twr_data_found = True
            if key == 'Rotor-Nacelle-Assembly':
                if twr_data_found == True:
                    print('RNA data given more than once!')
                    exit()
                else:
                    dict_rna_data = data[key]
                    rna_data_found = True
            if key == 'Floater':
                if floater_data_found == True:
                    print('Floater data given more than once!')
                    exit()
                else:
                    dict_floater_data = data[key]
                    floater_data_found = True

    rho_st = 7850
    #
    twr_data = model.TowerDataClass(dict_tower_data)
    floater_data = model.FloaterDataClass(dict_floater_data)
    rna_data = model.RNADataClass(dict_rna_data)

    # Set global z
    interface_point = floater_data.t_lf + floater_data.hgt + floater_data.t_uf
    twr_data.p = [0, 0, -(interface_point + twr_data.h / 2)]
    rna_data.p = [0, 0, -(interface_point + twr_data.h)]

    # Build floater
    floater_model = model.FloaterClass(floater_data, rho_st)  # CoGz will be set for BOS at z=0

    # Build tower and RNA
    wtg_model = model.WtgClass(twr_data, rna_data, rho_st)  # CoGz will be controlled by p

    # Assemble parts into complete unit
    unit_model = model.UnitClass([wtg_model, floater_model])
    return unit_model, floater_model, wtg_model


def launch_monopole(fio, settings):
    unit_model, floater_model, wtg_model = create_mass_models(settings)

    tb.write_gmsh(fio, floater_model)
    vertices, faces = mmio.load_MSH(fio.gmsh_dir.joinpath('MOLO_{}c.msh'.format(floater_model.nc)))
    start_mesh = Mesh(vertices, faces)
    start_mesh.merge_duplicates()
    start_mesh.heal_normals()
    start_mesh.heal_mesh()
    start_mesh.rotate_z(-np.pi / 2)  # IMPORTANT
    print('Before equilibrium calc')
    unit_model.print_vector_matrix_global()
    hs_floater = hs.Hydrostatics(start_mesh, verbose=True)
    hs_floater.gravity = 9.81
    hs_floater.rho_water = 1025.
    hs_floater.mass = unit_model.mass / 1000  # Give mass in tons
    print('Mass given to hydro is {:5.2f} t'.format(hs_floater.mass))
    hs_floater.gravity_center = -unit_model.inertias.reduction_point
    hs_floater.equilibrate()
    #
    tb.save_M_and_K(fio.data_io_dir, M=unit_model.inertias.mass_matrix_global,
                    MMK=hs_floater.hs_data['stiffness_matrix'])
    # Update model with calculated draft
    unit_model.set_new_reduction_point([0,0,hs_floater.hs_data['draught']])
    #unit_model.inertias.reduction_point = [0, 0, unit_model.inertias.reduction_point[2] + hs_floater.hs_data['draught']]
    print('\n\nEquilibrium calc gives {:5.2f} m draught'.format(hs_floater.hs_data['draught']))
    # hs_floater.show()
    # unit_model.print_vector_matrix_global()

    return unit_model, hs_floater

def launch_dipole(fio, model_input):
    unit_model, floater_model, wtg_model = create_mass_models(model_input)

    tb.write_gmsh(fio, floater_model)
    vertices, faces = mmio.load_MSH(fio.gmsh_dir.joinpath('MOLO_{}c.msh'.format(floater_model.nc)))
    start_mesh = Mesh(vertices, faces)
    start_mesh.merge_duplicates()
    start_mesh.heal_normals()
    start_mesh.heal_mesh()
    start_mesh.rotate_z(-np.pi / 2)  # IMPORTANT
    print('Before equilibrium calc')
    unit_model.print_vector_matrix_global()
    hs_floater = hs.Hydrostatics(start_mesh, verbose=True)
    hs_floater.gravity = 9.81
    hs_floater.rho_water = 1025.
    hs_floater.mass = unit_model.mass / 1000  # Give mass in tons
    print('Mass given to hydro is {:5.2f} t'.format(hs_floater.mass))
    hs_floater.gravity_center = -unit_model.inertias.reduction_point
    hs_floater.equilibrate()
    #
    tb.save_M_and_K(fio.data_io_dir, M=unit_model.inertias.mass_matrix_global,
                    MMK=hs_floater.hs_data['stiffness_matrix'])
    # Update model with calculated draft
    unit_model.set_new_reduction_point([0,0,hs_floater.hs_data['draught']])
    #unit_model.inertias.reduction_point = [0, 0, unit_model.inertias.reduction_point[2] + hs_floater.hs_data['draught']]
    print('\n\nEquilibrium calc gives {:5.2f} m draught'.format(hs_floater.hs_data['draught']))
    # hs_floater.show()
    # unit_model.print_vector_matrix_global()

    return unit_model, hs_floater

