__author__ = "Eivind Sonju"
__copyright__ = "Copyright (C) 2017-2019 Verbun AS. All rights reserved."
__version__ = "2.0"

import json
import os

import numpy as np
import pandas as pd
import pickle
import meshmagick.hydrostatics as hs
import meshmagick.mmio as mmio
import model
import tool_box as tb
from meshmagick.mesh import Mesh
import report

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

    twr_data = model.TowerDataClass(settings.job_data['tower'])
    floater_data = model.FloaterDataClass(settings.job_data['floater'])
    rna_data = model.RNADataClass(settings.job_data['rna'])

    # Set global z
    interface_point = floater_data.t_lf + floater_data.hgt + floater_data.t_uf
    twr_data.p = [0, 0, -(interface_point + twr_data.h / 2)]
    rna_data.p = [0, 0, -(interface_point + twr_data.h)]

    # Build floater
    floater_model = model.FloaterClass(floater_data, settings.rho_st)  # CoGz will be set for BOS at z=0

    # Build tower and RNA
    wtg_model = model.WtgClass(twr_data, rna_data, settings.rho_st)  # CoGz will be controlled by p

    # Assemble parts into complete unit
    unit_model = model.UnitClass([wtg_model, floater_model])
    return unit_model, floater_model, wtg_model


def init_models(settings):  #
    unit_model, floater_model, wtg_model = create_mass_models(settings)

    msh_file = tb.msh_file(settings, mesh_type='stability')
    stability_vertices, stability_panels = mmio.load_MSH(msh_file)

    stability_mesh = Mesh(stability_vertices, stability_panels)

    stability_mesh.merge_duplicates()
    stability_mesh.heal_normals()
    stability_mesh.heal_mesh()
    stability_mesh.rotate_z(-np.pi / 2)  #  # TODO: Re-orient mesh template

    # unit_model.print_vector_matrix_global()
    hs_floater = hs.Hydrostatics(stability_mesh, verbose=True) # TODO: Set mass, gravity and water density here
    hs_floater.gravity = abs(settings.gravity)
    hs_floater.rho_water = settings.rho_sw
    hs_floater.mass = unit_model.mass / 1000  # Give mass in tons
    print('\nMass given to hydro is {:5.2f} t'.format(hs_floater.mass))
    hs_floater.gravity_center = -unit_model.inertias.reduction_point
    # TODO: THIS IS TEMP SPEEDUP. FIX IT !!
    # if settings.do_equilibrate:
    #     hs_floater.equilibrate()
    hs_floater.set_displacement(hs_floater.mass)

    #
    # hs_floater.show()
    print(hs_floater.get_hydrostatic_report())
    settings._report.write_hydrostatic_report_latex_table(hs_floater)
    tb.save_M_and_K(settings.fio.data_io_dir, M=unit_model.inertias.mass_matrix_global,
                    MMK=hs_floater.hs_data['stiffness_matrix'])
    # Update model with calculated draft
    # print('\nMass matrix just before adjusting to draught')
    # for part in unit_model.get_all_parts():
    #     part.print_vector_matrix_global()
    unit_model.move_reduction_point([0, 0, hs_floater.hs_data['draught']])

    for part in unit_model.parts_list:
        part.print_vector_matrix_global()
    settings.draught = hs_floater.hs_data['draught']

    # unit_model.inertias.reduction_point = [0, 0, unit_model.inertias.reduction_point[2] + hs_floater.hs_data['draught']]
    print('\nEquilibrium calc gives {:5.2f} m draught'.format(hs_floater.hs_data['draught']))
    # hs_floater.show()
    # print('\nUpdated global mass matrix after adjusting to draught')
    # for part in unit_model.get_all_parts():
    #     part.print_vector_matrix_global()

    msh_file = tb.msh_file(settings, mesh_type='nemoh')
    nemoh_vertices, nemoh_panels = mmio.load_MSH(msh_file)

    if settings.use_dipols:
        nemoh_vertices, nemoh_panels, not_dipol_index, dipol_index = tb.prepare_dipol_mesh(nemoh_vertices,
                                                                                           nemoh_panels,
                                                                                           settings)
        if len(dipol_index) > 0:
            settings.thin_panels = dipol_index  # index must start with 0
        else:
            settings.thin_panels = '0'

    # print(msh_file.stem)
    nemoh_mesh = Mesh(nemoh_vertices, nemoh_panels)
    nemoh_mesh.heal_normals()
    # nemoh_mesh.show()

    # nemo_mesh_not_dipol= Mesh(nemoh_vertices, nemoh_panels[not_dipol_index])
    # nemo_mesh_dipol= Mesh(nemoh_vertices, nemoh_panels[dipol_index])
    # nemo_mesh_not_dipol.show()
    # nemo_mesh_dipol.show()

    if settings.use_symmmetri:
        nemoh_mesh = nemoh_mesh.de_symmetrize_xz()
        nemoh_mesh.merge_duplicates()

    # nemoh_mesh.show()

    mesh_dat = settings.fio.nemoh_root.joinpath('{}.dat'.format(msh_file.stem))
    settings.mesh_file = str(mesh_dat)

    with open(settings.fio.data_io_dir.joinpath('nemoh_mesh_vertices.pkl'), "wb") as f:
        pickle.dump(nemoh_mesh.vertices, f)

    with open(settings.fio.data_io_dir.joinpath('nemoh_mesh_faces.pkl'), "wb") as f:
        pickle.dump(nemoh_mesh.faces, f)

    mmio.write_MAR(settings.mesh_file, nemoh_mesh.vertices, nemoh_mesh.faces)

    if settings.use_symmmetri:
        with open(mesh_dat, 'r') as f:
            lines = f.readlines()
        # print(lines[:4])
        lines[0] = lines[0].replace('0', '1')
        # print(lines[:4])
        with open(mesh_dat, 'w') as f:
            f.writelines(lines)

    print('\nNemoh mesh written to {}'.format(settings.mesh_file))

    return unit_model, hs_floater
