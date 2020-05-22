__author__ = "Eivind Sonju"
__copyright__ = "Copyright (C) 2017-2019 Verbun AS. All rights reserved."
__version__ = "2.0"

import json
import os

import numpy as np
import pandas as pd
import pickle
import core.meshmagick.hydrostatics as hs
import core.meshmagick.mmio as mmio
from core import model
import core.tool_box as tb
from core.meshmagick.mesh import Mesh
from core import report


def rem_list(df):
    # Dirty :-(
    for item in df:
        try:
            df[item] = df[item][0]
        except:
            pass
    return df


def create_mass_models(settings):
    twr_data = model.TowerDataClass(settings.job_data['wtg'][settings.wtg_model]["Tower"])
    floater_data = model.FloaterDataClass(settings.job_data['floater'])
    rna_data = model.RNADataClass(settings.job_data['wtg'][settings.wtg_model])

    # Set global z
    interface_point = floater_data.t_lf + floater_data.hgt + floater_data.t_uf - settings.draught
    twr_data.p = [0, 0, -(interface_point + twr_data.h / 2)]
    rna_data.p = [0, 0, -(interface_point + twr_data.h)]

    # Build floater

    floater_model = model.FloaterClass(floater_data, settings)  # CoGz will be set for BOS at z=0

    # Build tower and RNA
    wtg_model = model.WtgClass(twr_data, rna_data, settings)  # CoGz will be controlled by p

    # Assemble parts into complete unit
    unit_model = model.UnitClass([wtg_model, floater_model])

    return unit_model, floater_model, wtg_model, floater_data


def init_models(settings):  #
    # Set ballast to zero
    nc=settings.job_data['floater']['Radial']['Number of columns']
    settings.job_data['floater']['Ballast filling ratio'] = 0.0
    unit_model, floater_model, wtg_model, floater_data = create_mass_models(settings)

    for part in unit_model.parts_list:
        part.print_vector_matrix_global()
    unit_model.print_vector_matrix_global()

    msh_file = tb.msh_file(settings, mesh_type='stability')
    stability_vertices, stability_panels = mmio.load_MSH(msh_file)

    stability_mesh = Mesh(stability_vertices, stability_panels)

    stability_mesh.merge_duplicates()
    stability_mesh.heal_normals()
    stability_mesh.heal_mesh()
    stability_mesh.rotate_z(-np.pi / 2)  # Re-orient mesh template

    # unit_model.print_vector_matrix_global()
    hs_floater = hs.Hydrostatics(stability_mesh, verbose=True)  # TODO: Set mass, gravity and water density here
    hs_floater.gravity = abs(settings.gravity)
    hs_floater.rho_water = settings.rho_sw
    hs_floater.mass = unit_model.mass / 1000  # Give mass in tons
    print('\nUn-ballasted mass given to hydro is {:5.2f} t'.format(hs_floater.mass))
    hs_floater.gravity_center = -unit_model.inertias.reduction_point
    # TODO: THIS IS TEMP SPEEDUP. FIX IT !!
    # if settings.do_equilibrate:
    #     hs_floater.equilibrate()
    hs_floater.set_displacement(hs_floater.mass)
    print('Un-ballasted draught is {:1.2f} m'.format(hs_floater.hs_data['draught']))
    print(hs_floater.get_hydrostatic_report())

    d_rc = floater_data.dia_rc
    d_hc = floater_data.dia_hc
    nr = 3
    nc = floater_data.nc
    a_wp = np.pi * (d_hc ** 2 + nr * nc * d_rc ** 2) / 4

    if settings.parameter_space.target_draught:
        target_draught= settings.parameter_space.target_draught
    else:
        target_draught= floater_data.hgt / 2

    m_ball = (target_draught - hs_floater.hs_data['draught']) * a_wp * floater_data.rho_bal

    fr = (m_ball / (a_wp * floater_data.rho_bal)) / floater_data.hgt
    print(
        '\nRequired ballast to {:1.2f} m is {:1.1f} ton\nFilling ratio is {:1.2f}'.format(target_draught, m_ball / 1000,
                                                                                         fr))
    settings.job_data['floater']['Ballast filling ratio'] = fr

    print("\nRecreate mass model with target ballast")
    unit_model, floater_model, wtg_model, floater_data = create_mass_models(settings)
    hs_floater = hs.Hydrostatics(stability_mesh, verbose=True)  # TODO: Set mass, gravity and water density here
    hs_floater.gravity = abs(settings.gravity)
    hs_floater.rho_water = settings.rho_sw
    hs_floater.mass = unit_model.mass / 1000  # Give mass in tons
    print('\nBallasted mass given to hydro is {:5.2f} t'.format(hs_floater.mass))
    hs_floater.gravity_center = -unit_model.inertias.reduction_point
    #hs_floater.equilibrate()
    hs_floater.set_displacement(hs_floater.mass)

    settings.job_data['floater']['Ballast filling ratio'] = fr # TODO: Why have to be set twice, see above. Check setter/getter
    settings.draught = hs_floater.hs_data['draught']
    print('\nRecreate mass model for ballasted draught = {:5.2f}m'.format(settings.draught))
    unit_model, floater_model, wtg_model, floater_data = create_mass_models(settings)

    #hs_floater.show()



    tb.save_M_and_K(settings, M=unit_model.inertias.mass_matrix_global,
                    MMK=hs_floater.hs_data['stiffness_matrix'])

    for part in unit_model.parts_list:
        part.print_vector_matrix_global()
    unit_model.print_vector_matrix_global()
    print(hs_floater.get_hydrostatic_report())
    settings._report.write_hydrostatic_report_latex_table(hs_floater)

    # unit_model.inertias.reduction_point = [0, 0, unit_model.inertias.reduction_point[2] + hs_floater.hs_data['draught']]
    # print('\nEquilibrium calc gives {:5.2f} m draught'.format(hs_floater.hs_data['draught']))
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
