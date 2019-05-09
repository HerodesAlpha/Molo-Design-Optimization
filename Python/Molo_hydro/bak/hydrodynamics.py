import numpy as np
from numpy import tan, arctan2, sqrt, pi, array, cross


def mmi_hollow_cylinder_xy(m, do, t, h):
    return (m / 12) * (3 * ((do / 2) ** 2 + (do / 2 - t) ** 2) + h ** 2)


def hydro_coefficients():
    pass


def hydro_loads():
    pass


def hydro_model(molo_model, wtg_model):
    """Calculate the added mass terms, restoring terms and wave excitation loads for MOLO Model 1, Three radial,
    Single layer

    Units are m, N, s, kg, Pa, J

    Model coordinate system is cartesian with z=0 at waterline and x,y at center of tower.
    Vertical axiz is z, positive upwards.
    x points outwards along radial. y follows right-hand rule.

    Args:
        dia_rc   : Diameter of buoyancy column
        thi_b   : Thickness of buoyancy column
        dia_hc   : Diameter of centre column
        thi_c   : Thickness of centre column
        hgt     : Distance from floater BOS to TOS
        gap     : Gap between columns given as ratio on radius. Default = 0.8
        ncolumns     : Total number of buoyancy elements

    Returns:
        M       : Mass matrix
        A       : Added mass matrix
        C       : Restoring matrix
        F       : Excitation force



    """
    # Defaults
    nr = 3  # Three equidistant radials

    # Init
    nc = molo_model['radial']['ncolumn']
    gap = molo_model['radial']['gap']
    hgt = molo_model['hgt']
    dia_rc = molo_model['radial']['column']['d']
    thi_rc = molo_model['radial']['column']['t']
    dia_hc = molo_model['hub']['column']['d']
    thi_hc = molo_model['hub']['column']['t']
    w_lf = molo_model['radial']['lower_flange']['w']
    t_lf = molo_model['radial']['lower_flange']['t']
    w_uf = molo_model['radial']['upper_flange']['w']
    t_uf = molo_model['radial']['upper_flange']['t']
    m_wtg = wtg_model['mass']
    m44_wtg = wtg_model['inertia'][3][3]
    cog_z_wtg = wtg_model['cog']['z']
    ballast = molo_model['ballast']

    M = []
    A = []
    C = []

    # -------------------------------------------------------------
    # PHYSICAL CONSTANTS
    # -------------------------------------------------------------
    #
    rho_air = 1.3
    rho_steel = 7850
    rho_sw = 1025
    rho_ballast = rho_sw
    grav = 9.806
    #
    # -------------------------------------------------------------
    # DERIVED GEOMETRICAL PROPERTIES
    # -------------------------------------------------------------
    #
    l_radial = nc * (1 + gap) * dia_rc
    #
    # -------------------------------------------------
    # MASS
    # -------------------------------------------------
    #
    # Calculate local mass
    m_rc = pi * dia_rc * thi_rc * hgt * rho_steel
    m_hc = pi * dia_hc * thi_hc * hgt * rho_steel
    m_lf = l_radial * w_lf * t_lf * rho_steel
    m_uf = l_radial * w_uf * t_uf * rho_steel
    m_radial = nc * m_rc + m_lf + m_uf
    m_hub = m_hc
    m_ballast = nr * nc * ballast * pi * (dia_rc - 2 * thi_rc) ** 2 / 4 * rho_ballast
    m_molo = m_hub + nr * m_radial + m_ballast
    #

    # Assemble total mass
    m_tot = m_molo + m_wtg
    #
    # -------------------------------------------------
    # COG
    # -------------------------------------------------
    #
    # Calculate vertical centre of gravity relative to BOS
    cogz = (m_wtg * cog_z_wtg + nr * (
            m_lf * t_lf / 2 + nc * m_rc * (hgt / 2) + m_uf * (hgt - t_lf / 2))) / m_tot
    #
    # -------------------------------------------------
    # WATERPLANE AREA
    # -------------------------------------------------
    #
    # Total waterplane area
    a_wp = nc * nr * pi * dia_rc ** 2 / 4 + pi * dia_hc ** 2 / 4
    c33 = a_wp * rho_sw * grav
    # Second area moment - sam
    sam_xy = pi * dia_hc ** 2 / 4
    sam_xy += nc * nr * pi * dia_rc ** 2 / 4
    r0 = dia_hc / 2 + (0.5 + gap) * dia_rc
    sam_xy += pi * dia_rc ** 2 / 4 * 3 * r0 ** 2 / 2
    for ic in range(nc - 1):
        r0 += (1 + gap) * dia_rc
        sam_xy += pi * dia_rc ** 2 / 4 * 3 * r0 ** 2 / 2
    #
    # -------------------------------------------------
    # DRAFT
    # -------------------------------------------------
    #
    z_bos = -(m_tot / (rho_sw * a_wp))
    #
    # -------------------------------------------------
    # INERTIA
    # -------------------------------------------------
    # Mass moment of inertia - mmi
    #
    # Flanges
    mmi_xy_lf = m_lf * l_radial / 2
    mmi_xy_uf = m_uf * l_radial / 2
    #
    # Assemble xy mass moment of inertia for columns
    # 1) Add center column
    mmi_xy_columns = mmi_hollow_cylinder_xy(m_hc, dia_hc, thi_hc, hgt)
    # 2) Add radial columns
    mmi_xy_columns += nr * nc * mmi_hollow_cylinder_xy(m_rc, dia_rc, thi_rc, hgt)
    # 3) Add off axis CoG contribution for each radial column
    r0 = dia_hc / 2 + (0.5 + gap) * dia_rc
    mmi_xy_columns += m_rc * 3 * r0 ** 2 / 2
    for ic in range(nc - 1):
        r0 += (1 + gap) * dia_rc
        mmi_xy_columns += m_rc * 3 * r0 ** 2 / 2
    #
    # Move axis to waterline and sum up for all parts
    mmi_xy_lf_wl = mmi_xy_lf + nr * m_lf * (-z_bos - t_lf / 2) ** 2
    mmi_xy_uf_wl = mmi_xy_uf + nr * m_uf * (-z_bos - hgt + t_uf / 2) ** 2
    mmi_xy_rc_wl = mmi_xy_columns + nr * nc * m_rc * (-z_bos - hgt / 2) ** 2
    mmi_xy_molo_wl = mmi_xy_lf_wl + mmi_xy_uf_wl + mmi_xy_rc_wl
    mmi_xy_tot_wl = mmi_xy_molo_wl + m44_wtg + m_wtg * (hgt - z_bos) ** 2
    #
    # Inertia tensor; small inbalance roll-yaw (i46) due to RNA eccentricity
    m44 = mmi_xy_tot_wl
    m55 = mmi_xy_tot_wl
    m66 = 0  # NOT IMPLEMENTED YET
    m46 = 0  # NOT IMPLEMENTED YET
    #
    # ------------------------------------------------
    # MASS PRODUCT OF INERTIA
    # -------------------------------------------------
    #

    #
    # ------------------------------------------------
    # ADDED MASS
    # -------------------------------------------------
    #
    # Two dimensional added mass for lower flange
    a33_2d_lf = rho_sw * pi * w_lf ** 2 / 4
    #
    # Added mass in heave from flange
    a33_lf = a33_2d_lf * l_radial * nr
    #
    a33 = a33_lf
    # Added mass in roll/pitch from flange
    a55_lf = a33_2d_lf * l_radial ** 2 / 2
    #
    a55 = a55_lf
    #
    # ------------------------------------------------
    # METACENTRIC HEIGHT
    # -------------------------------------------------
    #
    # Volume of displacement - vd
    vd = m_tot / rho_sw
    kb = hgt / 2
    bm = sam_xy / vd
    kg = cogz - z_bos
    km = kb + bm
    gm = km - kg
    #
    # ------------------------------------------------
    # RADIUS OF INERTIA
    # -------------------------------------------------
    rii_xy_wl = sqrt(mmi_xy_molo_wl / (m_molo))
    rii_xy_wl_tot = sqrt(mmi_xy_tot_wl / (m_tot))
    rii_wp = sqrt(sam_xy / a_wp)
    #
    # ------------------------------------------------
    # RESTORING COEFFICIENTS
    # -------------------------------------------------
    #
    c33 = rho_sw * grav * a_wp
    c55 = rho_sw * grav * vd * gm
    #
    # ------------------------------------------------
    # NATURAL PERIOD
    # -------------------------------------------------
    #
    tn3 = 2 * pi * sqrt((m_tot + a33) / c33)
    tn5 = 2 * pi * sqrt((m55 + a55) / c55)
    #
    #
    # ------------------------------------------------
    # MAX STATIC TRUST FROM WTG
    # -------------------------------------------------
    #
    d_rot = wtg_model['Rotor diameter']
    # Rotor disk area - a_rot
    a_rot = pi * d_rot ** 2 / 4
    ct = wtg_model['Trust coefficient']
    u_ref = wtg_model['Max trust windspeed']
    # Trust force - ft
    ft = 0.5 * rho_air * ct * a_rot * u_ref ** 2
    # Rated power - rp
    rp = a_rot * wtg_model['Rated power per unit area']
    # Overturning moment at waterline - otm
    wtg_hub_z = wtg_model['Hub height above tower bottom']
    otm = ft * (wtg_hub_z + hgt + z_bos)
    theta = otm / c55
    #
    # ------------------------------------------------
    # COLLECT PROPERTIES
    # -------------------------------------------------
    wtg_properties = dict()
    wtg_properties['d_rot'] = {'text': 'Rotor diameter', 'val': d_rot, 'unit': 'm'}
    wtg_properties['m_wtg'] = {'text': 'WTG Mass', 'val': m_wtg / 1000, 'unit': 'ton'}
    wtg_properties['cog_z_wtg'] = {'text': 'CoGz WTG above bottom of tower', 'val': cog_z_wtg, 'unit': 'm'}
    wtg_properties['wtg_hub_z'] = {'text': 'Hub height above tower bottom', 'val': wtg_hub_z, 'unit': 'm'}
    #
    molo_properties = dict()
    molo_properties['l_radial'] = {'text': 'Length of radial from face of center column', 'val': l_radial, 'unit': 'm'}
    molo_properties['m_ballast'] = {'text': 'Mass of ballast', 'val': m_ballast / 1000, 'unit': 'ton'}
    molo_properties['m_molo'] = {'text': 'Mass of MOLO', 'val': m_molo / 1000, 'unit': 'ton'}
    molo_properties['m_tot'] = {'text': 'Total mass', 'val': m_tot / 1000, 'unit': 'ton'}
    molo_properties['a_wp'] = {'text': 'Total waterplane area', 'val': a_wp, 'unit': 'm^2'}
    molo_properties['z_bos'] = {'text': 'Draft', 'val': z_bos, 'unit': 'm'}
    molo_properties['a33_lf'] = {'text': 'A33 for lower flange', 'val': a33_lf / 1000, 'unit': 'ton'}
    molo_properties['a55_lf'] = {'text': 'A55 for lower flange', 'val': a55_lf / 1000, 'unit': 'ton-m^2'}
    molo_properties['cogz'] = {'text': 'CoGz', 'val': cogz, 'unit': 'm'}
    molo_properties['rii_xy_wl'] = {'text': 'Horizontal axis radius of inertia for MOLO', 'val': rii_xy_wl, 'unit': 'm'}
    molo_properties['rii_xy_wl_tot'] = {'text': 'Total radius of inertia', 'val': rii_xy_wl_tot, 'unit': 'm'}
    molo_properties['gm'] = {'text': 'Metacentric height', 'val': gm, 'unit': 'm'}
    molo_properties['tn3'] = {'text': 'Natural period in heave', 'val': tn3, 'unit': 's'}
    molo_properties['tn5'] = {'text': 'Natural period in pitch', 'val': tn5, 'unit': 's'}
    molo_properties['otm'] = {'text': 'Overturning moment', 'val': otm / 1000 / 9.81, 'unit': 'ton-m'}
    molo_properties['theta'] = {'text': 'Static pitch', 'val': theta * 180 / pi, 'unit': 'deg'}
    molo_properties['sam_xy'] = {'text': 'Second area moment waterplane', 'val': sam_xy, 'unit': 'm^4'}
    molo_properties['rii_wp'] = {'text': 'Radius of waterplane area', 'val': rii_wp, 'unit': 'm'}
    #
    hydro_output = dict()
    hydro_output['wtg'] = {'name': 'Rotor-Nacelle Assembly and Tower', 'properties': wtg_properties}
    hydro_output['molo'] = {'name': 'Molo Floater', 'properties': molo_properties}
    #
    return hydro_output


def molo_model(d_rc=7, d_hc=7):
    out = {'hgt': 14,
           'radial': {
               'ncolumn': 4,
               'gap': 0.8,
               'column': {'d': d_rc, 't': 0.005 * d_rc},
               'lower_flange': {'w': d_rc, 't': 0.05},
               'upper_flange': {'w': d_rc, 't': 0.05}},
           'hub': {
               'column': {'d': d_hc, 't': 0.005 * d_hc},
               'lower_flange': {'w': 7, 't': 0.05},
               'upper_flange': {'w': 7, 't': 0.05}},
           'ballast': 3
           }
    return out


def wtg_model():
    d_rot = 167
    z_hub = (d_rot / 2 + 35)
    d_twr = 6
    t_twr = 0.05
    m_twr = pi * d_twr * t_twr * 7850 * z_hub
    m_rna = 500000
    m_wtg = m_twr + m_rna
    cog_z_wtg = (m_rna * z_hub + m_twr * z_hub / 2) / (m_wtg)
    inertia_xy_twr = mmi_hollow_cylinder_xy(m_twr, d_twr, t_twr, z_hub)
    inertia_xy_wtg_bos = inertia_xy_twr + m_twr * (z_hub / 2) ** 2 + m_rna * z_hub ** 2

    # print('WTG CoGz relative to bottom of tower is {:1.2f} m'.format(cog_z_wtg))
    # print('WTG radius of inertia for bottom of tower is {:1.1f} m'.format(
    #    sqrt(inertia_xy_wtg_bos / (m_wtg))))

    wtg_inertia = np.diagflat([m_wtg, m_wtg, m_wtg, inertia_xy_wtg_bos, inertia_xy_wtg_bos, m_wtg * 2 ** 2])

    return {'mass': m_twr + m_rna,
            'cog': {'x': 0, 'y': 0, 'z': cog_z_wtg},
            'inertia': wtg_inertia,
            'Rotor diameter': d_rot,
            'Trust coefficient': 0.7,
            'Max trust windspeed': 10,
            'Rated power per unit area': 420,
            'Hub height above tower bottom': z_hub
            }

# M, A, C = hydro_model(molo_model(), wtg_model())
