__author__ = "Eivind Sonju"
__copyright__ = "Copyright (C) 2017-2019 Verbun AS. All rights reserved."
__version__ = "2.0"

import os
import subprocess
from pathlib import Path

import numpy as np

import meshmagick.mesh_clipper as mc
import meshmagick.mmio as mmio
import tool_box as tb
from meshmagick.mesh import Mesh


def writeCalFile(dir, rho_sw, water_depth, omega, dof, aO, mesh_file, nrNode, nrPanel):
    # In case of array simulation, do stuff
    dirCheck = aO['dirCheck']
    irfCheck = aO['irfCheck']
    kochCheck = aO['kochCheck']
    fsCheck = aO['fsCheck']

    # Create the Nemoh calibration file
    with open('Nemoh.cal', 'w') as fid:
        fid.write('--- Environment ---\n')
        fid.write(str(rho_sw) + '				! RHO 			! KG/M**3 	! Fluid specific volume\n')
        fid.write('9.81				! G			! M/S**2	! Gravity\n')
        fid.write(str(water_depth) + '				! DEPTH			! M		! Water depth\n')
        fid.write('0.	0.			! XEFF YEFF		! M		! Wave measurement point\n')
        fid.write('--- Description of floating bodies ---\n')
        fid.write('{:d}				! Number of bodies\n'.format(1))
        fid.write('--- Body 1 ---------------------------\n')
        fid.write('{}\t			! Name of mesh file\n'.format(mesh_file))
        fid.write(str(nrNode) + '\t' + str(nrPanel) + '			! Number of points and number of panels\n')
        fid.write('{:d}				! Number of degrees of freedom\n'.format(sum(dof)))
        for iDof in range(len(dof)):
            if (iDof == 0 and dof[iDof] == 1):
                fid.write('1 1. 0.	0. 0. 0. 0.		! Surge\n')
            elif (iDof == 1 and dof[iDof] == 1):
                fid.write('1 0. 1.	0. 0. 0. 0.		! Sway\n')
            elif (iDof == 2 and dof[iDof] == 1):
                fid.write('1 0. 0. 1. 0. 0. 0.		! Heave\n')
            elif (iDof == 3 and dof[iDof] == 1):
                fid.write('2 1. 0. 0. 0. 0. 0.		! Roll about waterline\n')
            elif (iDof == 4 and dof[iDof] == 1):
                fid.write('2 0. 1. 0. 0. 0. 0.		! Pitch about waterline\n')
            elif (iDof == 5 and dof[iDof] == 1):
                fid.write('2 0. 0. 1. 0. 0. 0.		! Yaw about waterline\n')
        fid.write('{:d}				! Number of resulting generalised forces\n'.format(sum(dof)))
        for iDof in range(len(dof)):
            if (iDof == 0 and dof[iDof] == 1):
                fid.write('1 1. 0.	0. 0. 0. 0.		! Force in X direction\n')
            elif (iDof == 1 and dof[iDof] == 1):
                fid.write('1 0. 1.	0. 0. 0. 0.		! Force in Y direction\n')
            elif (iDof == 2 and dof[iDof] == 1):
                fid.write('1 0. 0. 1. 0. 0. 0.		! Force in Z direction\n')
            elif (iDof == 3 and dof[iDof] == 1):
                fid.write('2 1. 0. 0. 0. 0. 0.		! Roll moment about waterline\n')
            elif (iDof == 4 and dof[iDof] == 1):
                fid.write('2 0. 1. 0. 0. 0. 0.		! Pitch moment about waterline\n')
            elif (iDof == 5 and dof[iDof] == 1):
                fid.write('2 0. 0. 1. 0. 0. 0.		! Yaw moment about waterline\n')
        fid.write('0				! Number of lines of additional information\n')

        fid.write('--- Load cases to be solved ---\n')
        fid.write(str(omega[0]) + '\t' + str(omega[1]) + '\t' + str(
                omega[2]) + '		! Number of wave frequencies, Min, and Max (rad/s)\n')
        if dirCheck:
            fid.write(str(aO['dirStep']) + '\t' + str(aO['dirStart']) + '\t' + str(
                    aO['dirStop']) + '		! Number of wave directions, Min and Max (degrees)\n')
        else:
            fid.write(str(dir[0]) + '\t' + str(dir[1]) + '\t' + str(
                    dir[2]) + '		! Number of wave directions, Min and Max (degrees)\n')
        fid.write('--- Post processing ---\n')
        if irfCheck:
            fid.write('1' + '\t' + str(aO['irfStep']) + '\t' + str(aO[
                                                                       'irfDur']) + '\t\t! IRF 				! IRF calculation (0 for no calculation), time step and duration\n')
        else:
            fid.write(
                    '0' + '\t0.01\t20.\t\t! IRF 				! IRF calculation (0 for no calculation), time step and duration\n')
        fid.write('1				! Show pressure\n')
        if kochCheck:
            fid.write(str(aO['kochStep']) + '\t' + str(aO['kochStart']) + '\t' + str(aO[
                                                                                         'kochStop']) + '		! Kochin function 		! Number of directions of calculation (0 for no calculations), Min and Max (degrees)\n')
        else:
            fid.write(
                    '0	0.	180.		! Kochin function 		! Number of directions of calculation (0 for no calculations), Min and Max (degrees)\n')
        if fsCheck:
            fid.write(str(aO['fsDeltaX']) + '\t' + str(aO['fsDeltaY']) + '\t' + str(aO['fsLengthX']) + '\t' + str(aO[
                                                                                                                      'fsLengthY']) + '	! Free surface elevation 	! Number of points in x direction (0 for no calcutions) and y direction and dimensions of domain in x and y direction	\n')
        else:
            fid.write(
                    '50	50	100 100	! Free surface elevation 	! Number of points in x direction (0 for no calcutions) and y direction and dimensions of domain in x and y direction	\n')


def runNemoh(fio, hydromodel, mesh_file, dir, rho_sw, water_depth, omega_calc, dof):
    curDir = os.getcwd()
    os.chdir(fio.nemoh_root)
    with open('ID.dat', 'w') as f:
        f.write('1\n.\n')
    with open('input.txt', 'w') as f:
        f.write('Calculation parameters\n0\n')

    aO = {'dirCheck': False, 'irfCheck': False, 'kochCheck': False, 'fsCheck': False}
    nrNode = hydromodel.nb_vertices
    nrPanel = hydromodel.nb_faces

    writeCalFile(dir, rho_sw, water_depth, omega_calc, dof, aO, mesh_file, nrNode, nrPanel)

    nemoh_exe_path = Path(r'C:\Users\eison\OneDrive - Verbun AS\Divisions\Offshore Wind\Library\Software\Nemoh v2.03')
    nemoh_exe_path = Path(
            r'C:\Users\eison\OneDrive - Verbun AS\Divisions\Offshore Wind\Projects\P2017.001\Work\Fortran\nemoh\bin')
    subprocess.call(str(nemoh_exe_path.joinpath('preProc.exe')))
    subprocess.call(str(nemoh_exe_path.joinpath('solver.exe')))
    subprocess.call(str(nemoh_exe_path.joinpath('postProc.exe')))
    os.chdir(curDir)


def calcAlphaBeta(wdir, nSel, rho, dof):
    # ---------------------------------------------------------------------------
    # INPUT
    # ---------------------------------------------------------------------------

    # Open IRF file
    # rho = 1025.0
    # nSel = 10
    # dof = 'all'
    pathFile = os.path.join(wdir, 'Nemoh', 'IRF.tec')
    with open(pathFile, 'r') as f:
        irfRaw = f.readlines()
    if sum(dof) < 2:
        strPat = 'DoF    1'
        irfInd = 2
    elif dof[2] == 1:
        strPat = 'DoF    {:d}'.format(sum(dof[0::2]))
        irfInd = sum(dof[0::2]) * 2
    else:
        strPat = 'DoF    1'
        irfInd = 2
    for iL in range(0, len(irfRaw)):
        irfInfo = irfRaw[iL]
        test = irfInfo.find(strPat)
        if test != -1:
            indStart = iL
            irfInfo = irfInfo.split()
            irfNrlong = irfInfo[3]
            indComma = irfNrlong.index(',')
            irfNr = int(irfNrlong[0:indComma])
    # Arrange data in IRF file
    time = [0] * irfNr
    IRF = [0] * irfNr
    for iL in range(indStart + 1, irfNr + indStart + 1):
        irfLine = irfRaw[iL].split()
        time[iL - indStart - 1] = float(irfLine[0])
        IRF[iL - indStart - 1] = float(irfLine[irfInd])
    aInf = float(irfLine[irfInd - 1])

    # ---------------------------------------------------------------------------
    # CALCULATION
    # ---------------------------------------------------------------------------

    IRF = np.array(IRF)
    dt = time[2] - time[1]
    m = 200
    n = len(IRF) - m - 1

    # STEP1: make polynome R with values s_k
    # -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~

    A = np.zeros((n + 1, m))
    for i in range(0, n + 1):
        A[i, :] = IRF[i:i + m]
    C = -IRF[m:n + m + 2]
    R, res, ran, sing = np.linalg.lstsq(A, C)
    del res, ran, sing
    R = np.append(R, 1)
    R = np.flipud(R)

    # STEP2: find roots of polynomial R
    # -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~

    Q = np.roots(R)

    # STEP3: Calculate Beta
    # -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~

    beta = np.log(Q) / dt
    del A, C, R

    # STEP3: Calculate Alpha
    # -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~

    Q2 = np.zeros((n + 1, m));
    Q2 = np.array(Q2, dtype=complex)
    for i in range(0, n + 1):
        Q2[i, :] = np.transpose(Q) ** i
    F = IRF[0:n + 1]
    C, res, ran, sing = np.linalg.lstsq(Q2, F)
    del res, ran, sing
    alpha = C * np.exp(-beta * time[1])
    del F, Q, Q2, C

    # STEP5: only consider beta-values with a negative real part, so function goes
    # to zero
    # -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~

    mask = np.real(beta) < 0
    alpha = alpha[mask]
    beta = beta[mask]
    m = len(alpha)

    # ---------------------------------------------------------------------------
    # SELECTION of alpha en beta components
    # ---------------------------------------------------------------------------

    # STEP1: Sorting of values
    # -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~
    absVal = -np.abs(np.divide(alpha, beta))
    indSort = np.flipud(np.argsort(absVal))
    alpha = np.flipud(alpha[indSort])
    beta = np.flipud(beta[indSort])

    # STEP2: Complex conjugate checking
    # -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~

    test1 = np.abs(np.imag(alpha[nSel - 1]))
    test2 = np.abs(np.imag(alpha[nSel]))
    if np.abs(test1 - test2) < 0.001:
        nSel = nSel + 1
    alphaSel = alpha[0:nSel]
    betaSel = beta[0:nSel]

    # ---------------------------------------------------------------------------
    # CALCULATE Error
    # ---------------------------------------------------------------------------

    # STEP1: Reconstruction of the IRF
    # -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~

    mat1 = np.matrix(np.ones((len(IRF), 1)))
    mat2 = np.matrix(alphaSel)
    mat3 = np.transpose(np.matrix(time)) * np.matrix(betaSel)
    fMatrix = np.multiply(mat1 * mat2, np.exp(mat3))
    fRecons = np.sum(fMatrix, axis=1)

    # STEP2: Calculation of the error
    # -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~

    diff1 = fRecons[1:len(IRF) / 2]
    diff2 = np.transpose(np.matrix(IRF[1:len(IRF) / 2]))
    diff3 = np.matrix(IRF)
    relDiffE = np.abs(diff1 - diff2) / np.max(np.abs(diff3))
    error2E = np.abs(np.sum(relDiffE) / ((len(IRF)) / 2))

    # ---------------------------------------------------------------------------
    # PLOT Result
    # ---------------------------------------------------------------------------

    # plt.plot(time,IRF)
    # plt.plot(time,np.real(fRecons))

    return (alphaSel, betaSel, error2E, aInf)


def calcM(wdir, rho=1025.0, dof=[0, 0, 1, 0, 0, 0]):
    pathFile = os.path.join(wdir, 'results', 'KH.dat')
    dof = np.array(dof)
    with open(pathFile, 'r') as f:
        khRaw = f.readlines()
    KH = np.zeros((6, 6))
    for iL in range(len(khRaw)):
        KH[iL, :] = np.array([float(a) for a in khRaw[iL].split()])
    if sum(dof) > 1:
        c = KH
    else:
        dof = list(dof)
        c = np.dot(KH, dof)[dof.index(1)]

    M = np.diag([1., 1., 1., 1., 1., 1.])
    pathFile = os.path.join(wdir, 'Nemoh', 'Hydrostatics.dat')
    with open(pathFile, 'r') as f:
        mRaw = f.readlines()
    mass = rho * float(mRaw[3].split()[2])
    M[0:3, 0:3] *= mass
    zG = float(mRaw[2].split()[-1])
    M[0, 4] = mass * zG
    M[1, 3] = -mass * zG
    M[4, 0] = mass * zG
    M[3, 1] = -mass * zG

    I = np.zeros((3, 3))
    pathFile = os.path.join(wdir, 'Nemoh', 'Inertia_hull.dat')
    with open(pathFile, 'r') as f:
        mRaw = f.readlines()
    for iL in range(len(mRaw)):
        I[iL, :] = np.array([float(a) for a in mRaw[iL].split()])
    M[3::, 3::] *= I
    if sum(dof) > 1:
        None
    else:
        M = np.dot(M, dof)[dof.index(1)]
    return (M, c)


def writeOutputNemoh(wdir, M, Mainf, c, alpha, beta):
    pathFile = os.path.join(wdir, 'Nemoh', 'outNemoh.dat')
    fid = open(pathFile, 'w')
    inString = '{0:e}'
    fid.write('Mass:\t' + inString.format(M) + '\n')
    fid.write('Mainf:\t' + inString.format(Mainf) + '\n')
    fid.write('c:\t' + inString.format(c) + '\n')
    fid.write('=========================================================\n')
    fid.write('ALPHA\n')
    fid.write('=========================================================\n')
    for iA in range(0, len(alpha)):
        fid.write(str(np.real(alpha[iA])) + '\t' + str(np.imag(alpha[iA])) + '\n')
    fid.write('=========================================================\n')
    fid.write('BETA\n')
    fid.write('=========================================================\n')
    for iB in range(0, len(beta)):
        fid.write(str(np.real(beta[iB])) + '\t' + str(np.imag(beta[iB])) + '\n')
    fid.close()


def getOmega(nemoh_results_directory):
    pathFile = nemoh_results_directory.joinpath('Fe.dat')
    omega = list()
    with open(pathFile, 'r') as f:
        irfRaw = f.readlines()
        # Find all occurences of separation lines containing 'Zone t="Diffraction force - beta ='
        first_sep_found = False
        for line in irfRaw:
            if r'Zone t="Diffraction force - beta =' in line:
                if not first_sep_found:
                    first_sep_found = True
                else:
                    break
            elif first_sep_found:
                omega.append(float(line.split()[0]))
    omega = np.asfarray(omega)
    return omega


def getDirections(nemoh_results_directory):
    pathFile = nemoh_results_directory.joinpath('Fe.dat')
    Dir = list()
    with open(pathFile, 'r') as f:
        irfRaw = f.readlines()
        for line in irfRaw:
            if r'Zone t="Diffraction force - beta =' in line:
                Dir.append(float(line[35:].split()[0]))
    return Dir


def get_ab(fio, dof, omega, dir, nbody=1):
    ndir = len(dir)
    nomega = len(omega)
    ndof = sum(dof)
    # Open Radiation Coefficients
    pathFile = fio.nemoh_results.joinpath('RadiationCoefficients.tec')
    with open(pathFile, 'r') as f:
        irfRaw = f.readlines()
    # Remove header text
    irfRaw = irfRaw[1 + sum(dof)::]
    # Remove table labels, one per direction and dof
    del irfRaw[::nomega + 1]
    # Map to numpy array
    x = np.array([line.split() for line in irfRaw], dtype='float')
    # Remove freq. vector
    x = x[:, 1:]
    # Split in A and B
    a = x[:, ::2]
    b = x[:, 1:][:, ::2]
    # TODO: Reshape on ndir>1 not validated
    # TODO: Return 6x6
    a = a.reshape(ndof, nomega, ndof).swapaxes(0, 1)

    # print(a)
    b = b.reshape(ndof, nomega, ndof).swapaxes(0, 1)

    res_a = np.zeros([nomega, 6, 6], dtype=float)
    res_b = np.zeros([nomega, 6, 6], dtype=float)
    idof = np.array([i for i, x in enumerate(dof) if x])
    # print(idof)
    for freq in range(nomega):
        for i in range(len(idof)):
            for j in range(len(idof)):
                res_a[freq, idof[i], idof[j]] = a[freq, i, j]
                res_b[freq, idof[i], idof[j]] = b[freq, i, j]
    return res_a, res_b


def get_fe(fio, dof, omega, dir, nbody=1):
    ndir = len(dir)
    nomega = len(omega)
    ndof = sum(dof)
    # Open Radiation Coefficients
    pathFile = fio.nemoh_results.joinpath('ExcitationForce.tec')
    with open(pathFile, 'r') as f:
        irfRaw = f.readlines()
    irfRaw = irfRaw[1 + sum(dof)::]
    # Remove table labels, one per direction
    del irfRaw[::nomega + 1]
    # Map to numpy array
    x = np.asfarray([line.split() for line in irfRaw])
    # Remove freq. vector
    x = x[:, 1:]
    # Split in mod and angle
    FeMod = x[:, ::2]
    FeAng = x[:, 1:][:, ::2]
    # Reshape on directions, 2D to 3D
    FeMod = FeMod.reshape(ndir, nomega, ndof)
    FeAng = FeAng.reshape(ndir, nomega, ndof)
    res = np.zeros([ndir, nomega, 6], dtype=complex)
    idof = np.array([i for i, x in enumerate(dof) if x])
    Fe_cmplx = tb.P2R(FeMod, FeAng)
    res[:, :, idof] = Fe_cmplx[:, :, :]
    return res


def get_fk_pressure(fio, ndir, nomega, npanels):
    pathFile = fio.nemoh_results.joinpath('FKPressure.dat')
    with open(pathFile, 'r') as f:
        lines = f.readlines()
        if len(lines) - 1 == ndir * nomega * npanels:
            words = [line.split() for line in lines[1:]]
            # index = np.asarray(words[0:3, :], dtype='int')
            # # coordinates = np.asarray(words[3:6,:], dtype='float')
            # polar_pressure = np.asarray(words[6:, :], dtype='float')
            pressure_polar = np.loadtxt(pathFile, skiprows=1,
                                        usecols=[6, 7])
            pressure_cmplx = tb.P2R(pressure_polar[:, 0], pressure_polar[:, 1])
            # cmplx_pressure = tb.P2R(polar_pressure[:, 0], polar_pressure[:, 1])

            return pressure_cmplx.reshape(nomega, ndir, npanels)

        else:
            print('Number of lines on FKPressure.dat mismatch')


def get_poten_pressure(fio, npoints, ppanels, problem_number):
    npanels = ppanels.shape[0]
    pressure_cmplx_panels = np.zeros(npanels, dtype=np.complex)
    pathFile = fio.nemoh_results.joinpath('pressure.{:5d}.dat'.format(problem_number))
    with open(pathFile, 'r') as f:
        lines = f.readlines()
    # Number of vertices and number of panels
    ls = lines[1].split(',')
    if int(ls[1][3::]) == npanels:
        pressure_polar_vertices = np.asarray([line.split() for line in lines[2:npoints + 2]], dtype='float')
        pressure_cmplx_vertices = tb.P2R(pressure_polar_vertices[:, 3],
                                         pressure_polar_vertices[:, 4])  # pressure_panels - 1
        for i, panel in enumerate(ppanels):
            pressure_cmplx_panels[i] = sum([pressure_cmplx_vertices[j] for j in panel]) / 4
        return pressure_cmplx_panels
    else:
        print('Number of panels on pressurefile mismatch')


def mesh(fio, hs_floater, sym=None):
    MeshClipper_instance = mc.MeshClipper(hs_floater.mesh)

    hydro_mesh = Mesh(MeshClipper_instance.lower_mesh.vertices, MeshClipper_instance.lower_mesh.faces)
    hydro_mesh.merge_duplicates()
    hydro_mesh.heal_normals()
    hydro_mesh.heal_mesh()

    if sym == None or sym == 0:
        nemoh_mesh = hydro_mesh
    elif sym == 1:
        nemoh_mesh = hydro_mesh.de_symmetrize_xz()
        nemoh_mesh.merge_duplicates()
        nemoh_mesh.heal_normals()
        nemoh_mesh.heal_mesh()
    else:
        print('sym cannot be {}'.format(sym))
        exit()

    print('Number of faces: {}'.format(nemoh_mesh.nb_faces))
    print('Number of vertices: {}'.format(nemoh_mesh.nb_vertices))

    mesh_file = 'mesh.dat'
    mmio.write_MAR(os.path.join(fio.nemoh_root, mesh_file), nemoh_mesh.vertices, nemoh_mesh.faces)

    if sym == 1:
        with open(os.path.join(fio.nemoh_root, mesh_file), 'r') as f:
            lines = f.readlines()
        # print(lines[:4])
        lines[0] = lines[0].replace('0', '1')
        # print(lines[:4])
        with open(os.path.join(fio.nemoh_root, mesh_file), 'w') as f:
            f.writelines(lines)

    return nemoh_mesh, mesh_file


def p2f(p_cmplx, pd):
    npanels = pd.ppanels.shape[0]
    f_normal = np.zeros((npanels), dtype=np.complex)
    f = np.zeros((npanels, 3), dtype=np.complex)
    for i, panel in enumerate(pd.ppanels):
        f_normal[i] = p_cmplx[i] * pd.ppanel_areas[i]
        for j in range(3):
            f[i, j] = -f_normal[i] * pd.ppanel_normals[i, j]
    return f

    # sec_panel_force = pressure[ppanels[sec_panel_index,:]]


class PanelData(object):
    def __init__(self, fio, sym=None):
        _pathFile = fio.nemoh_results.joinpath('pressure.{:5d}.dat'.format(1))
        with open(_pathFile, 'r') as f:
            _lines = f.readlines()
        # Number of vertices and number of panels
        ls = _lines[1].split(',')
        self._npoints = int(ls[0][7::])
        self._npanel = int(ls[1][3::])
        self._ppoints = np.asarray([line.split() for line in _lines[2:self._npoints + 2]], dtype='float')[:, 0:3]
        self._ppanels = np.asarray(
                [line.split() for line in _lines[self._npoints + 2:self._npoints + 1 + self._npanel + 1]],
                dtype='int') - 1
        # --------------------------------------
        # Calculate section forces from pressure
        # --------------------------------------
        #
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

        if sym == None or sym == 0:
            pass
        elif sym == 1:
            self._ipanel = np.array([i for i, y in enumerate(self._ppanel_centers[:, 1]) if y < 0])
            self._ppanels[self._ipanel, :] = np.fliplr(self._ppanels[self._ipanel, :])
        else:
            print('sym cannot be {}'.format(sym))
            exit()

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





def diffraction_problem_number(iw, ibeta, Nbeta, Nradiation):
    return 1 + ibeta + iw * (Nbeta + Nradiation)


def radiation_problem_number(iw, iradiation, Nbeta, Nradiation):
    return 1 + iradiation + Nbeta + iw * (Nbeta + Nradiation)
