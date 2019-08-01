#!/usr/bin/env python
#  -*- coding: utf-8 -*-

import math

import numpy as np

import meshmagick.mesh_clipper as mc
import meshmagick.mmio as mmio
from meshmagick.mesh import Mesh, Plane


def test_clipper():
    vertices, faces = mmio.load_VTP('meshmagick/tests/data/SEAREV.vtp')
    searev = Mesh(vertices, faces)

    plane = Plane()
    clipper = mc.MeshClipper(searev, plane, assert_closed_boundaries=True, verbose=False)

    for iter in range(50):
        thetax, thetay = np.random.rand(2) * 2 * math.pi
        plane.rotate_normal(thetax, thetay)
        clipper.plane = plane
