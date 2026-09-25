import cmath
from math import pi
import numpy as np
import pytest
from bulbford.dynamics import MAIN2, DISK2, MAIN3, bulb, rho_on_path, orbit


def test_P1_period2_disk_exact():
    """|c_ant − c_root| = 1/2 for the 1/2-bulb: ρ = 4(c+1) = −1 at c = −5/4, root −3/4."""
    b = bulb(MAIN2, 1, 2)
    assert abs(b.ant.c - (-1.25)) < 1e-12
    assert abs(b.root - (-0.75)) < 1e-15
    assert abs(b.ant.rho + 1) < 1e-12 and abs(b.cen.rho) < 1e-12


def test_rho_closed_form_q2():
    """conjugate of λz+z² at λ=−e^ε: ρ = −λ² + 2λ + 4 (coarse path: exercises the collapse guard)."""
    us = np.array([0.3, -0.7 + 0.4j, 1.5j, 2.0])
    got = rho_on_path(MAIN2, 1, 2, us)
    lam = -np.exp(us / 4)
    assert np.allclose(got, -lam ** 2 + 2 * lam + 4, atol=1e-11)


def test_orbit_derivatives_finite_difference():
    c, z, n, h = -0.3 + 0.2j, 0.1 - 0.05j, 5, 1e-6
    o = orbit(MAIN2, c, z, n)
    fd_z = (orbit(MAIN2, c, z + h, n).z - orbit(MAIN2, c, z - h, n).z) / (2 * h)
    fd_c = (orbit(MAIN2, c + h, z, n).z - orbit(MAIN2, c - h, z, n).z) / (2 * h)
    fd_zz = (orbit(MAIN2, c, z + h, n).z_z - orbit(MAIN2, c, z - h, n).z_z) / (2 * h)
    fd_zc = (orbit(MAIN2, c + h, z, n).z_z - orbit(MAIN2, c - h, z, n).z_z) / (2 * h)
    assert abs(o.z_z - fd_z) < 1e-7 and abs(o.z_c - fd_c) < 1e-7
    assert abs(o.z_zz - fd_zz) < 1e-6 and abs(o.z_zc - fd_zc) < 1e-6


@pytest.mark.parametrize("fam", [MAIN2, DISK2, MAIN3])
def test_bulb_boundary_conditions(fam):
    b = bulb(fam, 2, 7)
    assert abs(b.cen.rho) < 1e-10 and abs(b.ant.rho + 1) < 1e-10
    assert 0.9 < b.G_ant < 1.3 and 0.9 < b.G_cen < 1.7   # main3 centre: 1.571 (legacy-confirmed)


def test_reproduces_legacy_q59_values():
    """V5 sample: (x*, G_ant) = (0.0169, 1.0262), (0.4237, 1.1471)."""
    assert abs(bulb(MAIN2, 59 - 1 and 1, 59).G_ant - 1.0262) < 2e-4   # p=1: x* = 1/59
    assert abs(bulb(MAIN2, 26, 59).G_ant - 1.1471) < 2e-4             # p=26: x* = 25/59
