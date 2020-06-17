# Licensed under a 3-clause BSD style license - see LICENSE.rst
# coding=utf-8

from ..telescope import kuiper_mont4k


def test_convergence_angle():
    ang = kuiper_mont4k.convergence_angle
    assert(ang > 0.0)


def test_offset_slope():
    slope = kuiper_mont4k.offset_slope
    assert(slope > 0.0)


def test_simple_focus():
    foc = kuiper_mont4k.simple_focus(pupsize=30., direction='intra')
    assert(foc != 0.)
