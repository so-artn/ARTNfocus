# Licensed under a 3-clause BSD style license - see LICENSE.rst
# coding=utf-8

from dataclasses import dataclass

import astropy.units as u


ARCSEC_PER_RADIAN = (1 * u.rad).to(u.arcsec).value

@dataclass
class Telescope:
    diameter: u.Quantity
    f_ratio: float
    pix_size: u.Quantity
    obscuration: float = 0.0
    binning: int = 1

    def __post_init__(self):
        self.focal_length = self.diameter * self.f_ratio
        self.radius = self.diameter / 2.0
        self.nmperrad = self.radius.to(u.nm).value
        self.nmperasec = self.nmperrad / ARCSEC_PER_RADIAN
        self.plate_scale = ARCSEC_PER_RADIAN * u.arcsec / self.focal_length.to(u.mm)


kuiper_mont4k = Telescope(
    diameter = 1.54 * u.m,
    f_ratio = 13.5,
    pix_size = 14 * u.um,
    binning = 3,
    obscuration = 0.266
)
