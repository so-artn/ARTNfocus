# Licensed under a 3-clause BSD style license - see LICENSE.rst
# coding=utf-8

from dataclasses import dataclass

import astropy.units as u


ARCSEC_PER_RADIAN = (1 * u.rad).to(u.arcsec).value

@dataclass
class Telescope:
    diameter: u.Quantity        # Primary mirror diameter
    f_ratio: float              # Focal ratio of optical system
    pix_size: u.Quantity        # Detector pixel size
    obscuration: float = 0.0    # Fractional obscuration due to secondary mirror/baffles
    binning: int = 1            # Detector binning
    counts_per_um: float = 1.0        # Reported focus counts per um of focal point movement

    def __post_init__(self):
        self.focal_length = self.diameter * self.f_ratio  # Focal length of optical system
        self.radius = self.diameter / 2.0  # Radius of primary mirror
        self.nmperrad = self.radius.to(u.nm).value  # nm of wavefront tilt per raadian
        self.nmperasec = self.nmperrad / ARCSEC_PER_RADIAN  # nm of wavefront tilt per arcsecond
        self.plate_scale = ARCSEC_PER_RADIAN * u.arcsec / self.focal_length.to(u.mm)  # Plate scale of focal plane

    def focus_offset(self, foc_delta: float):
        pass

kuiper_mont4k = Telescope(
    diameter = 1.54 * u.m,
    f_ratio = 13.5,
    pix_size = 14 * u.um,
    binning = 3,
    obscuration = 0.266,
    counts_per_um = 600. / 28633.5  # empirically determined...
)
