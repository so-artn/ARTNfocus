# Licensed under a 3-clause BSD style license - see LICENSE.rst
# coding=utf-8

from dataclasses import dataclass

import numpy as np
import astropy.units as u


ARCSEC_PER_RADIAN = (1 * u.rad).to(u.arcsec).value

@dataclass
class Telescope:
    diameter: u.Quantity        # Primary mirror diameter
    f_ratio: float              # Focal ratio of optical system
    pix_size: u.Quantity        # Detector pixel size
    obscuration: float = 0.0    # Fractional obscuration due to secondary mirror/baffles
    binning: int = 1            # Detector binning
    focus_slope: float = 1.0    # change in pupil diameter in pixels per focus readout unit

    def __post_init__(self):
        self.focal_length = self.diameter * self.f_ratio  # Focal length of optical system
        self.radius = self.diameter / 2.0  # Radius of primary mirror
        self.nmperrad = self.radius.to(u.nm).value  # nm of wavefront tilt per raadian
        self.nmperasec = self.nmperrad / ARCSEC_PER_RADIAN  # nm of wavefront tilt per arcsecond
        self.plate_scale = ARCSEC_PER_RADIAN * u.arcsec / self.focal_length.to(u.mm)  # Plate scale of focal plane

    @property
    def convergence_angle(self):
        """
        Angle of convergence of the telescope optics
        """
        return np.arctan2(self.radius, self.focal_length)

    @property
    def offset_slope(self):
        """
        Change in focal point per focus readout unit
        """
        foc_um_slope = self.focus_slope * self.pix_size * self.binning
        offset_slope = 0.5 * foc_um_slope / np.tan(self.convergence_angle)
        return offset_slope

kuiper_mont4k = Telescope(
    diameter = 1.54 * u.m,
    f_ratio = 13.5,
    pix_size = 14 * u.um,
    binning = 3,
    obscuration = 0.266,
    focus_slope = 0.06919
)
