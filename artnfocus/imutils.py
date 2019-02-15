# Licensed under a 3-clause BSD style license - see LICENSE.rst
# coding=utf-8

from typing import List

import numpy as np

from astropy import stats, wcs
from astropy.nddata import CCDData

import photutils
import ccdproc


def ARTNreduce(filename: str) -> List[CCDData]:
    reduced = []
    hdus = (1, 2)
    for h in hdus:
        im = CCDData.read(filename, hdu=h)
        oscansec = im.header['BIASSEC']
        trimsec = im.header['TRIMSEC']
        im = ccdproc.subtract_overscan(im, fits_section=oscansec, overscan_axis=None)
        im = ccdproc.trim_image(im, fits_section=trimsec)
        reduced.append(im)
    return reduced
