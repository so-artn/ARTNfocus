# Licensed under a 3-clause BSD style license - see LICENSE.rst
# coding=utf-8

from typing import List

import numpy as np

from astropy import wcs
from astropy.nddata import CCDData

import photutils
import ccdproc


def ARTNreduce(filename: str) -> CCDData:
    reduced = []
    hdus = (1, 2)
    for h in hdus:
        im = CCDData.read(filename, hdu=h)
        oscansec = im.header['BIASSEC']
        trimsec = im.header['TRIMSEC']
        im = ccdproc.subtract_overscan(im, fits_section=oscansec, overscan_axis=None)
        im = ccdproc.trim_image(im, fits_section=trimsec)
        reduced.append(im)

    # hacky, hard-coded way to make combined image for mont4k at 3x3 binning. result is E to left, N up.
    w = reduced[1].wcs.copy()
    w.array_shape = (1365, 1364)
    blank = np.zeros((1365, 1364))
    blank[:, :682] = reduced[1].data
    blank[:, 682:] = np.fliplr(reduced[0].data)
    stitched = CCDData(blank, wcs=w, unit=reduced[0].unit)
    stitched.data[:, 1088] = (stitched.data[:, 1087] + stitched.data[:, 1089]) / 2.

    return stitched


def sub_background(image: CCDData, filter_size: int = 9, box_size: int = 50) -> np.ndarray:
    bkg_estimator = photutils.MedianBackground()
    bkg = photutils.Background2D(
        image,
        (box_size, box_size),
        filter_size=(filter_size, filter_size),
        bkg_estimator=bkg_estimator
    )
    sub = image.data - bkg.background
    return sub
