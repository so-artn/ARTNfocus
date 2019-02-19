# Licensed under a 3-clause BSD style license - see LICENSE.rst
# coding=utf-8

from typing import Tuple

import numpy as np

from astropy import stats
from astropy.nddata import CCDData
from astropy.convolution import Gaussian2DKernel

import photutils
import ccdproc


def ARTNreduce(filename: str) -> CCDData:
    """
    Take a raw FITS image from an ARTN imager (currently only mont4k is supported), perform basic overscan subtraction
    and trimming, and then stitch the images into a single image with the correct sky orientation (E left, N up).
    """
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
    """
    Perform background subtraction using photutils' median background estimator over a 2D mesh.
    """
    bkg_estimator = photutils.MedianBackground()
    bkg = photutils.Background2D(
        image,
        (box_size, box_size),
        filter_size=(filter_size, filter_size),
        bkg_estimator=bkg_estimator
    )
    sub = image.data - bkg.background
    return sub


def find_donuts(
    image: CCDData,
    snr: float = 2.,
    fwhm: float = 10.,
    ksize: int = 15,
    npixels: int = 25
) -> Tuple[photutils.segmentation.core.SegmentationImage, photutils.segmentation.properties.SourceCatalog]:
    """
    Find extended sources in image with default parameters tuned for expected donut size.
    """
    threshold = photutils.detect_threshold(image, snr=snr)
    sigma = fwhm * stats.gaussian_fwhm_to_sigma
    kernel = Gaussian2DKernel(sigma, x_size=ksize, y_size=ksize)
    kernel.normalize()
    segm = photutils.detect_sources(image, threshold, npixels=npixels, filter_kernel=kernel)
    cat = photutils.source_properties(image, segm, wcs=image.wcs)
    return segm, cat
