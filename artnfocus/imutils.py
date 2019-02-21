# Licensed under a 3-clause BSD style license - see LICENSE.rst
# coding=utf-8

from typing import Tuple, List

import numpy as np

import astropy.units as u
from astropy import stats
from astropy.nddata import CCDData, Cutout2D
from astropy.convolution import Gaussian2DKernel

import photutils
import ccdproc

from astroscrappy import detect_cosmics


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
    cr_mask, clean_data = detect_cosmics(blank, sigclip=5., niter=5, cleantype='medmask', psffwhm=10.)
    stitched = CCDData(clean_data, wcs=w, unit=reduced[0].unit)

    # nuke the bad column
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
    snr: float = 2.,    # Threshold SNR for segmentation
    fwhm: float = 10.,  # Kernel FWHM for segmentation
    ksize: int = 15,    # Kernel size
    npixels: int = 25   # Number of connected pixels required to be considered a source
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


def cutout_donuts(
    image: CCDData,
    cat: photutils.segmentation.properties.SourceCatalog,
    size: int = 100,                # cutout size
    buffer: int = 5,                # edge buffer
    saturation: int = 60000,        # when saturation is reached
    min_area: float = 500.,         # minimum source size
    max_ellipticity: float = 0.15,  # maximum ellipticity
    nsigma: float = 10.0            # noise threshold for max value
) -> Tuple[photutils.segmentation.properties.SourceCatalog, List[Cutout2D], float]:
    """
    Create cleaned list of donuts and return list of cutouts for the donuts and their mean width
    """
    cutouts = []
    good = []
    fwhms = []
    minpos = u.pix * (int(size / 2.) + buffer)  # give a bit of a buffer at the edge
    maxpos = image.shape[0] * u.pix - minpos
    mean, median, stddev = stats.sigma_clipped_stats(image, sigma=2.0, maxiters=None)
    for s in cat:
        valid_pos = s.xcentroid > minpos and s.xcentroid < maxpos and s.ycentroid > minpos and s.ycentroid < maxpos
        unsaturated = s.max_value < saturation
        big_enough = s.area > min_area * u.pix ** 2
        is_round = s.ellipticity < max_ellipticity
        is_bright = s.max_value > mean + nsigma * stddev
        if valid_pos and unsaturated and big_enough and is_round and is_bright:
            cutout = Cutout2D(image, (s.xcentroid.value, s.ycentroid.value), (size, size), wcs=image.wcs)
            y_sum = np.sum(cutout.data, axis=0)
            x_sum = np.sum(cutout.data, axis=1)
            y_fwhm = len(np.where(y_sum > np.max(y_sum)/2.)[0])
            x_fwhm = len(np.where(x_sum > np.max(x_sum)/2.)[0])
            fwhm = 0.5 * (y_fwhm + x_fwhm)
            cutouts.append(cutout)
            good.append(s)
            fwhms.append(fwhm)
    fwhm = np.median(fwhms)
    clean_cat = photutils.SourceCatalog(good)
    return (clean_cat, cutouts, fwhm)
