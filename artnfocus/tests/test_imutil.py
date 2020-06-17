# Licensed under a 3-clause BSD style license - see LICENSE.rst
# coding=utf-8

import os
import pkg_resources

from ..imutils import ARTNreduce, sub_background, find_donuts, cutout_donuts
from ..telescope import kuiper_mont4k


def test_get_focus():
    fitsfile = pkg_resources.resource_filename("artnfocus", os.path.join("data", "focusinit.fits"))
    im = ARTNreduce(fitsfile)
    assert(im is not None)
    binning = im.header['BINNING']
    assert(binning == 3)
    im.data = sub_background(im)
    assert(im is not None)
    segm, cat = find_donuts(im)
    assert(len(cat) > 0)
    clean_cat, cutouts, fwhm = cutout_donuts(im, cat)
    assert(len(clean_cat) > 0)
    foc_corr = kuiper_mont4k.simple_focus(pupsize=fwhm, direction="intra", binning=binning)
    assert(foc_corr != 0.)
