# Licensed under a 3-clause BSD style license - see LICENSE.rst
# coding=utf-8

from astropy import stats, wcs
from astropy.nddata import CCDData
import photutils
import ccdproc


def imreduce(filename: str):
