#!/usr/bin/env python

import sys
import argparse
import logging

import numpy as np

from astropy.io import fits

from .telescope import kuiper_mont4k as tel
from .imutils import (sub_background, find_donuts, cutout_donuts, ARTNreduce)

__all__ = ["get_focus", "main"]

log = logging.getLogger("Mont4K Focus")
log.setLevel(logging.DEBUG)

ch = logging.StreamHandler(sys.stdout)
ch.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(levelname)s - %(message)s')
ch.setFormatter(formatter)
log.addHandler(ch)


def get_focus(file, foc_dir):
    """
    Perform data reduction and calculate focus offset from pupil size of detected stars
    """
    im = ARTNreduce(file)
    binning = im.header['BINNING']
    im.data = sub_background(im)
    segm, cat = find_donuts(im)
    clean_cat, cutouts, fwhm = cutout_donuts(im, cat)
    foc_corr = None
    if not np.isnan(fwhm):
        foc_corr = tel.simple_focus(pupsize=fwhm, direction=foc_dir, binning=binning)
    return foc_corr


def main():
    """
    Set up argument handling for command-line use and perform full analysis of pair of images
    """
    parser = argparse.ArgumentParser(
        description="Script to analyze Kuiper-Mont4K images of out-of-focus stars and calculate focus corrections."
    )

    parser.add_argument('--intra', '-i', type=str, nargs='+', help="Intra-focal FITS images to analyze")
    parser.add_argument('--extra', '-e', type=str, nargs='+', help="Extra-focal FITS images to analyze")
    parser.add_argument('--file', '-f', type=str, default=None, help="File to write focus value out to")

    args = parser.parse_args()

    foc_vals = []
    in_foc_corrs = []
    ex_foc_corrs = []
    bad_hdr = False

    if len(args.intra) == 0 and len(args.extra) == 0:
        log.error("Must specify at least one FITS image, either intra or extra focal.")
        return

    if len(args.intra) > 0:
        log.info(f"Averaging focus for {len(args.intra)} intra-focus {'file' if len(args.intra) == 1 else 'files'}:")
        for f in args.intra:
            with fits.open(f) as hlist:
                hdr = hlist[0].header
            if 'FOCUS' in hdr:
                focus = int(hdr['FOCUS'])
                log.debug(f"Got FOCUS={focus} from header for {f}")
            else:
                log.warning(f"No FOCUS in header for {f}")
                focus = None
                bad_hdr = True

            foc_corr = get_focus(f, 'intra')

            if foc_corr is not None:
                if focus is not None:
                    foc_ret = focus + foc_corr
                    log.info(f"{f} -- best_focus={foc_ret}")
                    foc_vals.append(foc_ret)
                else:
                    log.info(f"{f} -- focus correction={foc_corr}")

                in_foc_corrs.append(foc_corr)
            else:
                log.error(f"No focus determined for {f}!")

    if len(args.extra) > 0:
        log.info(f"Averaging focus for {len(args.extra)} extra-focus {'file' if len(args.extra) == 1 else 'files'}:")
        for f in args.extra:
            with fits.open(f) as hlist:
                hdr = hlist[0].header
            if 'FOCUS' in hdr:
                focus = int(hdr['FOCUS'])
                log.debug(f"Got FOCUS={focus} from header for {f}")
            else:
                log.warning(f"No FOCUS in header for {f}")
                focus = None
                bad_hdr = True

            foc_corr = get_focus(f, 'extra')

            if foc_corr is not None:
                if focus is not None:
                    foc_ret = focus + foc_corr
                    log.info(f"{f} -- best_focus={foc_ret}")
                    foc_vals.append(foc_ret)
                else:
                    log.info(f"{f} -- focus correction={foc_corr}")

                ex_foc_corrs.append(foc_corr)
            else:
                if foc_corr is not None:
                    log.error(f"No focus determined for {f}!")

    if not bad_hdr and len(foc_vals) > 1:
        best_foc = int(np.round(np.median(foc_vals)))
        log.info(f"Median best-focus position = {best_foc}")
        if args.file is not None:
            try:
                with open(args.file, "w") as fp:
                    fp.write(f"{best_foc}")
            except Exception as e:
                log.error(f"Error writing focus value to file {args.file}: {e}")
    else:
        if args.file is not None:
            try:
                with open(args.file, "w") as fp:
                    fp.write("N/A")
            except Exception as e:
                log.error(f"Error writing to file {args.file}: {e}")


if __name__ == "__main__":
    main()
