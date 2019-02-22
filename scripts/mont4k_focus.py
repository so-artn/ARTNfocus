#!/usr/bin/env python

import sys
import argparse
import logging
from pathlib import Path

from astropy.io import fits

from artnfocus.telescope import kuiper_mont4k as tel
from artnfocus.imutils import *


log = logging.getLogger("Mont4K Focus")
log.setLevel(logging.DEBUG)

ch = logging.StreamHandler(sys.stdout)
ch.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(levelname)s - %(message)s')
ch.setFormatter(formatter)
log.addHandler(ch)


def get_focus(file, foc_dir):
    im = ARTNreduce(file)
    im.data = sub_background(im)
    segm, cat = find_donuts(im)
    clean_cat, cutouts, fwhm = cutout_donuts(im, cat)
    foc_corr = tel.simple_focus(pupsize=fwhm, direction=foc_dir)
    return foc_corr


def main():
    parser = argparse.ArgumentParser(
        description="Script to analyze Kuiper-Mont4K images of out-of-focus stars and calculate focus corrections."
    )

    parser.add_argument('--intra', '-i', type=str, nargs='+', help="Intra-focal FITS images to analyze")
    parser.add_argument('--extra', '-e', type=str, nargs='+', help="Extra-focal FITS images to analyze")

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

            if focus is not None:
                foc_ret = focus + foc_corr
                log.info(f"{f} -- best_focus={foc_ret}")
                foc_vals.append(foc_ret)
            else:
                log.info(f"{f} -- focus correction={foc_corr}")

            in_foc_corrs.append(foc_corr)

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

            if focus is not None:
                foc_ret = focus + foc_corr
                log.info(f"{f} -- best_focus={foc_ret}")
                foc_vals.append(foc_ret)
            else:
                log.info(f"{f} -- focus correction={foc_corr}")

            ex_foc_corrs.append(foc_corr)

    if not bad_hdr:
        log.info(f"Median best-focus position = {np.median(foc_vals)}")

if __name__ == "__main__":
    main()
