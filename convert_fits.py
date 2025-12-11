#!/usr/bin/env python3
import numpy as np
from astropy.io import fits
import sys
import os
import pathlib


def gps_to_unix(gps_time):
    return (315964800 + gps_time - 18)


def covert_to_new_format(blink_fits):
    # Primary HDU (empty, just headers)
    primary_hdu = fits.PrimaryHDU()
    primary_hdu.header['TIME'] = blink_fits[0].header['TIME']
    primary_hdu.header['MILLITIM'] = 0 
    primary_hdu.header['MARKER'] = 0 # TODO: what is this field??
    #primary_hdu.header['OBSID'] = 
    #primary_hdu.header['PROJID'] = 
    
    hdu_list = [primary_hdu]
    for blink_hdu in blink_fits:
        comp_hdu = fits.CompImageHDU(data=blink_hdu.data)
        primary_hdu.header['TIME'] = blink_hdu.header['TIME']
        primary_hdu.header['MILLITIM'] = blink_hdu.header['MILLITIM']
        hdu_list.append(comp_hdu)

 
    # Combine into an HDUList
    hdul = fits.HDUList(hdu_list)
    return hdul


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} <gpubox file 1> <gpubox file 2> ...")
        exit(0)

    output_dir = 'converted'
    p = pathlib.Path(output_dir)
    if not p.exists():
        p.mkdir()
    
    for file in sys.argv[1:]:
        blink_fits = fits.open(file)
        filename = pathlib.Path(file).name
        new_fits = covert_to_new_format(blink_fits)
        new_fits.writeto(f"converted/{filename}", overwrite=True)

    
