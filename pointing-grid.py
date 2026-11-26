#!/usr/bin/env python3

import numpy as np
import argparse

from astropy.coordinates import SkyCoord
import astropy.units as u


def deg_to_radec(ra_deg, dec_deg):
    """
    Convert RA/Dec from decimal degrees to sexagesimal.

    Parameters
    ----------
    ra_deg : float
        Right Ascension in degrees
    dec_deg : float
        Declination in degrees

    Returns
    -------
    ra_hms : str
        RA in hh:mm:ss.ss
    dec_dms : str
        Dec in dd:mm:ss.s
    """
    coord = SkyCoord(ra=ra_deg * u.deg,
                     dec=dec_deg * u.deg,
                     frame="icrs")

    ra_hms = coord.ra.to_string(unit=u.hour,
                                sep=':',
                                precision=2,
                                pad=True)

    dec_dms = coord.dec.to_string(sep=':',
                                  precision=1,
                                  alwayssign=True,
                                  pad=True)

    return ra_hms, dec_dms

def generate_pointing_grid(ra_deg, dec_deg, delta_deg, size):
    """
    Generate a grid of pointings around a central RA/Dec.

    Parameters
    ----------
    ra_deg : float
        Central Right Ascension (degrees)
    dec_deg : float
        Central Declination (degrees)
    delta_deg : float
        Step size in degrees (applied to Dec, and RA corrected by cos(dec))
    n_ra : int
        Number of points in RA direction (must be odd)
    n_dec : int
        Number of points in Dec direction (must be odd)

    Returns
    -------
    coords : list of tuples
        List of (RA_deg, Dec_deg)
    """
    n_ra = size
    n_dec = size
    assert n_ra % 2 == 1 and n_dec % 2 == 1, "Use odd numbers to center the grid"

    ra_list = []
    dec_list = []

    # Convert central dec to radians for cosine correction
    dec0_rad = np.deg2rad(dec_deg)
    cos_dec = np.cos(dec0_rad)

    # Grid indices centered at zero
    ra_offsets = np.arange(-(n_ra//2), n_ra//2 + 1)
    dec_offsets = np.arange(-(n_dec//2), n_dec//2 + 1)

    coords = []

    for i in dec_offsets:
        for j in ra_offsets:
            d_dec = i * delta_deg
            d_ra = (j * delta_deg) / cos_dec

            ra = ra_deg + d_ra
            dec = dec_deg + d_dec

            # Wrap RA into [0, 360)
            ra = ra % 360.0

            coords.append((ra, dec))

    return coords

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--delta", type=float, default=0.1)
    parser.add_argument("--grid-size", type=int, default=3, help="Size of the grid (side)")
    parser.add_argument("COORDS", nargs=2, type=float, help="Central coordinates in degrees.")

    args = vars(parser.parse_args())

    coords = generate_pointing_grid(args["COORDS"][0], args["COORDS"][1], args["delta"], args["grid_size"])
    group_size = 10
    for i in range(int(len(coords) / group_size)):
        with open(f"pointings/pointing_{i:04d}.txt", "w") as fp:
            for c in coords[i*group_size:(i+1)*group_size]:
                ra_hms, dec_dms = deg_to_radec(c[0], c[1])
                fp.write(f"{ra_hms} {dec_dms}\n")
