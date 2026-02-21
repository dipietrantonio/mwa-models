#!/usr/bin/env python3


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


if __name__ == "__main__":
    import sys

    ra_deg = float(sys.argv[1])
    dec_deg = float(sys.argv[2])

    ra_hms, dec_dms = deg_to_radec(ra_deg, dec_deg)

    print(f"{ra_hms} {dec_dms}")
