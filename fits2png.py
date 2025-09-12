#!/usr/bin/env python3

from astropy.io import fits
import sys
import matplotlib.pyplot as plt


def dump_image(input, output):
    my_fits = fits.open(input)
    fits_img = my_fits[0].data
    plt.imsave(output, fits_img, cmap='gray')


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} <input fits> <output png>")
    dump_image(sys.argv[1], sys.argv[2])