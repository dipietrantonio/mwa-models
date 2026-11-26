#!/usr/bin/env python3
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
import numpy as np
from astropy.wcs import WCS
from argparse import ArgumentParser
import os
plt.rcParams.update({'font.size': 13})
from astropy.io import fits
import matplotlib.image as mpimg


def create_panel(fits_files, png = False, colormap = "gray", pixel = (115, 85), frb = 1):
    fits_files = [x.strip() for x in open(fits_files[0]).readlines()]
    #print(fits_files)
    print(png)
    #fits_files = sorted(fits_files)
    images = []
    wcs_array = []
    if not png:
        for fname in fits_files:
            with fits.open(fname) as hdul:
                img = hdul[-1].data.astype(float)
                wcs_array.append(WCS(hdul[-1].header))
                # If data has extra dimensions, squeeze them
                img = np.squeeze(img)
                images.append(img)
        images = np.array(images)
        vmin = np.nanmin(images)
        vmax = np.nanmax(images)
        n_images, _, __ = images.shape
    else:
        for fname in fits_files:
            img = mpimg.imread(fname)
            images.append(img.squeeze())
        n_images = len(images)

    
    
    # Use global scaling so brightness does not jump between frames
    
    # ------------------------------------------------------------------
    # Set up figure with two panels
    # ------------------------------------------------------------------
    col_size = 6
    if n_images < col_size:
        n_cols = n_images
        n_rows = 1
    else:
        n_cols = col_size
        n_rows = int(n_images / col_size)

    fig = plt.figure(figsize=(10, 10))
    for i in range(n_images):
        if png:
            ax = fig.add_subplot(n_rows, n_cols, i+1)
            im1 = ax.imshow(images[i], origin="lower")
        else:
            ax = fig.add_subplot(n_rows, n_cols, i+1, projection=wcs_array[i])
            im1 = ax.imshow(images[i], origin="lower", cmap=colormap, vmin=vmin, vmax=vmax)
        # ax.set_title(f"17:05:55.{start_t + 2 * i:2d} UTC")
        #ax.set_title(f"16:22:54.{start_t + 2 * i:2d} UTC")
        ax.set_ylabel(" ")
        ax.set_xlabel(" ")
   
    fig.supxlabel("Right Ascension", y=0.045)
    fig.supylabel("Declination", x=0.095, ha='right')
    if not png:
        #fig.tight_layout()
        # ------------------------------------------------------------------
        # Shared colorbar (attached to both axes)
        # ------------------------------------------------------------------
        cbar = fig.colorbar(im1, ax=fig.get_axes(), orientation='vertical',  fraction=0.02, pad=0.05)
        cbar.set_label("Beam strength")

    plt.show()
   


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("FITS", nargs='+', type=str, help="FITS files representing snapshots to animate. Files are ordered by name.")
    parser.add_argument("--png", action='store_true')
    args = vars(parser.parse_args())
    
    create_panel(args["FITS"], png=args["png"])



