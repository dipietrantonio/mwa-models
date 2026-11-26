#!/usr/bin/env python3
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
import numpy as np
from astropy.wcs import WCS
from argparse import ArgumentParser
import os
from matplotlib.patches import Rectangle, ConnectionPatch
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

plt.rcParams.update({'font.size': 13})
from astropy.io import fits

def create_wcs_header(frb = 1):
    hdr = fits.Header()

    # ------------------------------------------------------------
    # Basic FITS structure
    # ------------------------------------------------------------
    hdr["NAXIS"]  = 2
    hdr["NAXIS1"] = 256
    hdr["NAXIS2"] = 256

    # ------------------------------------------------------------
    # WCS definition (your values)
    # ------------------------------------------------------------
    hdr["CTYPE1"] = "RA---SIN"
    hdr["CTYPE2"] = "DEC--SIN"

    hdr["CUNIT1"] = "deg"
    hdr["CUNIT2"] = "deg"

    hdr["CRPIX1"] = 129.0
    hdr["CRPIX2"] = 129.0
    if frb == 1:
        hdr["CRVAL1"] = 92.3916173556868
        hdr["CRVAL2"] = -12.9818318398965

        hdr["CDELT1"] = 0.128770968630974
        hdr["CDELT2"] = 0.122176629285385

        # SIN projection parameters
        hdr["PV2_1"] = -0.112925713521288
        hdr["PV2_2"] = -0.247368856436011
    else:
        hdr["CRVAL1"] = 157.944872121852
        hdr["CRVAL2"] = 18.7431502478412 

        hdr["CDELT1"] = 0.133395581144591
        hdr["CDELT2"] = 0.180014434294879

        # SIN projection parameters
        hdr["PV2_1"] = 0.0987080112826317
        hdr["PV2_2"] = -1.01815278888421
       

    return hdr

def draw_rect(images, at, width, height):
    selected_coords = []
    start_row, start_col = at
    line_width = 1
    for row in range(start_row, start_row + height):# + line_width):
        for w in range(line_width):
            selected_coords.append((row, start_col + w))
            selected_coords.append((row, start_col + width + w))
    
    for col in range(start_col, start_col + width):
        for w in range(line_width):
            selected_coords.append((start_row + w, col))
            selected_coords.append((start_row + w + height, col))
    selected_coords.append((start_row + height, start_col + width))
    for p in range(images.shape[0]):
        for i, j in selected_coords:
            images[p, i, j] = 255


def create_panel(fits_files, colormap = "cividis", pixel = (115, 85), frb = 1):
    colormap = "Grays"
    fits_files = sorted(fits_files)
    wcs = WCS(create_wcs_header(frb))
    images = []
    for fname in fits_files:
        print(fname)
        with fits.open(fname) as hdul:
            img = hdul[-1].data.astype(float)
            # If data has extra dimensions, squeeze them
            img = np.squeeze(img)
            images.append(img)

    images = np.array(images) * 2.0 * 0.78
    # Use global scaling so brightness does not jump between frames
    vmin = np.nanmin(images)
    vmax = np.nanmax(images)
    n_images, _, __ = images.shape
    side = 20
    #draw_rect(images, pixel, side, side)

    # ------------------------------------------------------------------
    # Set up figure with two panels
    # ------------------------------------------------------------------
    if n_images < 3:
        n_cols = n_images
        n_rows = 1
    else:
        n_cols = 3
        n_rows = int(n_images / 3)

    fig, axes = plt.subplots(
        n_rows, n_cols, figsize=(10, 5),
        subplot_kw={'projection': wcs}
    )
    fig.supxlabel("Right Ascension", y=0.045)
    fig.supylabel("Declination", x=0.095, ha='right')

    axes = axes.flatten()
    for i, label in zip(range(n_images), ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)"]):
        # Plot images with the SAME normalization
        im1 = axes[i].imshow(images[i], origin="lower",
                        cmap=colormap, vmin=vmin, vmax=vmax)
        # add highlight rectangle around the source
        x0, y0 = 95, 125
        half_h = 10
        highlight_size = 20
        main_box = Rectangle(
            (x0 - half_h, y0 - half_h),
            highlight_size,
            highlight_size,
            fill=False,
            edgecolor="white",
            linewidth=1.2,
        )
        axes[i].add_patch(main_box)
        # Create inset in top-right corner
        axins = inset_axes(axes[i], width="29%", height="29%", loc="upper right", borderpad=0.5,)
        axins.imshow(images[i], origin="lower", interpolation='nearest', cmap=colormap, vmin=vmin, vmax=vmax)
        for spine in axins.spines.values():
            spine.set_edgecolor("white")
            spine.set_linewidth(1.5)
        
        
        zoom_size = 30
        # Zoom limits centred on the FRB
        half_z = zoom_size / 2
        xmin = max(x0 - half_z, -0.5)
        xmax = min(x0 + half_z, 256 - 0.5)
        ymin = max(y0 - half_z, -0.5)
        ymax = min(y0 + half_z, 256 - 0.5)

        axins.set_xlim(xmin, xmax)
        axins.set_ylim(ymin, ymax)

        axins.set_xticks([])
        axins.set_yticks([])
        #axins.set_title("FRB inset", fontsize=9)
        con1 = ConnectionPatch(
            xyA=(x0 - half_h, y0 + half_h),
            coordsA=axes[i].transData,
            xyB=(xmin, ymax),
            coordsB=axins.transData,
            linewidth=1.0,
            color="white",
        )

        con2 = ConnectionPatch(
            xyA=(x0 + half_h, y0 - half_h),
            coordsA=axes[i].transData,
            xyB=(xmax, ymin),
            coordsB=axins.transData,
            linewidth=1.0,
            color="white",
        )

        fig.add_artist(con1)
        fig.add_artist(con2)

        # Name the image (a), (b), (c), ..
        # bbox={
        #                 "facecolor": "black",
        #                 "alpha": 0.5,
        #                 "edgecolor": "none",
        #                 "pad": 3,
        #             }
        axes[i].text(0.03, 0.97, label, transform=axes[i].transAxes, ha="left",
            va="top", fontsize=14, fontweight="bold", color="white",
            )
    
    start_t = 58
    # start_t = 28
    for i, ax in enumerate(axes):
        #Tue Dec 03 2019 frb 1
        # Mar 01, 2020	16:22:54.28 FRB 2
        ax.set_title(f"17:05:55.{start_t + 2 * i:2d} UTC")
        #ax.set_title(f"16:22:54.{start_t + 2 * i:2d} UTC")
        ax.set_ylabel(" ")
        ax.set_xlabel(" ")
    # ------------------------------------------------------------------
    # Shared colorbar (attached to both axes)
    # ------------------------------------------------------------------
    cbar = fig.colorbar(im1, ax=axes, orientation='vertical',  fraction=0.02, pad=0.05)
    cbar.set_label("Pixel value (Jy)")

    plt.show()
   
def animate2(fits_files, colormap = "cividis", pixel = (115, 85), frb = 1):
    fits_files = sorted(fits_files)
    wcs = WCS(create_wcs_header(frb))
    images = []
    for fname in fits_files:
        with fits.open(fname) as hdul:
            img = hdul[-1].data.astype(float)
            # If data has extra dimensions, squeeze them
            img = np.squeeze(img)
            images.append(img)

    images = np.array(images) * 2.0 * 0.78
    # Use global scaling so brightness does not jump between frames
    vmin = np.nanmin(images)
    vmax = np.nanmax(images)
    n_images, _, __ = images.shape
    side = 20

    # ------------------------------------------------------------------
    # Set up figure
    # ------------------------------------------------------------------
    fig = plt.figure(figsize=(6, 6))

    # im = ax.imshow(images[0], origin="lower",, cmap=colormap, vmin=vmin, vmax=vmax)

    start_t = 58

   
    def update(i):
        # Plot images with the SAME normalization
        fig.clear()
        ax = fig.add_subplot(111, projection=wcs)
        ax.set_title(f"17:05:55.{start_t + 2 * i:2d} UTC")
        ax.set_xlabel("Right Ascension")
        ax.set_ylabel("Declination")
        im = ax.imshow(images[i], origin="lower", interpolation='nearest',
                        cmap=colormap, vmin=vmin, vmax=vmax)
        plt.colorbar(im, ax=ax, label="Pixel value (Jy)")
        # add highlight rectangle around the source
        x0, y0 = 95, 125
        half_h = 10
        highlight_size = 20
        main_box = Rectangle(
            (x0 - half_h, y0 - half_h),
            highlight_size,
            highlight_size,
            fill=False,
            edgecolor="white",
            linewidth=1.2,
        )
        ax.add_patch(main_box)
        # Create inset in top-right corner
        axins = inset_axes(ax, width="29%", height="29%", loc="upper right", borderpad=0.5,)
        axins.imshow(images[i], origin="lower", interpolation='nearest', cmap=colormap, vmin=vmin, vmax=vmax)
        for spine in axins.spines.values():
            spine.set_edgecolor("white")
            spine.set_linewidth(1.5)
        
        
        zoom_size = 30
        # Zoom limits centred on the FRB
        half_z = zoom_size / 2
        xmin = max(x0 - half_z, -0.5)
        xmax = min(x0 + half_z, 256 - 0.5)
        ymin = max(y0 - half_z, -0.5)
        ymax = min(y0 + half_z, 256 - 0.5)

        axins.set_xlim(xmin, xmax)
        axins.set_ylim(ymin, ymax)

        axins.set_xticks([])
        axins.set_yticks([])
        #axins.set_title("FRB inset", fontsize=9)
        con1 = ConnectionPatch(
            xyA=(x0 - half_h, y0 + half_h),
            coordsA=ax.transData,
            xyB=(xmin, ymax),
            coordsB=axins.transData,
            linewidth=1.0,
            color="white",
        )

        con2 = ConnectionPatch(
            xyA=(x0 + half_h, y0 - half_h),
            coordsA=ax.transData,
            xyB=(xmax, ymin),
            coordsB=axins.transData,
            linewidth=1.0,
            color="white",
        )

        fig.add_artist(con1)
        fig.add_artist(con2)
    
    

    # ------------------------------------------------------------------
    # Shared colorbar (attached to both axes)
    # ------------------------------------------------------------------
    fps = 3
    frame_interval_ms = (1 / fps) * 1e3 # playback speed
    anim = FuncAnimation(
        fig,
        update,
        frames=len(images),
        interval=frame_interval_ms,
        blit=False
    )
    # Save as GIF
    anim.save("animation.gif", writer=PillowWriter(fps=fps))
   


def animate(fits_files, fps, output_file, colormap = "gray", export_frames = True, format = ".png", dpi = 600, pixel = None, frb = 1):
    fits_files = sorted(fits_files)
    frame_interval_ms = (1 / fps) * 1e3 # playback speed

    wcs = WCS(create_wcs_header(frb))

    images = []
    for fname in fits_files:
        with fits.open(fname) as hdul:
            img = hdul[-1].data.astype(float)
            # If data has extra dimensions, squeeze them
            img = np.squeeze(img)
            images.append(img)

    images = np.array(images)
    # Use global scaling so brightness does not jump between frames
    vmin = np.nanmin(images)
    vmax = np.nanmax(images)

    side = 20
    draw_rect(images, pixel, side, side)
    frame_dir = "frames"
    if export_frames:
        os.makedirs(frame_dir, exist_ok=True)
    # ------------------------------------------------------------------
    # Set up figure
    # ------------------------------------------------------------------
    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(111, projection=wcs)
    im = ax.imshow(images[0], origin="lower", cmap=colormap, vmin=vmin, vmax=vmax)
    title = ax.set_title("Frame 0")
    ax.set_xlabel("Right Ascension")
    ax.set_ylabel("Declination")
    plt.colorbar(im, ax=ax, label="Pixel value (Jy)")

    # ------------------------------------------------------------------
    # Animation update function
    # ------------------------------------------------------------------
    def update(frame):
        im.set_data(images[frame])
        title.set_text(f"Frame {frame}  ({frame * 20} ms)")
        if export_frames:
            fname = os.path.join(frame_dir, f"frame_{frame:04d}.png")
            fig.savefig(fname, dpi=dpi, bbox_inches="tight")
        return [im, title]

    anim = FuncAnimation(
        fig,
        update,
        frames=len(images),
        interval=frame_interval_ms,
        blit=False
    )
    # Save as GIF
    anim.save(output_file, writer=PillowWriter(fps=fps))
    print(f"Saved animation to {output_file}")


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--fps", type=int, default=5, help="Frames (snapshots) per second.")
    parser.add_argument("--output", "-o", type=str, default="animation.gif", help="Output filename.")
    parser.add_argument("--panel", action='store_true', help="Create a panel of images instead.")
    parser.add_argument("--pixel", type=str, default="115,85", help="Pixel coordinates of the top left corner of the rectangle.")
    parser.add_argument("FITS", nargs='+', type=str, help="FITS files representing snapshots to animate. Files are ordered by name.")
    parser.add_argument("--frb", default=1, type=int, help="Index/ID of the FRB to decide which coordinates to display.")
    args = vars(parser.parse_args())
    x, y = args["pixel"].split(',')
    if args["panel"]:
        create_panel(args["FITS"], pixel=(int(x), int(y)), frb=args["frb"])
    else:
        animate2(args["FITS"],) # args["fps"], args["output"], pixel=(int(x), int(y)), frb=args["frb"])


