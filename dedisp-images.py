#!/usr/bin/env python3
from astropy.io import fits
import numpy as np
import argparse
import pathlib
import functools
from math import ceil
SPEED_OF_LIGHT = 299792458 # m/s
K = 4.15

@functools.lru_cache(maxsize=300)
def read_fits(input_filename):
    my_fits = fits.open(input_filename)
    fits_imgs = [my_fits[i].data for i in range(len(my_fits))]
    return fits_imgs


def get_data(input_dir, all_seconds, channel_set, n_timebins, n_channels, t, c):
    lt = int(t / n_timebins)
    lc = int(c / n_channels)
    sec_id = all_seconds[lt]
    chan_id = channel_set[lc]
    fname = f"{input_dir}/start_time_{sec_id}_coarse_{chan_id}_real.fits"
    images = read_fits(fname)
    rt = int(t % n_timebins)
    rc = int(c % n_channels)
    return images[rt * n_channels + rc]


def to_fits(data, output_filename):
    new_fits = fits.ImageHDU(data)
    new_fits.writeto(output_filename, overwrite=True)


def dispersive_delay_s(DM, f_low_ghz, f_high_ghz):
    return K * DM * (f_low_ghz**(-2) - f_high_ghz**(-2)) / 1000


def compute_frequency_list_ghz(central_freq_mhz, n_channels, channel_width_mhz):
    bottom_freq_mhz = central_freq_mhz - (n_channels * channel_width_mhz) / 2
    return [(bottom_freq_mhz + i * channel_width_mhz) / 1e3 for i in range(n_channels + 1)]


def compute_delay_table(frequencies, dm_list, int_time):
    delay_table = np.zeros((len(dm_list), len(frequencies)), dtype=int)
    top_freq_idx = len(frequencies) - 1
    for dm_idx, dm in enumerate(dm_list):
        delay_table[dm_idx, top_freq_idx] = 0
        for i in range(len(frequencies) - 2, -1, -1):
            delay_table[dm_idx, i] = int(round(dispersive_delay_s(dm, frequencies[i], frequencies[top_freq_idx]) / int_time))
    return delay_table


def incoherent_dedisp(input_dir, all_seconds, channel_set, n_channels, n_timebins, img_size, delay_table, time_range, gain, output_dir):
    n_total_timebins = len(all_seconds) * n_timebins
    n_total_channels =  len(channel_set) * n_channels
    delay_steps = [delay_table[0, f + 1] for f in range(n_total_channels)]
    if time_range is None:
        stime = 0
        endtime = n_total_timebins
    else:
        sstime, sendtime = time_range.split(',')
        stime = int(sstime)
        endtime = int(sendtime)

    for t in range(stime, endtime):
        print(f"Processing time bin {t} ..")
        avg_image = np.zeros((img_size, img_size))
        for c in range(n_total_channels):
            dt = t + delay_steps[c]
            avg_image += get_data(input_dir, all_seconds, channel_set, n_timebins, n_channels, int(dt % n_total_timebins), c) * gain
        avg_image /= n_total_channels

        to_fits(avg_image, f"{output_dir}/timebin_{t:05d}.fits")



def image_differencing(filenames, img_size, output_dir):
   
    all_images = np.zeros((len(filenames), img_size, img_size))
    for i, filename in enumerate(filenames):
        data = read_fits(filename)[-1]
        all_images[i, :, :] = data
    median_image = np.median(all_images, axis=0)

    for filename in filenames:
        diff_image = read_fits(filename)[-1] - median_image
        nfname = pathlib.Path(filename).name[:-5]
        to_fits(diff_image, f"{output_dir}/{nfname}_diff.fits")


def retrieve_fits_info(file_list, n_channels, n_timebins):
    fits_imgs = read_fits(file_list[0])
    n_images_per_block = len(fits_imgs)
    if n_images_per_block != (n_channels * n_timebins):
        print("ERROR: Number of images in a FITS file does not correpond to the expected value given the number of channels and time bins requested.")  
    # now see how many seconds are there
    curr_sec = None
    channel_set = []
    all_seconds = []
    for fname in sorted(file_list):
        flname = pathlib.Path(fname).name
        fname_split = flname.split('_')
        sec_id = int(fname_split[2])
        channel_id = int(fname_split[4])
        if curr_sec is None: curr_sec = sec_id
        if curr_sec == sec_id:
            channel_set.append(channel_id)
        else:
            # change in second
            if len(channel_set) != 24:
                print(f"ERROR: second {curr_sec} has only {len(channel_set)} coarse channels")
                raise Exception()
            all_seconds.append(curr_sec)
            curr_sec = sec_id
            channel_set = [channel_id]
    # last second still in the buffer
    if len(channel_set) != 24:
        print(f"ERROR: second {curr_sec} has only {len(channel_set)} coarse channels")
        raise Exception()
    all_seconds.append(curr_sec)
    return sorted(all_seconds), sorted(channel_set)
            





if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--dm", type=float, default=0, help="Dispersion measure used to dedisperse the dynamic spectrum.")
    parser.add_argument("--chan-avg", type=int, default=4, help="Channel averaging factor.")
    parser.add_argument("--time-avg", default=1, type=int, help="Number of contiguous time bins to average.")
    parser.add_argument("--freq", type=float, default=154.237, help="Central frequency (in MHz) of the central frequency channel.")
    parser.add_argument("--nchans", type=int, default=768, help="Number of frequency channels.")
    parser.add_argument("--chan-width", type=float, default=0.04, help="Frequency channel width in MHz")
    parser.add_argument("--time-res", type=float, default=0.02, help="Time resolution in seconds.")
    parser.add_argument("--output", type=str, default="out.fits", help="Output FITS filename.")
    parser.add_argument("--interp", action='store_true', help="Enable interpolation when plotting the dynamic spectrum.")
    parser.add_argument("--fpeaks", action='store_true', help="Only save dynamic spectra with actual peaks (SNR >= 4) in them.")
    parser.add_argument("--save", action='store_true', help="Save plots instead of displaying them.")
    parser.add_argument("--dark", action='store_true', help="Enable dark background mode.")
    parser.add_argument("--id", action='store_true', help="Try to identify the candidate using the Pulsar Survey Scraper")
    parser.add_argument("FITS FILES", nargs='+', type=str, help="FITS file containing the dynamic spectrum.")
    parser.add_argument("--time-range", type=str, default=None, help="Time range to dedisperse. Default: all of it.")
    parser.add_argument("--output-dir", type=str, required=True, help="Output directory where to save output files.")
    parser.add_argument("--dedisp", action='store_true', help="Dedisperse images.")
    parser.add_argument("--diff", action='store_true', help='Run difference imaging.')
    parser.add_argument("--gain", type=float, default=1, help="Calibration factor.")

    args = vars(parser.parse_args())
    frequencies = compute_frequency_list_ghz(args["freq"], args["nchans"], args["chan_width"])
    delays = compute_delay_table(frequencies, [args["dm"]], args["time_res"])
   
    output_dir = args["output_dir"]
    import os
    try:
        os.mkdir(output_dir)
    except:
        pass
    input_dir = pathlib.Path(args["FITS FILES"][0]).parents[0]

    if args["dedisp"]:
        all_seconds, channel_set = retrieve_fits_info(args["FITS FILES"], 32, 50)
        incoherent_dedisp(input_dir, all_seconds, channel_set, 32, 50, 256, delays, args["time_range"], args["gain"], output_dir)

    elif args["diff"]:
        image_differencing(args["FITS FILES"], 256, output_dir)

