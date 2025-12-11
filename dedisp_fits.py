#!/usr/bin/env python3

from math import ceil, floor
from astropy.io import fits
import argparse
import numpy as np

SPEED_OF_LIGHT = 299792458 # m/s
K = 4.15



def read_fits(input_filename):
    my_fits = fits.open(input_filename)
    fits_img = my_fits[0].data
    return fits_img


def to_fits(data, output_filename):
    new_fits = fits.ImageHDU(data)
    new_fits.writeto(output_filename, overwrite=True)


def average_channels(dyspec, avg_factor):
    orig_freq_dim, orig_ts_dim = dyspec.shape
    if orig_freq_dim % avg_factor != 0:
        raise Exception("Averaging factor is not a multiple of the number of channels.")
    new_freq_dim = int(orig_freq_dim / avg_factor)
    new_dyspec = np.ndarray((new_freq_dim, orig_ts_dim))
    for i in range(new_freq_dim):
        new_dyspec[i, :] = np.mean(dyspec[i*avg_factor:i*avg_factor + avg_factor, :], axis=0)
    
    return new_dyspec



def dispersive_delay_s(DM, f_low_ghz, f_high_ghz):
    return K * DM * (f_low_ghz**(-2) - f_high_ghz**(-2)) / 1000



def compute_delay_table(frequencies, dm_list, int_time):
    delay_table = np.zeros((len(dm_list), len(frequencies)), dtype=int)
    top_freq_idx = len(frequencies) - 1
    for dm_idx, dm in enumerate(dm_list):
        delay_table[dm_idx, top_freq_idx] = 0
        for i in range(len(frequencies) - 2, -1, -1):
            delay_table[dm_idx, i] = int(round(dispersive_delay_s(dm, frequencies[i], frequencies[top_freq_idx]) / int_time))
    return delay_table


def incoherent_dedisp(data, delay_table):
    n_channels = data.shape[0]
    delay_steps = [delay_table[0, f + 1] for f in range(n_channels)]
    return np.array([np.roll(row, -x) for row, x in zip(data, delay_steps)])



def compute_frequency_list_ghz(central_freq_mhz, n_channels, channel_width_mhz):
    bottom_freq_mhz = central_freq_mhz - (n_channels * channel_width_mhz) / 2
    return [(bottom_freq_mhz + i * channel_width_mhz) / 1e3 for i in range(n_channels + 1)]



if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--dm", type=float, default=0, help="Dispersion measure used to dedisperse the dynamic spectrum.")
    parser.add_argument("--chan-avg", type=int, default=4, help="Channel averaging factor.")
    parser.add_argument("--freq", type=float, default=153.5, help="Central frequency (in MHz) of the central frequency channel.")
    parser.add_argument("--nchans", type=int, default=768, help="Number of frequency channels.")
    parser.add_argument("--chan-width", type=float, default=0.04, help="Frequency channel width in MHz")
    parser.add_argument("--time-res", type=float, default=0.02, help="Time resolution in seconds.")
    parser.add_argument("--output", type=str, default="out.fits", help="Output FITS filename.")
    parser.add_argument("FITS FILE", type=str, help="FITS file containing the dynamic spectrum.")
    args = vars(parser.parse_args())


    if args["dm"] > 0:
        frequencies = compute_frequency_list_ghz(args["freq"], args["nchans"], args["chan_width"])
        delays = compute_delay_table(frequencies, [args["dm"]], args["time_res"])
        dyspec = incoherent_dedisp(read_fits(args["FITS FILE"]), delays)
    else:
        dyspec = read_fits(args["FITS FILE"])
    
    new_data = average_channels(dyspec, args["chan_avg"])
    to_fits(new_data, args["output"])


    
    