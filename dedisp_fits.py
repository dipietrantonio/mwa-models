#!/usr/bin/env python3

from math import ceil, floor
from astropy.io import fits
import argparse
import numpy as np
from matplotlib import pyplot as plt

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
    dyspec -= np.nanmedian(dyspec, axis=1)[:,np.newaxis]
    orig_freq_dim, orig_ts_dim = dyspec.shape
    if orig_freq_dim % avg_factor != 0:
        raise Exception("Averaging factor is not a multiple of the number of channels.")
    new_freq_dim = int(orig_freq_dim / avg_factor)
    new_dyspec = np.ndarray((new_freq_dim, orig_ts_dim))
    for i in range(new_freq_dim):
        new_dyspec[i, :] = np.mean(dyspec[i*avg_factor:i*avg_factor + avg_factor, :], axis=0)
    
    return new_dyspec


def compute_time_series(dyspec):
    return average_channels(dyspec, dyspec.shape[0])[0, :]

def compute_iqr(values):
    
    sorted_values = sorted(values)
    q75 = int(len(sorted_values) * 0.75)
    q25 = int(len(sorted_values) * 0.25)
    iqr = sorted_values[q75] - sorted_values[q25]
    stdev = iqr / 1.35
    mean = sorted_values[int(len(sorted_values) / 2)]
    return mean, stdev


def peak_finding(values, snr_threshold = 5):
    peak_idxs = []
    mean, stdev = compute_iqr(values)
    for i, val in enumerate(values):
        estimated_snr = (val - mean) / stdev
        if estimated_snr >= snr_threshold:
            peak_idxs.append(i)
    return peak_idxs


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



def plot_ts_and_dynspec(ds, ts, t, freq, title = None, interp = False, plot_peaks = False):
    """
    ds   : (nchan, nt) dynamic spectrum
    ts   : (nt,) time series
    t    : (nt,) time array [s]
    freq : (nchan,) frequency array [MHz]
    """
    cmap="viridis"

    fig = plt.figure(figsize=(10, 6))
    gs = fig.add_gridspec(
        nrows=2, ncols=1,
        height_ratios=[1, 3],   # TS smaller than dynspec
        hspace=0.05
    )

    ax_ts = fig.add_subplot(gs[0])
    ax_ds = fig.add_subplot(gs[1], sharex=ax_ts)

    # --- Time series ---
    ax_ts.plot(t, ts, color="k", lw=0.8)
    ax_ts.set_ylabel("Intensity")
    ax_ts.tick_params(labelbottom=False)
    ax_ts.grid(alpha=0.3)

    # get peaks
    if plot_peaks:
        peak_idxs = peak_finding(ts)
        peak_vals = [ts[i] for i in peak_idxs]
        ax_ts.scatter([t[i] for i in peak_idxs], peak_vals, color='k')

    # --- Dynamic spectrum ---
    if interp:
        im = ax_ds.imshow(
            ds,
            aspect="auto",
            origin="lower",
            extent=[t[0], t[-1], freq[0], freq[-1]],
            cmap=cmap
        )
    else:
        im = ax_ds.imshow(
        ds,
        aspect="auto",
        origin="lower",
        extent=[t[0], t[-1], freq[0], freq[-1]],
        cmap=cmap,
        interpolation='none'
    )

    ax_ds.set_xlabel("Time (s)")
    ax_ds.set_ylabel("Frequency (MHz)")
    if title is not None:
        fig.suptitle(title)
    #cbar = fig.colorbar(im, ax=ax_ds, pad=0.01)
    #cbar.set_label("Intensity")



def extract_filename_info(filename : str):
    # dynamic_spectrum_00043_00040_dm_397.0_offset_188_candID_6998.fits
    if not filename.startswith("dynamic_spectrum"): raise ValueError()
    components = filename.split('_')
    if len(components) != 10: raise ValueError()
    x, y = int(components[2]), int(components[3])
    dm = float(components[5])
    offset = int(components[7])
    cand_id = components[9][:-5]
    return x, y, dm, offset, cand_id



def analyse_spectrum(dyspec, frequencies, time_res, DM,  plot_title, channel_avg, interp, save_peaks_only):
    if DM > 0:
        delays = compute_delay_table(frequencies, [DM], time_res)
        dyspec = incoherent_dedisp(dyspec, delays)

    dyspec = average_channels(dyspec, channel_avg)
    time_series = compute_time_series(dyspec)
    if save_peaks_only:
        peak_idxs = peak_finding(time_series)
        if len(peak_idxs) == 0: return False

    plot_ts_and_dynspec(dyspec, time_series,
                    [x * time_res for x in range(len(time_series))],
                    [x*1e3 for x in frequencies[:-1]], title=plot_title, interp=interp, plot_peaks=save_peaks_only)
    return True



def process_followup_fits_list(filenames, frequencies, time_res, channel_avg, interp, save_plots, save_peaks_only):
    for filename in filenames:
        x, y, dm, offset, cand_id = extract_filename_info(filename)
        plot_title = f"Candidate {cand_id} - DM {dm} - location ({x}, {y})"
        dyspec = read_fits(filename)
        not_skipped = analyse_spectrum(dyspec, frequencies, time_res, dm, plot_title, channel_avg, interp, save_peaks_only)

        if not_skipped:
            if save_plots:
                plt.savefig(f"{filename}_postprocessed.png", dpi=800)
            else:
                plt.show()
            plt.close()



if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--dm", type=float, default=0, help="Dispersion measure used to dedisperse the dynamic spectrum.")
    parser.add_argument("--chan-avg", type=int, default=4, help="Channel averaging factor.")
    parser.add_argument("--freq", type=float, default=154.237, help="Central frequency (in MHz) of the central frequency channel.")
    parser.add_argument("--nchans", type=int, default=768, help="Number of frequency channels.")
    parser.add_argument("--chan-width", type=float, default=0.04, help="Frequency channel width in MHz")
    parser.add_argument("--time-res", type=float, default=0.02, help="Time resolution in seconds.")
    parser.add_argument("--output", type=str, default="out.fits", help="Output FITS filename.")
    parser.add_argument("--interp", action='store_true', help="Enable interpolation when plotting the dynamic spectrum.")
    parser.add_argument("--fpeaks", action='store_true', help="Only save dynamic spectra with actual peaks (SNR >= 4) in them.")
    parser.add_argument("--save", action='store_true', help="Save plots instead of displaying them.")
    parser.add_argument("FITS FILE", nargs='+', type=str, help="FITS file containing the dynamic spectrum.")

    args = vars(parser.parse_args())
    frequencies = compute_frequency_list_ghz(args["freq"], args["nchans"], args["chan_width"])

    try:
        extract_filename_info(args['FITS FILE'][0])
        process_followup_fits_list(args['FITS FILE'], frequencies, args["time_res"],
                                   args["chan_avg"], args["interp"], args["save"], args["fpeaks"])
    except ValueError:
        # Not the standard followp filename.. use standard processing
        dyspec = read_fits(args["FITS FILE"][0])
        analyse_spectrum(dyspec, frequencies, args["time_res"], args["dm"], args["chan_avg"], args["interp"])
        plt.show()


    
    