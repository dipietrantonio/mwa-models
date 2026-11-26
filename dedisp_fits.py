#!/usr/bin/env python3

from math import ceil, floor
from astropy.io import fits
import argparse
import numpy as np
from matplotlib import pyplot as plt
import struct
from statistics import stdev, median, mean
from scipy.stats import norm, skew, kurtosis, kstest
import os
import shutil

SPEED_OF_LIGHT = 299792458 # m/s
K = 4.148808 # 4.15


def extract_filename_info(filename : str):
    # dynamic_spectrum_00043_00040_dm_397.0_offset_188_candID_6998.fits
    if not filename.startswith("dynamic_spectrum"): raise ValueError()
    components = filename.split('_')

    if len(components) == 10:
        x, y = int(components[2]), int(components[3])
        dm = float(components[5])
        offset = int(components[7])
        cand_id = components[9][:-5]
        return x, y, dm, offset, cand_id
    
    elif len(components) == 16:
        x, y = int(components[2]), int(components[3])
        ra = float(components[5])
        dec = float(components[7])
        snr = float(components[9])
        dm = float(components[11])
        offset = int(components[13])
        cand_id = components[15][:-5]
        return x, y, ra, dec, snr, dm, offset, cand_id

    else:
        raise ValueError("Error while parsing the filename.")




def read_fits(input_filename):
    my_fits = fits.open(input_filename)
    fits_img = my_fits[0].data
    return fits_img



def read_filterbank(filename):

    string_keys = {
        "rawdatafile", "source_name"
    }

    int_keys = {
        "machine_id", "telescope_id", "data_type",
        "barycentric", "pulsarcentric", "nbits",
        "nsamples", "nchans", "nifs", "nbeams",
        "ibeam"
    }

    double_keys = {
        "az_start", "za_start", "src_raj", "src_dej",
        "tstart", "tsamp", "fch1", "foff",
        "refdm", "period"
    }

    def read_string(f):
        length = struct.unpack("<i", f.read(4))[0]

        if not 0 < length < 65536:
            raise ValueError("Invalid header string length")

        return f.read(length).decode("ascii")

    header = {}

    with open(filename, "rb") as f:

        if read_string(f) != "HEADER_START":
            raise ValueError("Not a SIGPROC filterbank")

        while True:

            key = read_string(f)

            if key == "HEADER_END":
                break

            if key in string_keys:
                value = read_string(f)

            elif key in int_keys:
                value = struct.unpack("<i", f.read(4))[0]

            elif key in double_keys:
                value = struct.unpack("<d", f.read(8))[0]

            elif key == "signed":
                value = struct.unpack("<b", f.read(1))[0]

            else:
                raise ValueError(f"Unknown header key: {key}")

            header[key] = value

        offset = f.tell()

    nchans = header["nchans"]
    nbits = header["nbits"]
    nifs = header.get("nifs", 1)
    signed = header.get("signed", 0)

    dtypes = {
        8: "i1" if signed else "u1",
        16: "<i2" if signed else "<u2",
        32: "<f4"
    }

    if nbits not in dtypes:
        raise ValueError(f"Unsupported nbits: {nbits}")

    dtype = np.dtype(dtypes[nbits])

    bytes_per_sample = nchans * nifs * dtype.itemsize
    payload_size = os.path.getsize(filename) - offset

    nsamples, remainder = divmod(
        payload_size, bytes_per_sample
    )

    if remainder:
        raise ValueError("Incomplete spectrum at end of file")

    data = np.memmap(
        filename,
        dtype=dtype,
        mode="r",
        offset=offset,
        shape=(nsamples, nifs, nchans)
    )

    header["nsamples"] = nsamples
    header["header_size"] = offset
    return header, data


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



def average_timesteps(dyspec, avg_factor, offset = 0):
    orig_freq_dim, orig_ts_dim = dyspec.shape
    new_ts_dim = int(ceil(orig_ts_dim / avg_factor))
    new_dyspec = np.ndarray((orig_freq_dim, new_ts_dim))
    # The offset is used to align the averaging operation with the start of a pulse,
    # so that the signal does not get diluted with noise. Ideally,
    # offset is the start of the pulse, averaging factor is the pulse width,
    # this will get maximum SNR.
    # Initial alignment is the initial number of bins to average and so that the subsequent
    # `avg_factor` bin groups can align with the pulse. If 0, bin groups are already aligned
    # with the pulse.
    initial_alignment = int(offset % avg_factor)
    dyspec = dyspec[:, initial_alignment:]
    for i in range(new_ts_dim):
        new_dyspec[:, i] = np.mean(dyspec[:, i*avg_factor:i*avg_factor + avg_factor], axis=1)
    return new_dyspec


def compute_time_series(dyspec):
    return average_channels(dyspec, dyspec.shape[0])[0, :]


def study_time_series_noise(time_series : np.ndarray, time_res):
    """
    Compute the noise level on progressively longer prefixes of the time series.
    The is a time series describing the noise level as a function of time series length.
    """
    # Temp: load marcins time series instead
    
    # time_series = np.array([float(x[1]) for x in [line.split() for line in open("/Users/dip008/Downloads/timeseries_timeres20.0ms_NO-FRB.txt").readlines()]])
    # print(time_series)
    # minimum length is 1 second, as is the incremental step size
    min_length = 2 #int(1 / time_res)
    step_size = 1 #int(1 / time_res)
    n_points = time_series.size
    n_steps = int((n_points - min_length + step_size - 1) / step_size + 1)
    noise_levels = np.zeros((n_steps,))
    mean_levels = np.zeros((n_steps,))
    for i in range(n_steps):
        i_e = min_length + i * step_size
        sub_series = time_series[0 : i_e]
        mean_levels[i], noise_levels[i] = compute_iqr(sub_series)
    return noise_levels

def study_time_series_noise2(time_series : np.ndarray, time_res):
    """
    Compute the noise level on progressively longer prefixes of the time series.
    The is a time series describing the noise level as a function of time series length.
    """
    # minimum length is 1 second, as is the incremental step size
    length = int(2 / time_res)
    n_points = time_series.size
    n_steps = n_points - length + 1
    noise_levels = np.zeros((n_steps,))
    mean_levels = np.zeros((n_steps,))
    for i in range(n_steps):
        sub_series = time_series[i : i + length]
        mean_levels[i], noise_levels[i] = compute_iqr(sub_series)
    return noise_levels


def fit_gaussian_to_timeseries(fig, time_series: np.ndarray, title):
    mu, sigma = norm.fit(time_series)
    N = len(time_series)
    sk = skew(time_series)
    kurt = kurtosis(time_series)
    ks_stat, ks_p = kstest(time_series, 'norm', args=(mu, sigma))
    gs = fig.add_gridspec(1, 2, width_ratios = [3, 1])
    ax_hist = fig.add_subplot(gs[0])
    ax_text = fig.add_subplot(gs[1])

    count, bins, _ = ax_hist.hist(time_series, bins=60, density=True, alpha=0.6)
    x = np.linspace(bins[0], bins[-1], 1000)
    ax_hist.plot(x, norm.pdf(x, mu, sigma), 'r', lw=2)

    ax_hist.set_title(title)
    ax_hist.set_xlabel("Flux")
    ax_hist.set_ylabel("Probability Density")

    stats_text = (
        f"N = {N}\n\n"
        f"Mean (μ) = {mu:.4g}\n"
        f"Std (σ) = {sigma:.4g}\n\n"
        f"Skew = {sk:.4g}\n"
        f"Excess Kurtosis = {kurt:.4g}\n\n"
        f"KS stat = {ks_stat:.4g}\n"
        f"KS p-value = {ks_p:.4g}"
    )

    ax_text.text(0.05, 0.95, stats_text, va='top', fontsize=11)
    ax_text.axis('off')
    fig.tight_layout()



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
            print(f"Peak no. {i+1}: est. SNR = {estimated_snr}, value = {val}, mean = {mean}, stdev = {stdev}")    
            peak_idxs.append((i, estimated_snr))
    return peak_idxs



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



def incoherent_dedisp(data, delay_table):
    n_channels = data.shape[0]
    dedispersed_dyspectra = []
    for dm_idx in range(delay_table.shape[0]):
        delay_steps = [delay_table[dm_idx, f + 1] for f in range(n_channels)]
        dedispersed_dyspectra.append(np.array([np.roll(row, -x) for row, x in zip(data, delay_steps)]))
    return dedispersed_dyspectra



def transform_spectrum(dyspec, frequencies, time_res, DM, channel_avg, time_avg, time_avg_offset, disable_norm = False, gain = 1):
    dm_time = None
    dm_list = None
    #plt.imshow(dyspec, origin='lower')
    #plt.show()
    if DM > 0:
        dm_radius = 10
        dm_list = sorted([DM - i for i in range(1, dm_radius + 1)] + [DM + i for i in range(dm_radius)])
        delays = compute_delay_table(frequencies, dm_list, time_res)
        dyspec_list = incoherent_dedisp(dyspec, delays)
        dyspec = dyspec_list[dm_radius]
        dm_time = np.zeros((len(dm_list), dyspec.shape[1]))
        for i, ds in enumerate(dyspec_list):
            dm_time[i, :] = np.mean(ds, axis=0)

    dyspec *= gain
    #medians = np.nanmedian(dyspec, axis=1)[:,np.newaxis]
    #dyspec -= medians

    if not disable_norm:
        n_freqs = dyspec.shape[0]
        stds = np.zeros((n_freqs,))
        for i in range(n_freqs):
            _, stdev = compute_iqr(dyspec[i, :])
            stds[i] = stdev
        stds = stds[:, np.newaxis]
        dyspec /= stds
    
    dyspec = average_channels(dyspec, channel_avg)
    if time_avg > 1:
        dyspec = average_timesteps(dyspec, time_avg, floor(time_avg_offset / time_res))

    time_series = compute_time_series(dyspec)
    median_series = np.median(dyspec, axis=0)
    peak_idxs = peak_finding(time_series)
    #print(peak_idxs)
    # print(time_series[peak_idxs[0][0]-4: peak_idxs[0][0]+15])
    # NOTE: test to compute SNR on smaller range around peak
    # max_peak = sorted(peak_idxs, key=lambda x: x[1])[-1][0]
    # width = int(3 / 0.08)
    # print("SNR on SMALLER series")
    # peak_finding(time_series[max_peak-width:max_peak+width])

    return dyspec, dm_time, dm_list, time_series, median_series, peak_idxs



def transform_spectrum_split(n_splits, dyspec, frequencies, time_res, DM, channel_avg, time_avg, time_avg_offset, disable_norm = False, gain = 1):
    fs_size = int((len(frequencies) + n_splits - 1) / n_splits)
    ds_size = int((dyspec.shape[0] + n_splits - 1) / n_splits)
    frequencies_split = []
    for i in range(n_splits):
        start_idex = i * fs_size - i
        end_index = start_idex + fs_size
        frequencies_split.append(frequencies[start_idex:end_index])
        # print(ds_size)
    dyspec_split = [dyspec[i * ds_size:(i + 1) * ds_size, 0:6000] for i in range(n_splits)]
    results = []
    for i in range(n_splits):
        #print(len(dyspec_split[i]), len(frequencies_split[i]))
        #print(frequencies_split[i][0], frequencies_split[i][-1])
        results.append(transform_spectrum(dyspec_split[i], frequencies_split[i], time_res, DM, channel_avg, time_avg, time_avg_offset, disable_norm, gain))
    
    print("OK")
    fig = plt.figure()
    split_panel = fig.add_gridspec(nrows=3, ncols=1)
    for i, r in enumerate(results):
        ax = fig.add_subplot(split_panel[2 - i])
        ax.plot(r[3])
        #ax.imshow(r[0], origin='lower')
    fig.tight_layout()
    plt.show()
    print("Done")
    return results
    



def plot_ts_and_dynspec(fig, ds, dm_time, dm_list, ts, median, peak_idxs, t, freq, title = None, interp = False, paper_mode = False, timeseries_only = False):
    """
    ds   : (nchan, nt) dynamic spectrum
    ts   : (nt,) time series
    t    : (nt,) time array [s]
    freq : (nchan,) frequency array [MHz]
    """
    cmap="viridis"
    if paper_mode:
        # do not plot median time series.
        gs = fig.add_gridspec(
            nrows=3, ncols=1,
            height_ratios=[1, 3, 1],   # TS smaller than dynspec
            #width_ratios=[4, 1],
            hspace=0.05,
            #wspace=0.05
        )
        ax_ts = fig.add_subplot(gs[0, 0])
        ax_ds = fig.add_subplot(gs[1, 0], sharex=ax_ts)
        ax_dm_time = fig.add_subplot(gs[2, 0], sharex=ax_ts)
        #ax_median_ds = fig.add_subplot(gs[1, 1], sharey=ax_ds)
    elif timeseries_only:
        gs = fig.add_gridspec(
            nrows=1, ncols=1,
            #height_ratios=[1, 3, 1],   # TS smaller than dynspec
            #width_ratios=[4, 1],
            #hspace=0.05,
            #wspace=0.05
        )
        ax_ts = fig.add_subplot(gs[0, 0])
    else:
        gs = fig.add_gridspec(
            nrows=4, ncols=1,
            height_ratios=[1, 1, 3, 1],   # TS smaller than dynspec
            #width_ratios=[4, 1],
            hspace=0.05
        )
        ax_ts = fig.add_subplot(gs[0])
        ax_median_ts = fig.add_subplot(gs[1], sharex=ax_ts)
        ax_ds = fig.add_subplot(gs[2], sharex=ax_ts)
        ax_dm_time = fig.add_subplot(gs[3], sharex=ax_ts)

    # if paper_mode:
    #     median = np.max(ds, axis=1)
    #     # if len(peak_idxs) > 0:
    #     #     i = sorted(peak_idxs, key=lambda x: x[1], reverse=True)[0][0]
    #     #     peak_spectrum = ds[:, i+300].flatten()
    #     #     print(len(peak_spectrum))
    #     #     ax_median_ds.plot(peak_spectrum, np.arange(len(peak_spectrum)), c='r', lw=0.8)
    #     ax_median_ds.plot(median, np.arange(len(median)), lw=0.8)
    #     ax_median_ds.set_xlabel("Median channel intensity")
    #     ax_median_ds.tick_params(labelleft=False)
   
    # --- Time series ---
    ax_ts.plot(t, ts, lw=0.8)
    ax_ts.set_ylabel(f"Mean intensity{' (Jy)' if paper_mode else ''}")
    if not timeseries_only:
        ax_ts.tick_params(labelbottom=False)
    if not paper_mode: ax_ts.grid(alpha=0.3)

    # get peaks
    if len(peak_idxs) > 0 and not paper_mode:
        peak_vals = [ts[i[0]] for i in peak_idxs]
        ax_ts.scatter([t[i[0]] for i in peak_idxs], peak_vals)

    # ---- Median time series ---
    if not paper_mode and not timeseries_only:
        ax_median_ts.plot(t, median, lw=0.8)
        ax_median_ts.set_ylabel("Median intensity")
        ax_median_ts.grid(alpha=0.3)
        ax_median_ts.tick_params(labelbottom=False)
   
    # --- Dynamic spectrum ---
    if not timeseries_only:
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
        ax_ds.tick_params(labelbottom=False)
        ax_ds.set_ylabel("Frequency (MHz)")
        
        # DM TIME plot
        ax_dm_time.imshow(
            dm_time,
            aspect="auto",
            origin="lower",
            extent=[t[0], t[-1], dm_list[0], dm_list[-1]],
            cmap=cmap,
            interpolation='none'
        )
        ax_dm_time.set_xlabel("Time (s)")
        ax_dm_time.set_ylabel("DM (pc$\,$cm$^{-3}$)")

    if title is not None and not paper_mode:
        fig.suptitle(title)
    #cbar = fig.colorbar(im, ax=ax_ds, pad=0.01)
    #cbar.set_label("Intensity")


def find_time_index(time_range, time):
    for i, v in enumerate(time_range):
        if v >= time:
            return i
    raise Exception()

def centre_on_brightest_burst(dyspec, dm_time, time_series, median_series, peak_idxs, dt):
    """
    Centre arrays on the brightest burst.

    Returns:
        dyspec_c        : centred/cropped dynamic spectrum
        dm_time_c       : centred/cropped DM-time plane
        time_series_c   : centred/cropped time series
        median_series_c : centred/cropped median time series
        peak_idxs_c     : peak indices shifted to the new cropped frame
        time_range_c    : time axis with brightest burst at t=0
        peak_idx_c      : index of brightest burst in cropped arrays
    """

    nt = len(time_series)

    # If peak_finding found peaks, use the brightest one in the time series.
    # Otherwise fall back to the global maximum.
    if len(peak_idxs) > 0:
        peak_idx = max(peak_idxs, key=lambda p: time_series[p[0]])[0]
    else:
        peak_idx = int(np.nanargmax(time_series))

    # Crop symmetrically around the peak so that the peak is visually centred.
    n_left = peak_idx
    n_right = nt - peak_idx - 1
    half_width = int(min(n_left, n_right) / 2)

    i0 = peak_idx - half_width
    i1 = peak_idx + half_width + 1

    dyspec_c = dyspec[:, i0:i1]
    time_series_c = time_series[i0:i1]
    median_series_c = median_series[i0:i1]

    if dm_time is not None:
        dm_time_c = dm_time[:, i0:i1]
    else:
        dm_time_c = None

    # Shift peak indices after cropping
    peak_idxs_c = []
    for p in peak_idxs:
        old_idx = p[0]
        if i0 <= old_idx < i1:
            new_peak = list(p)
            new_peak[0] = old_idx - i0
            peak_idxs_c.append(tuple(new_peak))

    peak_idx_c = peak_idx - i0

    # Time axis centred on brightest burst
    n_c = len(time_series_c)
    time_range_c = (np.arange(n_c) - peak_idx_c) * dt

    return (
        dyspec_c,
        dm_time_c,
        time_series_c,
        median_series_c,
        peak_idxs_c,
        time_range_c,
        peak_idx_c,
    )


def plot_spectrum(main_fig, mean_fig, median_fig, dyspec, time_offset, frequencies, time_res, dm, channel_avg,
                  time_avg,  time_avg_offset, plot_title, interp, disable_norm, gain, fit_mean, fit_median,
                  paper_mode = False, timeseries_only = False):
    #print("Called")
    #dyspec, dm_time, dm_list, time_series, median_series, peak_idxs = transform_spectrum_split(3, dyspec, frequencies, time_res, dm, channel_avg, time_avg, time_avg_offset - time_offset, disable_norm, gain,)
    dyspec, dm_time, dm_list, time_series, median_series, peak_idxs = transform_spectrum(dyspec, frequencies, time_res, dm, channel_avg, time_avg, time_avg_offset - time_offset, disable_norm, gain,)

    dt = time_res * time_avg if time_avg > 1 else time_res

    # (dyspec, dm_time, time_series, median_series, peak_idxs, time_range, peak_idx ) = centre_on_brightest_burst(
    #     dyspec,
    #     dm_time,
    #     time_series,
    #     median_series,
    #     peak_idxs,
    #     dt,
    # )
    time_range = [time_offset + x * time_res * time_avg for x in range(len(time_series))]
    freq_range = [x*1e3 for x in frequencies[:-1]]
    plot_ts_and_dynspec(main_fig, dyspec, dm_time, dm_list, time_series, median_series, peak_idxs,
                    time_range, freq_range, title=plot_title, interp=interp, paper_mode=paper_mode, timeseries_only=timeseries_only)

    def parse_time_range(range_spec):
        if range_spec == "": return 0, len(time_range)
        comp = range_spec.split(',')
        if len(comp) == 1:
            low = find_time_index(time_range, float(comp[0]))
            high = len(time_range)
        else:
            low = find_time_index(time_range, float(comp[0]))
            high = find_time_index(time_range, float(comp[1]))
        return low, high

    if fit_mean is not None and mean_fig is not None:
        try:
            low, high = parse_time_range(fit_mean)
        except:
            print("WARNING: could not parse time range.")
            low = 0
            high = len(time_series)
        fit_gaussian_to_timeseries(mean_fig, time_series[low:high], "Mean Flux Density")
    
    if fit_median is not None and median_fig is not None:
        try:
            low, high = parse_time_range(fit_median)
        except:
            print("WARNING: could not parse time range.")
            low = 0
            high = len(median_series)
            
        fit_gaussian_to_timeseries(median_fig, median_series[low:high], "Median Flux Density")




def process_followup_fits_list(filenames, frequencies, time_res, channel_avg, time_avg, time_avg_offset, args_dm, interp, disable_norm, gain,
                               save_plots, identities, fit_mean, fit_median, paper_mode, timeseries_only):
    
    main_fig = plt.figure(figsize=(10, 6))
    mean_fig = None
    median_fig = None
    if fit_mean is not None:
        mean_fig = plt.figure(figsize=(10, 4))
    if fit_median is not None:
        median_fig = plt.figure(figsize=(10, 4))
    current_file_idx = 0

    def process_fits(idx, filename):
        filename_info = extract_filename_info(filename)
        if len(filename_info) == 5: # old format
            x, y, dm, offset, cand_id = filename_info
            plot_title = f"[{idx + 1} / {len(filenames)}] Candidate {cand_id} - DM {dm} - location ({x}, {y})"
        else:
            x, y, ra, dec, snr, dm, offset, cand_id = filename_info
            ident = "ID Unknown"
            if identities is not None:
                ident = identities[filename]
            plot_title = f"[{idx + 1} / {len(filenames)}] Candidate {cand_id} - SNR {snr} - DM {dm} - Location (RA = {ra}, DEC = {dec}) \n\n{ident}"

        dyspec = read_fits(filename)
        if args_dm > 0:
            dm = args_dm
        plot_spectrum(main_fig, mean_fig, median_fig, dyspec, int(offset), frequencies, time_res, dm, channel_avg,
                      time_avg, time_avg_offset, plot_title, interp, disable_norm, gain, fit_mean, fit_median, paper_mode, timeseries_only)
    
    if save_plots:
        for idx, file in enumerate(filenames):
            main_fig.clear()
            process_fits(idx, file)
            main_fig.savefig(f"{file}_postprocessed.png", dpi=800)
        return
    

    def on_keypress(event, main_fig):
        nonlocal current_file_idx
        if event.key in ["up", "pageup", "left"]:
            # go one slide backwards
            current_file_idx -= 1
            if current_file_idx == -1:
                current_file_idx = len(filenames) - 1
        elif event.key in ["down", "pagedown", "right"]:
            # go forward
            current_file_idx = (current_file_idx + 1) % len(filenames)
        else:
            # command not recognised: no nothing
            return
        
        main_fig.clear()
        if mean_fig is not None: mean_fig.clear()
        if median_fig is not None: median_fig.clear()

        process_fits(current_file_idx, filenames[current_file_idx])
        main_fig.canvas.draw()
        if mean_fig is not None: mean_fig.canvas.draw()
        if median_fig is not None: median_fig.canvas.draw()

    
    main_fig.canvas.mpl_connect('key_press_event', lambda event: on_keypress(event, main_fig))

    # Display the initial candidate
    process_fits(0, filenames[0])
    plt.show()



def identify_candidate(filename):
    info = extract_filename_info(filename)
    ident = "ID Unknown"
    if len(info) > 5:
        x, y, ra, dec, snr, dm, offset, cand_id = info
        import scraper
        cands = scraper.query_scraper(ra, dec)
        cands_with_close_dm = cands.loc[(cands['DM'] >= dm - 1) & (cands['DM'] <= dm + 1)]
        if len(cands_with_close_dm) > 0:
            first_cand = cands_with_close_dm.iloc[0]
            ident = f"ID - {first_cand['Name']} - DM {first_cand['DM']} pc cm^-3 - Distance {first_cand['Distance']:.2f} deg."
    return ident


def identify_candidates(filenames):
    print("Running candidates identification...")
    from concurrent.futures import ThreadPoolExecutor
    # Use a ThreadPoolExecutor to manage the threads
    with ThreadPoolExecutor(max_workers=16) as executor: # You can adjust max_workers
        # The map method automatically distributes the urls to the fetch_url function
        results = dict(zip(filenames, list(executor.map(identify_candidate, filenames))))
        all_values = set(results.values())
        if 'ID Unknown' in all_values: all_values.remove("ID Unknown")
    print(f"Candidates identification done. {len(all_values)} known sources.")
    return results



if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--dm", type=float, default=0, help="Dispersion measure used to dedisperse the dynamic spectrum.")
    parser.add_argument("--chan-avg", type=int, default=4, help="Channel averaging factor.")
    parser.add_argument("--time-avg", default=1, type=int, help="Number of contiguous time bins to average.")
    parser.add_argument("--time-avg-offset", default=0.0, type=float, help="Offset in seconds to the position where to align the time averaging operation.")
    parser.add_argument("--freq", type=float, default=154.237, help="Central frequency (in MHz) of the central frequency channel.")
    parser.add_argument("--nchans", type=int, default=768, help="Number of frequency channels.")
    parser.add_argument("--chan-width", type=float, default=0.04, help="Frequency channel width in MHz")
    parser.add_argument("--time-res", type=float, default=0.02, help="Time resolution in seconds.")
    parser.add_argument("--output", type=str, default="out.fits", help="Output FITS filename.")
    parser.add_argument("--interp", action='store_true', help="Enable interpolation when plotting the dynamic spectrum.")
    parser.add_argument("--fpeaks", action='store_true', help="Only save dynamic spectra with actual peaks (SNR >= 5) in them.")
    parser.add_argument("--snoise", action='store_true', help="Study the noise of the time series.")
    parser.add_argument("--save", action='store_true', help="Save plots instead of displaying them.")
    parser.add_argument("--paper", action="store_true", help="Enable paper-ready visualisation mode.")
    parser.add_argument("--dark", action='store_true', help="Enable dark background mode.")
    parser.add_argument("--id", action='store_true', help="Try to identify the candidate using the Pulsar Survey Scraper")
    parser.add_argument("--max-cands", type=int, default=-1, help="Maximum number of candidates a dynamic spectrum can contain not to be discarded in filtering. This helps with RFI.")
    parser.add_argument("--gain", type=float, default=1, help="Apply calibration gain to the dynamic spectrum, before averaging.")
    parser.add_argument("--disable-norm", action='store_true', help="Disable median subtraction used to normalise power across channels.")
    parser.add_argument("--timeseries-only", action='store_true', help="Display only the dedispersed time series.")
    parser.add_argument("--fit-median", nargs='?', type=str, const="", default=None, help="Fit a Gaussian model to the median time series.\n"
                        "Optionally, a time range 'start,end' in seconds can used to select a section of the time series. E.g. --fit-median 674-684.")
    parser.add_argument("--fit-mean", nargs='?', type=str, const="", default=None, help="Fit a Gaussian model to the mean time series.\n"
                        "Optionally, a time range 'start,end' in seconds can used to select a section of the time series. E.g. --fit-mean 674-684.")
    parser.add_argument("FITS FILE", nargs='+', type=str, help="FITS file containing the dynamic spectrum.")

    args = vars(parser.parse_args())
    frequencies = compute_frequency_list_ghz(args["freq"], args["nchans"], args["chan_width"])

    if args["dark"]:
        plt.style.use('dark_background')
    if args["paper"]:
        plt.rcParams.update({'font.size': 13})
    
    if args["fpeaks"]:
        if not os.path.exists("filtered"): os.mkdir("filtered")
        for file in args['FITS FILE']:
            try:
                file_info = extract_filename_info(file)
                if len(file_info) == 5:
                    dm = file_info[2]
                else:
                    dm = file_info[5]
            except:
                print(f"Could not parse DM information from the filename '{file}'. This is necessary for filtering. Exiting..")
                exit(1)
            dyspec = read_fits(file)
            dyspec, dm_time, dm_list, time_series, median_series, peak_idxs = transform_spectrum(
                dyspec, frequencies, args["time_res"], dm, args["chan_avg"], args["time_avg"], args["time_avg_offset"])
            if len(peak_idxs) > 0 and (args["max_cands"] < 0 or len(peak_idxs) <= args["max_cands"]):
                shutil.copy2(file, f"filtered/{file}")
    elif args["snoise"]:
        # Script to study noise levels in dynamic spectra
        median_stdevs = []
        max_stdevs = []
        mean_stdevs = []
        min_stdevs = []
        n_series = []
        max_len = 0
        min_max = []

        snr_scatter = []
        for file in args['FITS FILE']:
            try:
                file_info = extract_filename_info(file)
                if len(file_info) == 5:
                    dm = file_info[2]
                else:
                    snr = file_info[4]
                    dm = file_info[5]
            except:
                print(f"Could not parse DM information from the filename '{file}'. This is necessary for filtering. Exiting..")
                exit(1)
            dyspec = read_fits(file)
            dyspec, dm_time, dm_list, time_series, median_series, peak_idxs = transform_spectrum(
                dyspec, frequencies, args["time_res"], dm, args["chan_avg"], args["time_avg"], args["time_avg_offset"])
            noise_levels = study_time_series_noise(time_series, 0.02)
            if len(peak_idxs):
                brightest_snr = peak_idxs[0][1] #sorted(peak_idxs, key=lambda x : x[1], reverse=True)[0][1]
            else:
                continue
            brightest_peak_idx = sorted(peak_idxs, key=lambda x : x[1], reverse=True)[0][0]
            brightest_val = time_series[brightest_peak_idx]
            #print(brightest_val)
            snr_scatter.append((snr, brightest_snr))
            plt.plot(noise_levels)
            snr_series = [brightest_val / x for x in noise_levels]
            # plt.plot(snr_series)
            # plt.xlabel("Offset")
            # plt.ylabel("Noise")
            # plt.show()
            # continue
            max_len = max(max_len, len(noise_levels))
            n_series.append(noise_levels)
            min_max.append(abs(noise_levels[-1] - noise_levels[1])/noise_levels[-1] * 100 )
            #min_max.append((max(noise_levels) - min(noise_levels)) / max(noise_levels) * 100)
            #stdevs.append(stdev(noise_levels))
        
        # for i in range(2, max_len + 1):
        #     all_i_noises = []
        #     for s in n_series:
        #         all_i_noises.append(stdev(s[:i]))
        #     median_stdevs.append(median(all_i_noises))
        #     mean_stdevs.append(mean(all_i_noises))
        #     max_stdevs.append(max(all_i_noises))
        #     min_stdevs.append(min(all_i_noises))

        #print(snr_scatter)
        #plt.scatter([x[0] for x in snr_scatter], [x[1] for x in snr_scatter])
        #plt.xlabel("2s SNR")
        #plt.ylabel("Full series SNR")
        plt.show()   
        # plt.hist(min_max)
        # plt.ylabel("Count")
        # plt.xlabel("Variation in noise value from 20+s to 2s time series (%)")
        # plt.show()
        # #plt.xlabel("Time (s)")
        # #plt.ylabel("Noise level")
        # #plt.yscale("log")
        # plt.ylabel("Noise deviation")
        # plt.xlabel("Time series length")
        # #plt.yscale('log')
        # plt.plot(median_stdevs)
        # plt.plot(mean_stdevs)
        # plt.plot(min_stdevs)
        # #plt.plot(max_stdevs)
        # plt.legend(["Median", "Mean", "Min", "Max"])
        # plt.show()

    else:
        try:
            extract_filename_info(args['FITS FILE'][0])
            identities = None
            if args["id"]:
                identities = identify_candidates(args["FITS FILE"])
            process_followup_fits_list(args['FITS FILE'], frequencies, args["time_res"], args["chan_avg"],
                                       args["time_avg"],  args["time_avg_offset"], args["dm"], args["interp"], args["disable_norm"],
                                       args["gain"], args["save"], identities, args["fit_mean"], args["fit_median"], args["paper"], args["timeseries_only"])
        except ValueError:
            # Not the standard followp filename.. use standard processing
            main_fig = plt.figure(figsize=(10, 6))
            mean_fig = None
            median_fig = None
            if args["fit_mean"] is not None:
                mean_fig = plt.figure(figsize=(10, 4))
            if args["fit_median"] is not None:
                median_fig = plt.figure(figsize=(10, 4))

            filename = args['FITS FILE'][0]
            if filename.endswith(".fits"):
                dyspec = read_fits()
            else:
                # .fil file format support (experimental)
                header, data = read_filterbank(filename)
                for k, v in header.items():
                    print(k, v)
                #print(header)
                dyspec = data[:, 0, :].T
                print(dyspec.shape)
                dyspec = dyspec[::-1, :]
            plot_spectrum(main_fig, mean_fig, median_fig, dyspec, 0, frequencies, args["time_res"], args["dm"],  args["chan_avg"],
                          args["time_avg"],  args["time_avg_offset"], args["FITS FILE"][0], args["interp"], args["disable_norm"],
                          args["gain"], args["fit_mean"], args["fit_median"], args["paper"], args["timeseries_only"])
            plt.show()


    
    