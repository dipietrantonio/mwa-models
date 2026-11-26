#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt


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

if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
            print("Usage: visualise_bf_ts.py <numpy dump>")
            exit(0)
    plot_separately = True
    final_series = None
    max_peak = 0
    peak_file = None
    for dumpfile in sys.argv[1:]:
        data = np.load(dumpfile) / 2 # because of XX and YY convention
        avg_factor = 2
        new_data = np.zeros((int(data.shape[0] / avg_factor),))
        if final_series is None:
            final_series = np.zeros((int(data.shape[0] / avg_factor),))
        for i in range(new_data.shape[0]):
            new_data[i] = np.mean(data[i*avg_factor:(i+1)*avg_factor])
        data = new_data
        if plot_separately:        
            peaks = peak_finding(data, 4)
            for p in peaks:
                if p[1] > max_peak:
                    print(p[1])
                    max_peak = p[1]
                    peak_file = dumpfile
            if len(peaks) >0: plt.plot([i * 0.02 * avg_factor for i in range(len(data))], data)
        else:
            final_series += data
    
    if not plot_separately:
        final_series /= (len(sys.argv) - 1)
        print(np.std(final_series))
        plt.plot([i * 0.02 * avg_factor for i in range(len(final_series))], final_series)

    plt.ylabel("Intensity (arb units)")
    plt.xlabel("Time (s)")
    plt.show()

    print(peak_file)