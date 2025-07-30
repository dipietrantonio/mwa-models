import numpy as np
import matplotlib.pyplot as plt
from models import dispersive_delay



def generate_sweep(DM, f_low_mhz, f_high_mhz, freq_res_mhz, int_time_s, start_time, signal):
    n_channels = int((f_high_mhz - f_low_mhz) / freq_res_mhz)
    freqs_ghz = np.linspace(f_low_mhz / 1000, f_high_mhz / 1000, n_channels)
    delays_s = [dispersive_delay(DM, freqs_ghz[i], freqs_ghz[i + 1]) / 1000
                for i in range(len(freqs_ghz) - 1)]
    total_delay_s = sum(delays_s)
    n_timesteps = int(total_delay_s / int_time_s) * 10
    background = np.random.normal(0.5, 0.05, (n_channels, n_timesteps))
    current_time_s = start_time
    for i in range(n_channels - 1):
        D = delays_s[-i]
        current_time_s += D
        idx = int(current_time_s/int_time_s)
        if idx < n_timesteps:
            background[i, idx] = signal()
        else:
            break
    return background


def dedisperse_with_dm(background, freqs_ghz, time_res_s, DM):
    delays_s = [dispersive_delay(DM, freqs_ghz[i], freqs_ghz[i + 1]) / 1000
            for i in range(len(freqs_ghz) - 1)]
    time_series = []
    C, T = background.shape
    for t in range(T):
        sum = 0
        current_time_s = time_res_s * t
        for i in range(C):
            current_time_s += delays_s[-i]
            idx = int(current_time_s / time_res_s)
            if idx >= T:
                break
            sum += background[i, idx]
        time_series.append(sum)
    return time_series


def peak_detection(time_series, SN = 2):
    if len(time_series) < 1:
        return None
    peaks = []
    means = []
    running_mean = time_series[0]
    w = 0.7
    for i, v in enumerate(time_series[1:], 1):
        if v >= SN * running_mean:
            peaks.append(i)
        means.append(running_mean)
        running_mean = w * running_mean + (1 - w) * v
    return means, peaks


print(dispersive_delay(600, 0.138, 0.169)/1000)
exit(0)


background = generate_sweep(600, 135, 165, 0.01, 0.05, 12, lambda : 24)

plt.imshow(background)
plt.show()

f_high_mhz = 165
f_low_mhz = 135
freq_res_mhz = 0.01

n_channels = int((f_high_mhz - f_low_mhz) / freq_res_mhz)
freqs_ghz = np.linspace(f_low_mhz / 1000, f_high_mhz / 1000, n_channels)

ts = dedisperse_with_dm(background, freqs_ghz, 0.05, 600)

tp = [i * 0.05 for i  in range(len(ts))]


plt.plot(tp, ts)

means, peaks = peak_detection(ts)
plt.plot(tp[:-1], means )
peak_times = [ p * 0.05 for p in peaks]
print(peak_times)
plt.show()