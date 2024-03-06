import numpy as np
import matplotlib.pyplot as plt
from models import dispersive_delay



def generate_sweep(DM, f_low_mhz, f_high_mhz, freq_res_mhz, int_time_s, signal):
    n_channels = int((f_high_mhz - f_low_mhz) / freq_res_mhz)
    freqs_ghz = np.linspace(f_low_mhz / 1000, f_high_mhz / 1000, n_channels)
    delays_s = [dispersive_delay(DM, freqs_ghz[i], freqs_ghz[i + 1]) / 1000
                for i in range(len(freqs_ghz) - 1)]
    total_delay_s = sum(delays_s)
    centre_freq_mhz = (f_high_mhz + f_low_mhz ) / 2
    background = np.random.normal(0, 0.5, (n_channels, int(total_delay_s / int_time_s) + 1))
    current_time_s = 0
    for i in range(n_channels - 1):
        D = delays_s[-i]
        current_time_s += D
        background[i, int(current_time_s/int_time_s)] = signal()
    return background


background = generate_sweep(600, 135, 165, 0.01, 0.05, lambda : 4)


plt.imshow(background)
plt.show()