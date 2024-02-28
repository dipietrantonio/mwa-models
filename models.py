SPEED_OF_LIGHT = 299792458 # m/s

class Interferometer:
    def __init__(self, n_antennas, n_pols, n_channels, sampling_freq, bit_per_sample, FoV, longest_baseline, time_res):
        self.n_antennas = n_antennas
        self.n_pols = n_pols
        self.n_channels = n_channels
        self.sampling_freq = sampling_freq
        self.bit_per_sample = bit_per_sample
        self.FoV = FoV
        self.longest_baseline = longest_baseline
        self.time_res = time_res
    
    @property
    def data_rate(self):
        return self.n_antennas * self.n_pols * self.n_channels * self.sampling_freq * self.bit_per_sample
    

class Correlator:

    def __init__(self, interferometer : Interferometer, integration_time, channels_to_avg):
        self.interferometer = interferometer
        self.integration_time = integration_time
        self.channels_to_avg = channels_to_avg
        self.n_channels = interferometer.n_channels / channels_to_avg
        self.n_baselines = (interferometer.n_antennas / 2) * (interferometer.n_antennas + 1)
        self.n_intervals = (1 / self.integration_time)
        self.bits_per_sample = 2 * 32 # 1 complex number, float32 for real and imag
    
    @property
    def data_rate(self):
        return self.n_intervals * self.n_channels * \
            self.n_baselines * self.interferometer.n_pols ** 2 * self.bits_per_sample


class Imager:

    def __init__(self, correlator : Correlator):
        self.correlator = correlator
        self.image_side = 1024
        self.bits_per_pixel = 32
    
    @property
    def data_rate(self):
        return self.image_side ** 2 * self.bits_per_pixel \
            * self.correlator.n_channels * self.correlator.n_intervals

MWA_PHASE_1 = Interferometer(128, 2, 3092, 1e4, 8, 610, 2864, 1e-4)