SPEED_OF_LIGHT = 299792458 # m/s
K = 4.15

def dispersive_delay(DM, f1_ghz, f2_ghz):
    return K * DM * (f1_ghz**(-2) - f2_ghz**(-2))


class Interferometer:
    def __init__(self, n_antennas, n_pols, n_channels, bit_per_sample, FoV, longest_baseline, time_res):
        self.n_antennas = n_antennas
        self.n_pols = n_pols
        self.n_channels = n_channels
        self.sampling_freq = 1 / time_res
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

    def __init__(self, image_side, bits_per_pixel, correlator : Correlator):
        self.correlator = correlator
        self.image_side = image_side
        self.bits_per_pixel = bits_per_pixel
    
    @property
    def data_rate(self):
        return self.image_side ** 2 * self.bits_per_pixel \
            * self.correlator.n_channels * self.correlator.n_intervals
    

class Dedispersion:

    def __init__(self, n_DMs, imager : Imager):
        self.n_DMs = n_DMs
        self.imager = imager
    
    @property
    def data_date(self):
        return self.n_DMs * self.imager.correlator.n_intervals * self.imager.image_side ** 2 * self.imager.bits_per_pixel


MWA_PHASE_1 = Interferometer(128, 2, 3072, 8, 610, 2864, 1e-4)

EDA2 = Interferometer(256, 2, None, None, 12000, 35, 1e-6)