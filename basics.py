from math import sqrt, pi, pow

SPEED_OF_LIGHT = 299792458 # m/s

def dispersive_delay_ms(DM, f1_ghz, f2_ghz):
    return 4.15 * DM * (f1_ghz**(-2) - f2_ghz**(-2))


def freq_to_wavelength_m(freq_hz):
    return SPEED_OF_LIGHT / freq_hz

#
# SENSITIVITY STUDY
#

def SEFD_jy(freq_hz):
    """
    System Equivalent Flux Density

    Note this implementation is not entirely correct, as SEFD is dependent on direction and polarization.

    See MWA Full Embedded Element Beam Model (Marcin Sokolowski 2017)
    """
    lamda = freq_to_wavelength_m(freq_hz)
    A_eff = 4.75 #m^2
    T_sky = 60*lamda**2.25 # K
    T_rc = 180 # K
    k_b = 1.380649e-23
    T = T_sky + T_rc
    val = 2 * k_b * T / A_eff # in W / m2 / Hz
    val_in_jy = val / 1e-26
    return val_in_jy

def sensitivity_jy(frequency_hz, integration_time_s, n_antennas, bandwidth_hz):
    # SEFD = 28001.7919 #16242.1012 (216MHz)
    eff = 1
    B = (n_antennas / 2) * (n_antennas - 1) # Should we include the number of polarizations
    return SEFD_jy(frequency_hz) /(eff * sqrt(bandwidth_hz * B * integration_time_s))


def frb_min_fluence_jyms(SNR, frequency_hz, integration_time_s, n_antennas, bandwidth_hz):
    return SNR * sensitivity_jy(frequency_hz, integration_time_s, n_antennas, bandwidth_hz) * integration_time_s * 1000


def pixsize_deg(longest_baseline_m, freq_hz):
    lmbda = freq_to_wavelength_m(freq_hz)
    u_max = longest_baseline_m / lmbda
    pixsize_rad = 1 / (2 * u_max)
    return  pixsize_rad * (180 / pi)
    

def npixels(FoV_deg2, longest_baseline_m):
    resolution = pixsize_deg(longest_baseline_m)
    return (FoV_deg2 / resolution**2)


def expected_frb_rate(ref_rate, ref_freq_hz, ref_fluence, freq_hz, fluence, alpha):
    return ref_rate * pow(freq_hz /ref_freq_hz, alpha) * pow(fluence / ref_fluence, -3.0/2.0)
