from basics import dispersive_delay_s
from matplotlib import pyplot as plt
plt.rcParams.update({'font.size': 20})

def compute_required_bandwidth_khz(centre_fq_mhz, dm, tres):
    """
    Given a time resolution in seconds, compute the fine channel bandwidth
    whose associated dispersive delay matches the specified time resolution.
    """
    # at the moment do not do anything for scattering
    # TODO: compute for the highest frequency channel in the band?
    cf_ghz = centre_fq_mhz / 1000
    bwghz = 0.03072
    # One can either compute it from the centre frequency of the band,
    # or from the bottom one. Let's go with the centre frequency.
    # start_freq_ghz = cf_ghz - bwghz / 2
    start_freq_ghz = cf_ghz

    delta_ghz = 0.000001 # 10 kHz
    f_hi_ghz = start_freq_ghz + delta_ghz
    
    delay = dispersive_delay_s(dm, start_freq_ghz, f_hi_ghz)
    while delay < tres:
        f_hi_ghz += delta_ghz
        delay = dispersive_delay_s(dm, start_freq_ghz, f_hi_ghz)
    
    return (f_hi_ghz - start_freq_ghz) * 1e6



def plot_freq_res(delta_t = 0.001, DM=600, central_freq_mhz = 150):
    T =  [i * delta_t for i in range(1, 50)]
    BW = [compute_required_bandwidth_khz(central_freq_mhz, DM, t) for t in T]
    plt.plot(T, BW)
    plt.title("Required bandwidth resolution as a function of integration time.")
    plt.xlabel("Integration time (s)")
    plt.ylabel("Bandwidth (kHz)")
    plt.show()


def plot_freq_as_dm(central_freq_mhz = 150):
    for t in [0.01, 0.02, 0.05]:
        DM =  [600 + 10 * i for i in range(0, 40)]
        BW = [compute_required_bandwidth_khz(central_freq_mhz, d, t) for d in DM]
        plt.plot(DM, BW)

    plt.legend(["10 ms", "20 ms", "50 ms"])
    plt.title("Required bandwidth resolution as a function of DM.")
    plt.xlabel("Dispersion Measure (pc cm-3)")
    plt.ylabel("Bandwidth (kHz)")
    plt.show()

if __name__ == "__main__":
    plot_freq_res()
    plot_freq_as_dm()
