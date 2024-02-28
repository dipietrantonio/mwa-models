import seaborn as so
import matplotlib.pyplot as plt
from math import log, pi
from models import MWA_PHASE_1 as MWA_MODEL, SPEED_OF_LIGHT


############################################################################
#                      COMPUTATIONAL COST STUDY
############################################################################

TIME_RES = MWA_MODEL.time_res # 1e-4 s
FOV = 610  # Deg^2
FREQ = 150 * 1e6

def freq_to_wavelength(freq):
    return SPEED_OF_LIGHT / FREQ


def npixels(Bmax):
    lmbda = freq_to_wavelength(FREQ)
    resolution = (lmbda/Bmax) * 180 / pi
    return (FOV / resolution**2)


def beamforming_cost_per_time_sample(Np):
    return Np * MWA_MODEL.n_antennas


def correlation_cost_per_time_sample():
    N_a = MWA_MODEL.n_antennas
    B = (N_a / 2) * (N_a - 1)
    return B


def fft_cost_per_time_sample(Np, integration_time):
    n_samples = integration_time / TIME_RES
    return Np * log(Np, 2) / n_samples


def gridding_cost_per_time_sample(integration_time):
    n_samples = integration_time / TIME_RES
    N_a = MWA_MODEL.n_antennas
    B = (N_a / 2) * (N_a - 1)
    return B / n_samples


def imaging_cost_per_time_sample(Np, integration_time):
    return correlation_cost_per_time_sample() + \
            gridding_cost_per_time_sample(integration_time) + \
            fft_cost_per_time_sample(Np, integration_time)


def plot_imaging_costs_frb_case():
    integration_time = 0.05
    Np = npixels(MWA_MODEL.longest_baseline)
    ccost = correlation_cost_per_time_sample()
    gcost = gridding_cost_per_time_sample(integration_time)
    fcost = fft_cost_per_time_sample(Np, integration_time)
    costs = [ccost, gcost, fcost]
    m = sum(costs)
    costs_perc = [x / m * 100 for x in costs]
    # print(costs)
    # print(costs_perc)
    # plt.barplot(total_costs)
    plt.bar(["Correlation", "Gridding", "FFT"], costs_perc)
    plt.xlabel("Phases")
    plt.ylabel("Fraction of total cost (%)")
    plt.title("Computational cost (per time sample) distribution over imaging steps (MWA Phase 1, freq 150MHz, int. time = 50ms)")
    plt.show()



def plot_imaging_costs_as_function_of_int_time():
    STEP_TIME = 0.005 # s
    MAX_TIME = 1 # s
    INT_TIMES = [i * 0.005 for i in range(1, int(MAX_TIME/STEP_TIME))]

    ccosts = []
    gcosts = []
    fcosts = []
    Np = npixels(MWA_MODEL.longest_baseline)
    
    for integration_time in INT_TIMES:
        ccost = correlation_cost_per_time_sample()
        gcost = gridding_cost_per_time_sample(integration_time)
        fcost = fft_cost_per_time_sample(Np, integration_time)
        costs = [ccost, gcost, fcost]
        m = sum(costs)
        costs_perc = [x / m * 100 for x in costs]

        ccosts.append(costs_perc[0])
        gcosts.append(costs_perc[1])
        fcosts.append(costs_perc[2])
    
    plt.plot(INT_TIMES, ccosts)
    plt.plot(INT_TIMES, gcosts)
    plt.plot(INT_TIMES, fcosts)
    
    plt.xlabel("Integration time (s)")
    plt.ylabel("Fraction of total cost (%)")
    plt.title("Computational cost per time sample of every imaging step as a function of integration time (MWA Phase I, freq = 150MHz).")
    plt.legend(["Correlation", "Gridding", "FFT"])
    plt.show()


def plot_beamforming_vs_imaging_frb_case():
    INTEGRATION_TIME = 0.05
    Np = npixels(MWA_MODEL.longest_baseline)
    bf_cost = beamforming_cost_per_time_sample(Np)
    img_cost = imaging_cost_per_time_sample(Np, INTEGRATION_TIME)
    print(f"Beamforming cost = {bf_cost}\nImaging cost = {img_cost}")
    plt.bar(["Beamforming", "Imaging"], [bf_cost, img_cost], width=0.2)
    plt.xlabel("Strategy")
    plt.yscale("log")
    plt.ylabel("Number of computational steps")
    plt.title("Computational cost of beamforming and imaging for MWA Phase I (freq = 150MHz, int_time = 50ms)")
    plt.show()


def plot_beamforming_vs_imaging_as_number_of_pixes():
    INTEGRATION_TIME = 0.05
    relative_cost = []
    X = list(range(1, int(1e2)))
    for Np in X:
        bf_cost = beamforming_cost_per_time_sample(Np)
        img_cost = imaging_cost_per_time_sample(Np, INTEGRATION_TIME)
        relative_cost.append(img_cost / bf_cost * 100)
    
    plt.plot(X, relative_cost)
    plt.plot(X, len(X) * [100])
    plt.xlabel("Number of pixels")
    plt.ylabel("Percentage (%)")
    plt.title("Computational cost of imaging w.r.t. beamforming as a function of pixels (MWA Phase I, int_time = 50ms)")
    plt.legend(["Relative cost", "100% mark"])
    plt.show()



if __name__ == "__main__":
    # plot_imaging_costs_frb_case()
    #plot_imaging_costs_as_function_of_int_time()
    plot_beamforming_vs_imaging_frb_case()
    # plot_beamforming_vs_imaging_as_number_of_pixes()