import seaborn as so
import matplotlib.pyplot as plt
from math import log, pi
from models import MWA_PHASE_1
############################################################################
#                      COMPUTATIONAL COST STUDY
############################################################################

TIME_RES = MWA_PHASE_1.time_res

def npixels(Bmax):
    return (610/((1.9/Bmax) * 180 / pi )**2) # convert to degree

def beamforming_cost(Np, total_time):
    return Np * MWA_PHASE_1.n_antennas * (total_time / TIME_RES)

def correlation_cost(n_timesteps):
    N_a = MWA_PHASE_1.n_antennas
    B = (N_a / 2) * (N_a - 1)
    return B * n_timesteps

def fft_cost(Np):
    return Np * log(Np, 2)

def gridding_cost():
    N_a = MWA_PHASE_1.n_antennas
    B = (N_a / 2) * (N_a - 1)
    return B

def imaging_cost(Np, int_time, total_time):
    intervals = total_time / int_time
    timesteps = int_time / TIME_RES
    return (correlation_cost(timesteps) + gridding_cost() + fft_cost(Np)) * intervals

def beamforming_vs_imaging():
    INT_TIME = 0.05
    MAX_TIME = 60 * 60
    BLINES = [3000, 6000, 10000]
    for max_baseline_length in BLINES:
        Np = npixels(max_baseline_length)

        bc = []
        ic = []
        TIME = [INT_TIME * i for i in range(1, int(MAX_TIME / INT_TIME))]
        for t in TIME:
            bc.append(beamforming_cost(Np, t))
            ic.append(imaging_cost(Np, INT_TIME, t))
        
        diff = [ y / x * 100 for x, y in zip(bc, ic)]
        print(bc[0], ic[0])
        plt.plot(TIME, diff)
        
        # plt.plot(TIME, bc)
        # plt.plot(TIME, ic)
        # plt.legend(["Beamforming", "Imaging"])
        # plt.show()
        plt.ylabel("Percentage %")
        plt.xlabel("Observation time")
        plt.legend(["Cost of imaging w.r.t. beamforming"])
    
    plt.legend(["Bmax = {}k".format(bmax) for bmax in BLINES])
    plt.show()



if __name__ == "__main__":
    beamforming_vs_imaging()