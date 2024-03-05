import matplotlib.pyplot as plt
from math import sqrt
from models import MWA_PHASE_1


#########################################################################
#                      SENSITIVITY STUDY
#########################################################################

def sensitivity(t_int, delta_v):
    SEFD = 28001.7919 #16242.1012 (216MHz)
    N_a = MWA_PHASE_1.n_antennas
    eff = 1
    B = (N_a / 2) * (N_a - 1)
    return SEFD /(eff * sqrt(delta_v * B * t_int))


def frb_min_fluence(N, int_time, delta_v):
    return N * sensitivity(int_time, delta_v) * int_time 


def sensitivity_study():
    N = 10
    T =  [i for i in range(1, 100)]
    V = [sensitivity(t_int, 1.28e6) for t_int in T]
    F = [frb_min_fluence(N, t_int, 1.28e6) for t_int in T]

    fig, (ax1, ax2) = plt.subplots(1, 2)

    ax1.plot(T, V)
    ax1.set_title("Noise / Sensitivity")
    ax1.set_xlabel("Integration time (ms)")
    ax1.set_ylabel("Flux (Jy)")
    ax2.plot(T, F)
    ax2.set_xlabel("Integration time (ms)")
    ax2.set_ylabel("Fluence (Jy ms)")
    ax2.set_title("Minimum detectable fluence")
    
    plt.show()


if __name__ == "__main__":
    sensitivity_study()