import seaborn as so
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


def sensitivity_study():
    T =  [i * 0.01 for i in range(1, 100)]
    V = [sensitivity(t_int, 1.28e6) for t_int in T]

    plt.plot(T, V)
    plt.xlabel("Time (s)")
    plt.ylabel("Sensitivity (Jy)")
    plt.show()


if __name__ == "__main__":
    sensitivity_study()