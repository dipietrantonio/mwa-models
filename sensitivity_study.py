import matplotlib.pyplot as plt
from basics import sensitivity_jy as sensitivity, frb_min_fluence_jyms as frb_min_fluence
from argparse import ArgumentParser
#########################################################################
#                      SENSITIVITY STUDY
#########################################################################



def sensitivity_study(freq, delta_t):
    SNR = 10
    T =  [i * delta_t for i in range(1, 50)]
    V = [sensitivity(freq, t_int, 128, 30.72e6) for t_int in T]
    F = [frb_min_fluence(SNR, freq, t_int, 128, 30.72e6) for t_int in T]

    fig, (ax1, ax2) = plt.subplots(1, 2)

    ax1.plot(T, V)
    ax1.set_title("Noise / Sensitivity")
    ax1.set_xlabel("Integration time (s)")
    ax1.set_ylabel("Flux (Jy)")
    ax2.plot(T, F)
    ax2.set_xlabel("Integration time (s)")
    ax2.set_ylabel("Fluence (Jy ms)")
    ax2.set_title("Minimum detectable fluence")
    
    plt.show()


if __name__ == "__main__":

    parser = ArgumentParser()
    parser.add_argument("-b", "--bandwidth", type=float, default=30.72, help="Bandwidth (in MHz)")
    parser.add_argument("-a", "--antennas", type=int, default=128, help="Number of antennas")
    parser.add_argument("-t", "--inttime", default=0.01, type=float, help="Integration time (in seconds)")
    parser.add_argument("-f", "--frequency", default=150.0, type=float, help="Central frequency (in MHz)")
    parser.add_argument("--snr", type=float, default=10.0, help="Minimum Signal-to-Noise Ratio (SNR)")

    args = vars(parser.parse_args())
    for t in [0.010, 0.020, 0.050]:
        S = sensitivity(args['frequency'] * 1e6, t, args['antennas'], args['bandwidth'] * 1e6)
        F = frb_min_fluence(args['snr'], args['frequency'] * 1e6, t, args['antennas'], args['bandwidth'] * 1e6)
        print(f"Sensitivity is {S:.2f}, minimum detectable fluence is {F:.2f}")
    sensitivity_study(args['frequency'] * 1e6, args['inttime'])