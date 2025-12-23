#!/usr/bin/env python3
from basics import dispersive_delay_s
from argparse import ArgumentParser


if __name__ == "__main__":

    parser = ArgumentParser()
    parser.add_argument("--dm", type=float, required=True, help="Dispersion measure")
    parser.add_argument("--freq", type=float, default=154.237, help="Central frequency (in MHz) of the central frequency channel.")
    parser.add_argument("--bandwidth", type=float, default=30.72, help="Width of the band (in MHz) the signal propagates through.")

    args = vars(parser.parse_args())
    low_freq_mhz = args["freq"] - args["bandwidth"] / 2
    high_freq_mhz = args["freq"] + args["bandwidth"] / 2
    
    delay = dispersive_delay_s(args["dm"], low_freq_mhz / 1e3, high_freq_mhz /1e3)
    print(delay)
