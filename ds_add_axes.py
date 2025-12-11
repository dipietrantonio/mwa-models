#!/usr/bin/env python3

# Source - https://stackoverflow.com/a
# Posted by ImportanceOfBeingErnest
# Retrieved 2025-11-14, License - CC BY-SA 3.0

import matplotlib.pyplot as plt
import numpy as np
from argparse import ArgumentParser
from PIL import Image

plt.rcParams.update({'font.size': 16})

def add_axis(filename, title, freq_range, time_range):
    img_pil = Image.open(filename)
    img_np = np.asarray(img_pil)

    fig, ax = plt.subplots()
    ax.imshow(img_np, extent=[time_range[0], time_range[1], freq_range[0], freq_range[1]], aspect='auto')
    #ax.set_ylim(freq_range)
    #ax.set_xlim(time_range)
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Frequency (MHz)")
    ax.set_title(title + "\n")
    plt.show()


if __name__ == "__main__":

    parser = ArgumentParser()
    parser.add_argument('--freq', default="138,169", type=str, help="Frequency range in MHz (e.g. 138,169)")
    parser.add_argument('--time', default="0,5", type=str, help="Time range in seconds (e.g. 0,5)")
    parser.add_argument('--title', default="Dynamic spectrum", type=str, help="Plot title.")
    parser.add_argument('IMAGE', nargs=1)
    args = vars(parser.parse_args())

    to_float_list = lambda x : [float(a) for a in x.split(',')]
    freq_range = to_float_list(args['freq'])
    time_range = to_float_list(args['time'])
    
    add_axis(args['IMAGE'][0], args['title'], freq_range, time_range)