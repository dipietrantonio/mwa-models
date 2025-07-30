import sys
from astropy.io import fits
import numpy
from statistics import mean, stdev, median, mode
import matplotlib.pyplot as plt


def read_image(filename):
    my_fits = fits.open(filename)
    return my_fits[0].data


def process_image_list(image_list):
    all_stats = []
    for i, image_file in enumerate(image_list):
        print(f"processing file {i}/{len(image_list)}")
        img = read_image(image_file)
        vals = numpy.reshape(img, newshape=(img.shape[0] * img.shape[1], 1))
        plt.hist(vals, 50)
        plt.show()
        
        stats = (image_file, numpy.mean(vals), numpy.std(vals), numpy.median(vals), numpy.min(vals), numpy.max(vals))
        all_stats.append(stats)
    return all_stats


def print_stats_table(stats):
    header = ["filename", "mean", "stdev", "median", "min", "max"]
    print("{:35s} {:10s} {:10s} {:10s}".format(*header))
    for s in stats:
        print("{:35s} {:10f} {:10f} {:10f} {:10f} {:10f}".format(*s))


if __name__ == "__main__":
    all_stats = process_image_list(sys.argv[1:])
    print_stats_table(all_stats)
