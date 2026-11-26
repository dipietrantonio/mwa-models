from math import log
import matplotlib.pyplot as plt
import numpy as np

def nondetection_rate(p, t):
    return - log(p) / t


if __name__ == "__main__":

    print(nondetection_rate(0.05, 1))
    probs = np.arange(0.1, 1, 0.01)
    ndr = np.array([nondetection_rate(p, 1) for p in probs])

    plt.plot(probs, ndr)
    plt.show()