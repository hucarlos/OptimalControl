import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':

    data = np.loadtxt("Wins.txt")
    times = [data[:,0], data[:, 1]]

    iters = [data[:,4], data[:, 5]]

    # multiple box plots on one figure

    f, (ax1, ax2) = plt.subplots(1, 2, sharex=False, sharey=False)

    bp = ax1.boxplot(iters)
    ax1.set_xticklabels(['SELQR', 'QRSELQR'])
    ax1.set_ylabel("Iterations")

    ax2.boxplot(times, 0, '')
    ax2.set_xticklabels(['SELQR', 'QRSELQR'])
    ax2.set_ylabel("Time (ms)")

    plt.show()