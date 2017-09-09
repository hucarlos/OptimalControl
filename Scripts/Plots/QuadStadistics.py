import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':

    data = np.loadtxt("QR_ELLIP4.txt")

    times = [data[:,1], data[:, 2]]

    iters = [data[:,5], data[:, 6]]

    cost = [data[:,3], data[:, 4]]

    # multiple box plots on one figure

    f, (ax1, ax2, ax3) = plt.subplots(1, 3, sharex=False, sharey=False)

    bp = ax1.boxplot(iters)
    ax1.set_xticklabels(['SELQR', 'QRSELQR'])
    ax1.set_ylabel("Iterations")

    ax2.boxplot(times, 0, '')
    ax2.set_xticklabels(['SELQR', 'QRSELQR'])
    ax2.set_ylabel("Time (s)")

    ax3.boxplot(cost, 0, '')
    ax3.set_xticklabels(['SELQR', 'QRSELQR'])
    ax3.set_ylabel("Accumulated cost")

    print np.mean(data[:, 1:7], axis=0)

    plt.show()