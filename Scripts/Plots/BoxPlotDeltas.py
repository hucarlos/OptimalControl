import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':

    data1 = np.loadtxt('Results2.txt')
    data2 = np.loadtxt('Results3.txt')


    # multiple box plots on one figure

    f, (ax1, ax2) = plt.subplots(1, 2, sharex=False, sharey=False)

    cost = [data1[:,3], data2[:,3]]
    bp = ax1.boxplot(cost, 0, '')
    ax1.set_xticklabels(['SELQR', 'QRSELQR'])
    ax1.set_ylabel("Cost")

    times = [data1[:,1], data2[:,1]]
    ax2.boxplot(times, 0, '')
    ax2.set_xticklabels(['SELQR', 'QRSELQR'])
    ax2.set_ylabel("Time (ms)")

    plt.show()
