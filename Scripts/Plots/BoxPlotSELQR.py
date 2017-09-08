import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':

    data = np.loadtxt("WinsrSELQRChol.txt")

    # Filter data using np.isnan
    # mask = ~np.isinf(data)
    # data = [data[m] for d, m in zip(data.T, mask.T)]
    # print data.shape

    times = [data[:,0], data[:, 1]]

    costSELQR = data[:,2]
    costSELQR = costSELQR[~np.isnan(costSELQR)]

    costRSELQR = data[:, 3]
    costRSELQR = costRSELQR[~np.isnan(costSELQR)]

    cost  = [costSELQR, costRSELQR]

    iters = [data[:,4], data[:, 5]]

    # multiple box plots on one figure

    f, (ax1, ax2) = plt.subplots(1, 2, sharex=False, sharey=False)

    bp = ax1.boxplot(cost, 0, '')
    ax1.set_xticklabels(['SELQR', 'QRSELQR'])
    ax1.set_ylabel("Cost")

    ax2.boxplot(times, 0, '')
    ax2.set_xticklabels(['SELQR', 'QRSELQR'])
    ax2.set_ylabel("Time (ms)")

    print np.mean(data, axis=0)

    plt.show()