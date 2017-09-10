import numpy as np
import matplotlib.pyplot as plt


def reject_outliers(data, m=2):
    return data[abs(data - np.mean(data)) < m * np.std(data)]


if __name__ == '__main__':
    data = np.loadtxt('WinsrSELQRChol.txt')

    # print data.shape
    # data = [data[m] for d, m in zip(data.T, mask.T)]
    # mask = ~np.isinf(data)
    # Filter data using np.isnan
    times = [data[:, 0], data[:, 1]]

    mask = ~np.isinf(data[:, 2])

    costSELQR = data[:, 2]
    costSELQR = costSELQR[mask]

    costRSELQR = data[:, 3]
    costRSELQR = costRSELQR[mask]

    cost = [costSELQR, costRSELQR]

    iters = [data[:, 4], data[:, 5]]

    # multiple box plots on one figure

    f, (ax1, ax2) = plt.subplots(1, 2, sharex=False, sharey=False)

    bp = ax1.boxplot(cost, 0, '')
    ax1.set_xticklabels(['SELQR', 'QRSELQR'])
    ax1.set_ylabel("Cost")

    ax2.boxplot(times, 0, '')
    ax2.set_xticklabels(['SELQR', 'QRSELQR'])
    ax2.set_ylabel("Time (ms)")

    print 'Mean iters: ', np.mean(data[:, 4]), np.mean(data[:, 5])
    print 'Mean costs: ', np.mean(costSELQR), np.mean(costRSELQR)
    print 'Mean time: ', np.mean(data[:, 0]), np.mean(data[:, 1])

    print 'Std iters: ', np.std(data[:, 4]), np.std(data[:, 5])
    print 'Std costs: ', np.std(costSELQR), np.std(costRSELQR)
    print 'Std time: ', np.std(data[:, 0]), np.std(data[:, 1])

    plt.show()
