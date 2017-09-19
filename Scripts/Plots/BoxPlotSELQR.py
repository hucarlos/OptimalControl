import numpy as np
import matplotlib.pyplot as plt
import matplotlib


if __name__ == '__main__':
    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['ps.fonttype'] = 42

    data = np.loadtxt('Results4.txt')

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

    bp = ax1.boxplot(cost, showmeans=True, showfliers=False)
    ax1.set_xticklabels(['E-LQR', 'RE-LQR'])
    ax1.set_title("Cost")
    ax1.grid()

    ax2.boxplot(times, showmeans=True, showfliers=False)
    ax2.set_xticklabels(['E-LQR', 'RE-LQR'])
    ax2.set_title("Time (ms)")
    ax2.grid()


    print 'Mean iters: ', np.mean(data[:, 4]), np.mean(data[:, 5])
    print 'Mean costs: ', np.mean(costSELQR), np.mean(costRSELQR)
    print 'Mean time: ', np.mean(data[:, 0]), np.mean(data[:, 1])

    plt.show()
