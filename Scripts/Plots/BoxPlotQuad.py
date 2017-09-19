import numpy as np
import matplotlib.pyplot as plt
import matplotlib

if __name__ == '__main__':

    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['ps.fonttype'] = 42

    data = np.loadtxt('QR8.txt')

    # print data.shape
    # data = [data[m] for d, m in zip(data.T, mask.T)]
    # mask = ~np.isinf(data)
    # Filter data using np.isnan
    times = [data[:, 1], data[:, 2]]

    mask = ~np.isinf(data[:, 3])

    costSELQR = data[:, 3]
    costSELQR = costSELQR[mask]

    costRSELQR = data[:, 4]
    costRSELQR = costRSELQR[mask]

    timeSELQR = data[:, 1]
    timeSELQR = timeSELQR[mask]

    timeRSELQR = data[:, 2]
    timeRSELQR = timeRSELQR[mask]

    itersSELQR = data[:, 5]
    itersSELQR = itersSELQR[mask]

    itersRSELQR = data[:, 6]
    itersRSELQR = itersRSELQR[mask]

    cost  = [costSELQR, costRSELQR]
    times = [timeSELQR, timeRSELQR]


    # multiple box plots on one figure

    f, (ax1, ax2) = plt.subplots(1, 2, sharex=False, sharey=False)

    bp = ax1.boxplot(cost, showfliers=False, showmeans=True)
    ax1.set_xticklabels(['E-LQR', 'RE-LQR'])
    ax1.set_title("Cost")
    ax1.grid()

    ax2.boxplot(times, showfliers=False, showmeans=True)
    ax2.set_xticklabels(['E-LQR', 'RE-LQR'])
    ax2.set_title("Time (s)")
    ax2.grid()


    print 'Mean iters: ', np.mean(itersSELQR), np.mean(itersRSELQR)
    print 'Mean costs: ', np.mean(costSELQR), np.mean(costRSELQR)
    print 'Mean time: ', np.mean(timeSELQR), np.mean(timeRSELQR)

    plt.show()
