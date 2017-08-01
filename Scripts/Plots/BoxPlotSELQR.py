import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':

    data = np.loadtxt("Wins.txt")
    times = [data[:,0], data[:, 1]]
    iters = [data[:,4], data[:, 5]]

    # multiple box plots on one figure
    plt.figure()
    plt.boxplot(iters, 0, '')

    plt.show()