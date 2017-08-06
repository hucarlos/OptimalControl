import numpy as np
import matplotlib.pyplot as plt
import itertools

if __name__ == '__main__':

    data = np.loadtxt("Points.txt")
    x = data[0, :]
    y = data[1, :]

    # multiple box plots on one figure
    plt.plot(x, y, 'ro')

    plt.show()

    a = list(itertools.permutations([1, -1, 1, -1, 1, -1]))
    print len(a)