import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':

    data = np.loadtxt("SELQR1503423184.txt")
    print np.mean(data, axis=0)