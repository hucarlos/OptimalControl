import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':

    data = np.loadtxt("SELQR1503903180.txt")
    print np.mean(data, axis=0)
    print np.var(data, axis=0)