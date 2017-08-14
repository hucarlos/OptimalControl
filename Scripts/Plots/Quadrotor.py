import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':

    dataSELQR  = np.loadtxt("SELQR.txt")
    dataRSELQR = np.loadtxt("iQRSELQR.txt")

    # multiple box plots for the controls
    f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex=True, sharey=False)

    rows = 150-1
    ax1.plot(dataSELQR[0:rows,12], label='SELQR')
    ax1.plot(dataRSELQR[0:rows, 12], label = 'RSELQR')
    ax1.set_xlabel('$k$')
    ax1.set_ylabel('$f_1$')
    ax1.legend()
    ax1.grid()

    ax2.plot(dataSELQR[0:rows, 13], label='SELQR')
    ax2.plot(dataRSELQR[0:rows, 13], label='RSELQR')
    ax2.set_xlabel('$k$')
    ax2.set_ylabel('$f_2$')
    ax2.grid()

    ax3.plot(dataSELQR[0:rows, 14], label='SELQR')
    ax3.plot(dataRSELQR[0:rows, 14], label='RSELQR')
    ax3.set_xlabel('$k$')
    ax3.set_ylabel('$f_3$')
    ax3.grid()

    ax4.plot(dataSELQR[0:rows, 15], label='SELQR')
    ax4.plot(dataRSELQR[0:rows, 15], label='RSELQR')
    ax4.set_xlabel('$k$')
    ax4.set_ylabel('$f_4$')
    ax4.grid()

    plt.show()