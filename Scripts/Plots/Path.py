import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":


    path         = "/Volumes/Aslan/OptimalControl/Configurations/Map1/"
    filename1    = "SELQR1.txt"
    filename2    = "iQRLQR1.txt"
    filename3    = "iQRSELQR1.txt"

    file_name_selqr     = path + filename1
    file_name_iqrlqr    = path + filename2
    file_name_iqrselqr  = path + filename3

    selqr_data      = np.loadtxt(file_name_selqr)
    iqrlqr_data     = np.loadtxt(file_name_iqrlqr)
    iqrselqr_data   = np.loadtxt(file_name_iqrselqr)

    positions = [[0, -13.5], [10, -5.0], [-9.5, -5.0], [-2, 3], [8, 7], [11, 20], [-12, 8], [-11, 21], [-1, 16],
                 [-11, -19], [10 + np.sqrt(2.0), -15 - np.sqrt(2.0)]]

    # positions = [[-26, 17], [-11, 9], [-14, -1], [-4, -6], [-21, -10], [6, 4], [1, 17],
    #               [14, -1], [25, 0], [1, -20], [10, 20], [-12, -20], [-11, 22], [18, -13],
    #               [15, 15], [-25, 5], [8, -10]]

    f, ((fig1, fig2), (fig3, fig4)) = plt.subplots(2, 2, sharex=False, sharey=False)

    # ax1 = fig1.gca()

    for position in positions:
        circle = plt.Circle((position[0], position[1]), 2.0, color='r')
        fig1.add_artist(circle)

    fig1.plot(selqr_data[:, 0], selqr_data[:, 1], label='SELQR')
    fig1.plot(iqrlqr_data[:, 0], iqrlqr_data[:, 1], label='iQRLQR')
    fig1.plot(iqrselqr_data[:, 0], iqrselqr_data[:, 1], label='iQRSELQR')
    fig1.set_xlim([-30, 30])
    fig1.set_ylim([-30, 30])
    fig1.grid()
    #fig1.legend(loc=0)

    # Plot linear velocity
    fig2.plot(selqr_data[:, 3], label='SELQR')
    fig2.plot(iqrlqr_data[:, 3], label='iQRLQR')
    fig2.plot(iqrselqr_data[:, 3], label='iQRSELQR')
    fig2.set_xlabel("$k$")
    fig2.set_ylabel("$v$")
    fig2.grid()
    fig2.legend(loc=0)

    # Plot angular velocity
    fig3.plot(selqr_data[:, 4], label='SELQR')
    fig3.plot(iqrlqr_data[:, 4], label='iQRLQR')
    fig3.plot(iqrselqr_data[:, 4], label='iQRSELQR')
    fig3.set_xlabel("$k$")
    fig3.set_ylabel("$\omega$")
    fig3.grid()
    fig3.legend(loc=0)


    plt.show()