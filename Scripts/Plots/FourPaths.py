import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib

if __name__ == "__main__":

    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['ps.fonttype'] = 42

    path         = "/Volumes/Aslan/OptimalControl/Configurations/Map3/"
    filename1    = "SELQREpsilon8.txt"
    filename2    = "iQRSELQREpsilon7.txt"
    filename3    = "iQRSELQREpsilon8.txt"

    file_name_selqr     = path + filename1
    file_name_iqrlqr    = path + filename2
    file_name_iqrselqr  = path + filename3

    selqr_data      = np.loadtxt(file_name_selqr)
    iqrlqr_data     = np.loadtxt(file_name_iqrlqr)
    iqrselqr_data   = np.loadtxt(file_name_iqrselqr)

    positions = [[0, -13.5], [10, -5.0], [-9.5, -5.0], [-2, 3], [8, 7], [11, 20], [-12, 8], [-11, 21], [-1, 16],
                 [-11, -19], [10 + np.sqrt(2.0), -15 - np.sqrt(2.0)]]

    positions = [[-26, 17], [-11, 9], [-14, -1], [-4, -6], [-21, -10], [6, 4], [1, 17],
                 [14, -1], [25, 0], [1, -20], [10, 20], [-12, -20], [-11, 22], [18, -13],
                 [15, 15], [-25, 5], [8, -10]]

    ax = plt.gca()

    for position in positions:
        circle = plt.Circle((position[0], position[1]), 2.0, color='r')
        ax.add_artist(circle)

    ax.plot(selqr_data[:, 0],   selqr_data[:, 1],    label='E-LQR')
    ax.plot(iqrlqr_data[:, 0],  iqrlqr_data[:, 1],   label='RE-LQR ($\\epsilon=0.1$)')
    ax.plot(iqrselqr_data[:, 0],iqrselqr_data[:, 1], label='RE-LQR ($\\epsilon=0.001$)')
    ax.set_xlim([-31, 31])
    ax.set_ylim([-31, 31])
    ax.legend(loc=4, prop={'size': 10})
    ax.grid()

    ax.add_patch(patches.Rectangle(
        (-30, -30),
        60,
        60,
        fill=False,  # remove background,
        edgecolor="red"
    )
    )

    ax.set_aspect('equal')
    plt.savefig("Map3_e3.pdf", bbox_inches='tight')
    plt.show()