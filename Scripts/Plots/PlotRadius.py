import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import Ellipse
import matplotlib

if __name__ == '__main__':

    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['ps.fonttype'] = 42

    path0   = np.loadtxt("SELQRPathsRadius6.txt")
    path1   = np.loadtxt("iQRSELQRPathsRadius6.txt")
    radius1 = np.loadtxt("RadiusPathsRadius6.txt")


    positions = [[0, -13.5], [10, -5.0], [-9.5, -5.0], [-2, 3], [8, 7], [11, 20], [-12, 8], [-11, 21], [-1, 16],
                 [-11, -19], [10 + np.sqrt(2.0), -15 - np.sqrt(2.0)]]

    positions = [[0, -13.5], [10, -5.0], [-9.5, -5.0], [-2, 3], [8, 7], [-12, 8], [-11, 21], [-1, 16],
                 [11, 20]]

    # positions = [[-26, 17], [-11, 9], [-14, -1], [-4, -6], [-21, -10], [6, 4], [1, 17],
    #               [14, -1], [25, 0], [1, -20], [10, 20], [-12, -20], [-11, 22], [18, -13],
    #               [15, 15], [-25, 5], [8, -10]]

    ax = plt.gca()

    for position in positions:
        circle = plt.Circle((position[0], position[1]), 2.0, color='r')
        ax.add_artist(circle)

    ax.plot(path0[:, 0], path0[:, 1], label='E-LQR', c='b', lw=2)
    ax.plot(path1[:, 0], path1[:, 1], label='RE-LQR', c='black', lw=2)
    ax.set_xlim([-15, 14])
    ax.set_ylim([-3, 25])
    ax.set_xlabel("$x$(m)")
    ax.set_ylabel("$y$(m)")
    ax.legend(loc=4)

    ax.grid()

    # ax.add_patch(patches.Rectangle(
    #     (-30, -30),
    #     60,
    #     60,
    #     fill=False,  # remove background,
    #     edgecolor="red"
    # )
    # )

    color1 = [0, 0.1, 1]#np.random.rand(3)
    color2 = [0, 1, 0]#np.random.rand(3)
    i =0
    for rad in radius1:
        ellipse = Ellipse(xy=(path1[i, 0], path1[i, 1]), width=2*np.sqrt(rad[0]), height=2*np.sqrt(rad[1]),
                           edgecolor='black', fc='black', alpha=0.2, lw=1)
        ax.add_patch(ellipse)

        i += 1

    ax.set_aspect('equal')
    ax.set_rasterized(True)

    plt.show()