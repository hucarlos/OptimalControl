import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib

if __name__ == '__main__':

    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['ps.fonttype'] = 42

    positions1 = [[0, -13.5], [10, -5.0], [-9.5, -5.0], [-2, 3], [8, 7], [11, 20], [-12, 8], [-11, 21], [-1, 16],
                 [-11, -19], [10 + np.sqrt(2.0), -15 - np.sqrt(2.0)]]


    positions3 = [[-26, 17], [-11, 9], [-14, -1],[-4,-6],[-21,-10], [6,4], [1,17],
    [14,-1],[25,0], [1,-20], [10,20], [-12,-20], [-11,22], [18,-13],
    [15,15], [-25,5],[8,-10]]

    ax = plt.gca()

    for position in positions3:
        circle = plt.Circle((position[0], position[1]), 2.0, color='r')
        ax.add_artist(circle)

    # for position in positions3:
    #     circle = plt.Circle((position[0], position[1]), 2.0, color='r')
    #     fig3.add_artist(circle)

    # ax.add_patch(patches.Rectangle(
    #     (-25, 25),
    #     60,
    #     60,
    #     fill=False,      # remove background,
    #     edgecolor="red"
    # )
    # )


    ax.set_xlim([-30, 28])
    ax.set_ylim([-25, 25])
    ax.set_title('Map 2')
    ax.set_xlabel("$x$ (cm)")
    ax.set_ylabel("$y$ (cm)")


    ax.set_aspect('equal')
    plt.savefig("Map2.pdf", bbox_inches='tight')

    plt.show()

