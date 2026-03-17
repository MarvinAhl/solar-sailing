import numpy as np
from numpy.linalg import norm
from scipy.integrate import solve_ivp
import numba
import spiceypy as spice
import pykep as pk

from solsail import solsail
from asteroids import asteroids

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

if __name__ == '__main__':
    spice.kclear()
    spice.furnsh("kernels/metakernel.tm")

    # Initialize solar sail
    I0 = 1361.  # W/m^2 = kg/s^3
    Am = 32 / 5 * 1e-6  # km^2/kg, data from lightsail 2
    Am *= 5  # Scaled a bit to make it more manoeuvreable
    c = 299792.458  # km/s
    au = pk.AU / 1e3
    mu_sun = pk.MU_SUN / 1e9

    ss = solsail(mu_sun, I0, au, Am, c)

    # Initialize asteroid database
    neas = asteroids("NEA_selection.csv")

    # Find best solution in search tree
    search_tree = np.load("search_tree.npy")

    max_level = search_tree['level'][-1]
    max_level_num = np.sum(search_tree['level'] == max_level)
    select_sol = 0  # 0 for best, 1 for second best, ...
    opt_idx = np.where(search_tree['level'] == 3)[0][select_sol]  # Tree levels ordered by time

    print("Maximum number of asteroid flybys: {} ({} solutions with that number)".format(max_level, max_level_num))

    t0 = search_tree[0][1]
    tf_opt = search_tree[opt_idx][1]
    tof_opt = tf_opt - t0

    print("Departure: " + spice.timout(t0, 'YYYY-MM-DD HR:MN:SC'))
    print("Optimal Arrival: " + spice.timout(tf_opt, 'YYYY-MM-DD HR:MN:SC'))
    print("Optimal ToF: {} years".format(tof_opt / (pk.DAY2SEC * 365.25)))

    # Backtrack optimal solution
    opt_branch = []
    opt_branch.append(search_tree[opt_idx])
    while opt_branch[-1][5] > 0:  # While there is a parent
        opt_branch.append(search_tree[opt_branch[-1][4]])  # Append parent to list
    
    # Reverse list order
    opt_branch = opt_branch[::-1]

    print("Asteroids passed:")

    # Plot the transfer
    fig = plt.figure(figsize=(6, 5))
    ax = fig.add_subplot(111, projection='3d', proj_type='ortho')

    ax.plot(0, 0, 0, 'y.', markersize=20)

    xxs = []
    for nodei, nodeip1 in zip(opt_branch[:-1], opt_branch[1:]):
        xx0i = nodei[0]
        tofi = nodeip1[1] - nodei[1]
        uusi = nodeip1[2]
        xxsi, _ = ss.propagate(xx0i, tofi, uusi, n=20)
        if np.size(xxs) == 0:
            xxs = xxsi.copy()
        else:
            xxs = np.vstack((xxs, xxsi))

        # Plot asteroid orbit
        ast_idx_i = nodeip1[3]
        Ti = neas.T[ast_idx_i]
        ast_t0 = nodeip1[1]
        ast_t1 = ast_t0 + Ti
        ast_ts = np.linspace(ast_t0, ast_t1, 1000)
        ast_xxs = neas.get_state(idx=ast_idx_i, et=ast_ts)
        ast_plot = ax.plot(ast_xxs[:, 0], ast_xxs[:, 1], ast_xxs[:, 2], linewidth=1, label=neas.names[ast_idx_i])
        ax.plot(ast_xxs[0, 0], ast_xxs[0, 1], ast_xxs[0, 2], '.', markersize=7, color=ast_plot[0].get_color(), zorder=100)

        print(neas.names[ast_idx_i] + ", idx: {}".format(ast_idx_i))
        print("  Distance: {} km".format(norm(ast_xxs[0, 0:3] - xxsi[-1, 0:3])))
        print("  Rel. velocity: {} km/s".format(norm(ast_xxs[0, 4:6] - xxsi[-1, 4:6])))

    ax.plot(xxs[:, 0], xxs[:, 1], xxs[:, 2], 'k-')

    ax.axis('equal')
    ax.set_xlabel('x / km')
    ax.set_ylabel('y / km')
    ax.set_zlabel('z / km')

    ax.legend()

    plt.show()

    pass