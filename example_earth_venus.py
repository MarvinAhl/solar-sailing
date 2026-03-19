import numpy as np
from numpy.linalg import norm
import spiceypy as spice
import pykep as pk

from solsail import solsail

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Lead spice kernels
spice.kclear()
spice.furnsh("./kernels/metakernel.tm")

if __name__ == '__main__':
    I0 = 1361.  # W/m^2 = kg/s^3
    Am = 32 / 5 * 1e-6  # km^2/kg, data from lightsail 2
    Am *= 5  # Make it a bit stronger
    c = 299792.458  # km/s
    au = pk.AU / 1e3
    mu_sun = spice.bodvrd("SUN", "GM", 1)[1][0]

    # Start at Earth (v_inf = 0) in the beginning of 2031
    t0 = spice.str2et("2031-01-01T00:00:00")
    xx0, _ = spice.spkezr("EARTH", t0, "ECLIPJ2000", "NONE", "SSB")
    xxf_fun = lambda tof: spice.spkezr("VENUS", tof + t0, "ECLIPJ2000", "NONE", "SSB")[0]
    # Maximum time of flight is 5 years
    tof_max = pk.DAY2SEC * 365.25 * 4
    # Number of constant control segments
    N = 10

    # Create solar sailing object
    ss = solsail(mu_sun, I0, au, Am, c)

    tof_guess = 87520508.8860183
    uus_guess = np.array([[ 0.00275677, 1.00138217],
                          [-1.05295792, 0.21831993],
                          [ 0.61175096, 0.35996422],
                          [-0.48579461, 1.33658835],
                          [-0.53209007, 0.72154463],
                          [-0.53819903, 0.21944811],
                          [-0.60633416, 0.2750143 ],
                          [-0.24237431, -0.78986399],
                          [-0.19269383, -0.91480272],
                          [-0.97670402, 0.71398152]])

    # Solve Solar Sailing shooting problem
    sol_tof, sol_uus, _ = ss.target(xx0, xxf_fun, tof_max, N, tof_guess=tof_guess, uus_guess=uus_guess)
    print("Solution:")
    print(sol_tof)
    print(sol_uus)

    # Propagate the solution for plotting
    sol_xxs, _ = ss.propagate(xx0, sol_tof, sol_uus, n=100)
    sol_xxf = sol_xxs[-1]
    
    # Errors
    xx1_venus, _ = spice.spkezr("VENUS", t0 + sol_tof, "ECLIPJ2000", "NONE", "SSB")
    print("Time of flight: {} years".format(sol_tof / 3600/24/365.25))
    print("Position error: {} km".format(norm(sol_xxf[0:3] - xx1_venus[0:3])))
    print("Velocity error: {} km/s".format(norm(sol_xxf[3:6] - xx1_venus[3:6])))

    # Plot result 3D
    fig3d = plt.figure(figsize=(6, 5))
    ax3d = fig3d.add_subplot(111, projection='3d', proj_type='ortho')

    # Plot Sun, Earth, and Venus
    ax3d.plot(0, 0, 0, 'y.', markersize=20)
    ts_earth = np.linspace(t0, t0 + pk.DAY2SEC/pk.DAY2YEAR, 300)
    rrs_earth, _ = spice.spkpos("EARTH", ts_earth, "ECLIPJ2000", "NONE", "SSB")
    ts_venus = np.linspace(t0, t0 + 225*pk.DAY2SEC, 600)
    rrs_venus, _ = spice.spkpos("VENUS", ts_venus, "ECLIPJ2000", "NONE", "SSB")
    ax3d.plot(rrs_earth[:, 0], rrs_earth[:, 1], rrs_earth[:, 2], 'b--')
    ax3d.plot(rrs_venus[:, 0], rrs_venus[:, 1], rrs_venus[:, 2], 'r--')
    ax3d.plot(rrs_earth[0, 0], rrs_earth[0, 1], rrs_earth[0, 2], 'b.', markersize=8)
    ax3d.plot(xx1_venus[0], xx1_venus[1], xx1_venus[2], 'r.', markersize=8)

    # Plot trajectory
    ax3d.plot(sol_xxs[:, 0], sol_xxs[:, 1], sol_xxs[:, 2], 'k-')

    ax3d.axis('equal')
    ax3d.set_xlabel('x / km')
    ax3d.set_ylabel('y / km')
    ax3d.set_zlabel('z / km')

    # Plot results 2D

    fig2d = plt.figure(figsize=(5, 4))
    ax2d = fig2d.add_subplot(111)
    
    ax2d.plot(0, 0, 'y.', markersize=20)
    ax2d.plot(rrs_earth[:, 0], rrs_earth[:, 1], 'b--')
    ax2d.plot(rrs_venus[:, 0], rrs_venus[:, 1], 'r--')
    ax2d.plot(rrs_earth[0, 0], rrs_earth[0, 1], 'b.', markersize=8)
    ax2d.plot(xx1_venus[0], xx1_venus[1], 'r.', markersize=8)

    # Plot trajectory
    ax2d.plot(sol_xxs[:, 0], sol_xxs[:, 1], 'k-')

    ax2d.axis('equal')
    ax2d.set_xlabel('x / km')
    ax2d.set_ylabel('y / km')

    # Plot control
    figcon = plt.figure(figsize=(5, 4))
    axphi = figcon.add_subplot(211)

    sol_t_con = np.linspace(0, sol_tof, N+1) / 3600 / 24
    
    axphi.step(sol_t_con, np.hstack((sol_uus[:, 0], sol_uus[-1, 0])) * 180 / np.pi, 'k-', where='post')
    axphi.set_ylabel("phi / deg")
    axphi.set_xlim((sol_t_con[0], sol_t_con[-1]))
    axphi.set_ylim((-90, 90))

    axth = figcon.add_subplot(212)

    axth.step(sol_t_con, np.hstack((sol_uus[:, 1], sol_uus[-1, 1])) * 180 / np.pi, 'k-', where='post')
    axth.set_xlabel("t / d")
    axth.set_ylabel("theta / deg")
    axth.set_xlim((sol_t_con[0], sol_t_con[-1]))
    axth.set_ylim((-90, 90))

    plt.show()

    pass