import numpy as np
from scipy.optimize import root_scalar
import pandas as pd
import spiceypy as spice
import pykep as pk

class asteroids:
    # Loads an asteroid database, see function definition for more details
    # Watch out: Requires a spice leapsecond kernel to be loaded!
    def __init__(self, file_name="NEA_selection.csv"):
        self.mu_sun = pk.MU_SUN / 1e9  # km^3 / s^2

        # Should be CSV file of asteroid data in this order per row
        raw_data = pd.read_csv(file_name, header=None,
                    names=[
                        "id",
                        "name",
                        "year",
                        "month",
                        "day",
                        "MA",
                        "AOP",
                        "RAAN",
                        "INC",
                        "ECC",
                        "SMA"
                        ],
                    dtype={
                        "id": int,
                        "name": str,
                        "year": int,   # Year, month, and day in TDT
                        "month": int,  # are the epoch of the osculating
                        "day": int,    # orbital elements in the following
                        "MA": float,   # deg, Mean Anomaly
                        "AOP": float,  # deg (J2000), Argument of Periapsis
                        "RAAN": float, # deg (J2000), Right Ascension of the Ascending Node
                        "INC": float,  # deg (J2000), Inclination
                        "ECC": float,  # Eccentricity
                        "SMA": float   # AU, semi-major axis
                        })

        self.ids = raw_data["id"].to_list()
        self.N = np.size(self.ids)
        self.names = raw_data["name"].to_list()

        # Convert TDT dates to ephemeris time
        date2et = lambda data: spice.str2et("{}-{}-{} TDT".format(data.year, data.month, data.day))
        self.et0s = raw_data.apply(date2et, axis=1)

        # Transform orbital elements to km and rad, and from J2000 to ECLIPJ2000
        # Actually turns out it's already in ECLIPJ2000 and the Lowell Observatory lied to me
        oe_idx = ["SMA", "ECC", "INC", "RAAN", "AOP", "MA"]
        oe = raw_data[oe_idx].to_numpy()
        oe[:, 0] = oe[:, 0] * (pk.AU / 1e3)
        oe[:, 2:6] = np.deg2rad(oe[:, 2:6])
        self.kep = oe

        # Orbital period
        self.T = 2*np.pi * np.sqrt(self.kep[:, 0]**3 / self.mu_sun)
        # Mean motion
        self.n = 2*np.pi / self.T

    # Returns the state of one or more asteroids at one or more epochs
    # The propagation is Keplerian so this should only be used for very preliminary trajectory design
    # idx: Single index or list of indices of asteroids whose states should be returned
    # et: Ephemeris time at which asteroid state is returned
    #     If scalar: The output is a 2d array of each idx asteroid's state for this epoch
    #     If 1d: Output is a 3d array of each idx asteroid's state for each epoch in et
    #     If 2d: Output is a 3d array of each idx asteroid's states for the epochs in its row in et.
    #            First dim must have same size as idx, second dim is number of epochs for each asteroid
    # Output: Cartesian asteroid states in ECLIPJ2000. 1, 2, or 3 dim array. If 3d then first dim is same as idx
    #         second dim is same as second dim of et, and third dim is 6
    def get_state(self, idx, et):
        N_idx = np.size(idx)
        N_t = 1
        et_mat = et

        if np.ndim(et) == 0:
            et_mat = np.tile(et, (N_idx, 1))
        elif np.ndim(et) == 1:
            et_mat = np.tile(et, (N_idx, 1))
            N_t = np.size(et)
        elif np.ndim(et) == 2:
            N_t = np.shape(et)[1]

        et0 = np.tile(self.et0s[idx], (N_t, 1)).T
        dt = et_mat - et0

        keps = np.atleast_3d(self.kep[idx])
        keps = np.moveaxis(keps, 2, 1)
        keps = np.tile(keps, (1, N_t, 1))
        keps[:, :, 5] += dt * np.tile(self.n[idx], (N_t, 1)).T  # Update mean anomaly
        keps[:, :, 5] = keps[:, :, 5] % (2 * np.pi)  # Wrap to 2 pi

        xxs = np.empty((N_idx, N_t, 6))
        for t_idx in range(N_t):
            # Solve kepler's equation to get the eccentric anomaly
            keps[:, t_idx, 5] = asteroids.__mean2eccentric(keps[:, t_idx, 5], keps[:, t_idx, 1])

            for i in range(N_idx):
                rri, vvi = pk.par2ic(keps[i, t_idx, :], self.mu_sun)
                xxs[i, t_idx, :] = np.hstack((rri, vvi))

        return xxs.squeeze()

    # Solves the Kepler equation
    def __mean2eccentric(M, e):
        N = M.size

        E_guess = M + e * np.sin(M)
        E = np.empty(N)

        # Scalar loop is faster than computing it in vectorized form
        for i in range(N):
            kep_eq = lambda E: E - e[i] * np.sin(E) - M[i]
            sol = root_scalar(kep_eq, x0=E_guess[i], method="secant")
            E[i] = sol.root

        return E % (2 * np.pi)