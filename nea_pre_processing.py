import pandas as pd
import spiceypy as spice
import numpy as np

if __name__ == '__main__':
    spice.kclear()
    spice.furnsh("kernels/naif0012.tls")

    asteroid_file = "astorb.dat"
    colspecs = [
            (0, 6),
            (7, 25),
            (74, 78),
            (82, 86),
            (106, 110),
            (110, 112),
            (112, 114),
            (115, 125),
            (126, 136),
            (137, 147),
            (147, 157),
            (158, 168),
            (168, 181)
        ]
    asteroid_data = pd.read_fwf(asteroid_file, colspecs=colspecs, header=None, nrows=887103,
                        names=[
                            "id",
                            "name",
                            "code2",
                            "code4",
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
                            "code2": int,
                            "code4": int,
                            "year": int,
                            "month": int,
                            "day": int,
                            "MA": float,
                            "AOP": float,
                            "RAAN": float,
                            "INC": float,
                            "ECC": float,
                            "SMA": float
                        })

    # Remove faulty or bad measurements
    asteroid_data = asteroid_data[asteroid_data.code2 == 0]

    # Remove critical asteroid (lost or something like that)
    asteroid_data = asteroid_data[asteroid_data.code4 == 0]

    asteroid_data = asteroid_data.drop(columns=["code2", "code4"])

    # Remove that never pass by close to earth
    asteroid_data = asteroid_data[asteroid_data.SMA * (1 - asteroid_data.ECC) < 1.1]
    asteroid_data = asteroid_data[asteroid_data.SMA * (1 + asteroid_data.ECC) > 0.9]

    # Orbital period should be below 5 years so it passes by at least once
    year = 3600 * 24 * 365.25
    asteroid_data = asteroid_data[2*np.pi * np.sqrt(asteroid_data.SMA**3 / 3.964e-14) < 5 * year]

    # Save remaining asteroids in csv file
    asteroid_data.to_csv('NEA_selection.csv', header=False, index=False)