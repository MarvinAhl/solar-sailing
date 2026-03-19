import numpy as np
from numpy.linalg import norm
from scipy.integrate import solve_ivp
import numba
from mbh import mbh

class solsail:
    def __init__(self, mu, I0, r0, Am, c):
        self.mu = mu
        self.I0 = I0
        self.r0 = r0
        self.Am = Am
        self.c = c
    
    # Propagate solar sailing with piecewise constant control
    # xx0: initial state [r0, r1, r2, v0, v1, v2]
    # tof: time of flight
    # uus: Controls [[phi0, th0], [phi1, th1], ..., [phiN-1, thN-1]]
    # n:   Number of states to output for each control segment, default is 1,
    #      at least 10 recommended for plotting
    def propagate(self, xx0, tof, uus, n=1):
        N = uus.shape[0]  # Get number of constant control segments

        # Starting times of the control sections and final time
        tus = np.linspace(0, tof, N+1)

        # Output times
        ts = np.empty(n * N + 1)
        ts[0] = 0

        # Output states
        xxs = np.empty((n * N + 1, 6))
        xxs[0] = xx0

        # Propagate each constant control segment forward and save states
        for i in range(N):
            t0i = tus[i]
            t1i = tus[i+1]
            tsi = np.linspace(t0i, t1i, n+1)
            ts[i*n+1 : (i+1)*n+1] = tsi[1:]

            uui = uus[i]
            xx0i = xxs[i*n]

            if not np.all(np.isfinite(xx0i)):
                print("Warning: Non-finite y0 of solve_ivp in solsail.propagate()")
                return None, None

            # Numerical propagation
            soli = solve_ivp(lambda t, xx: solsail.__dynamics(t, xx, uui, self.mu, self.I0, self.r0, self.Am, self.c),
                t_span=(t0i, t1i), y0=xx0i, method='LSODA', t_eval=tsi, rtol=1e-6, atol=1e-6)
            
            xxs[i*n+1 : (i+1)*n+1] = soli.y[:, 1:].T
        
        # Return times and states
        return xxs, ts
    
    def target(self, xx0, xxf_fun, tof_max, N=5, tof_guess=None, uus_guess=None, verbose=True):
        # Scaling factors
        tof_unit = tof_max
        pos_unit = norm(xx0[0:3])
        vel_unit = norm(xx0[3:6])

        # Initial guess creation (if anything supplied)
        x_guess = None
        if tof_guess is not None and uus_guess is not None:
            tof_guess_sc = tof_guess / tof_unit
            uus_guess_vec = np.reshape(uus_guess, (2*N))
            x_guess = np.hstack((tof_guess_sc, uus_guess_vec))
        
        # Define bounds
        lb_tof = [1 / tof_unit]
        ub_tof = [tof_max / tof_unit]
        lb_uu = [-np.pi/2] * 2 * N
        ub_uu = [np.pi/2] * 2 * N
        lb = np.hstack((lb_tof, lb_uu))
        ub = np.hstack((ub_tof, ub_uu))

        bounds = [(l, u) for l, u in zip(lb, ub)]

        # Optimize
        f = lambda x: self.fitness(x, xx0, xxf_fun, N, tof_unit, pos_unit, vel_unit)
        sol_x, sol_f = mbh(f, bounds, x_guess, stop=5, pert=0.5, local_ftol=1e-5, verbose=verbose)

        # Extract and save results
        sol_tof = sol_x[0] * tof_unit
        sol_uus = sol_x[1:].reshape((N, 2))

        return sol_tof, sol_uus, sol_f
    
    def fitness(self, x, xx0, xxf_fun, N, tof_unit, pos_unit, vel_unit):
        tof = x[0] * tof_unit  # tof is scaled
        uus = x[1:].reshape((N, 2))

        # Propagate all piecewise constant solar sailing segments and only save final state
        ts = np.linspace(0, tof, N+1)
        xxf = xx0
        for i in range(N):
            uui = uus[i]
            soli = solve_ivp(lambda t, xx: solsail.__dynamics(t, xx, uui, self.mu, self.I0, self.r0, self.Am, self.c),
                                t_span=(ts[i], ts[i+1]), y0=xxf, method='LSODA', t_eval=[ts[i+1]], rtol=1e-6, atol=1e-6)
            xxf = soli.y.squeeze().copy()
            del soli

        xxf_targ = np.array(xxf_fun(tof))

        # Target only position, or position and velocity, depending on what the target function outputs
        constr = 0
        if xxf_targ.size == 3:
            constr = norm(xxf_targ - xxf[0:3]) / pos_unit
        elif xxf_targ.size == 6:
            constr = norm(xxf_targ[0:3] - xxf[0:3]) / pos_unit + \
                    norm(xxf_targ[3:6] - xxf[3:6]) / vel_unit

        # This is the optimization objective
        return constr


    # Constant control Solar Sailing dynamics considering sun 2-body gravity and SRP
    @numba.njit(numba.float64[:](numba.float64, numba.float64[:],
                                numba.float64[:], numba.float64,
                                numba.float64, numba.float64,
                                numba.float64, numba.float64))
    def __dynamics(t, xx, uu, mu, I0, r0, Am, c):
        phi = uu[0]
        th = uu[1]

        rr = xx[0:3]
        vv = xx[3:6]
        r = norm(rr)

        # Compute 2 body acceleration
        aa_2b = -mu * rr / r**3

        # Compute SRP in RTH frame
        a_srp = 2 * I0/c * (r0 / r)**2 * Am * np.cos(phi) * np.cos(th)
        aa_srp_rth = np.array([np.cos(phi) * np.cos(th), np.sin(phi) * np.cos(th), np.sin(th)]) * a_srp
        
        # Compute transformation from RTH to inertial frame
        r_hat = rr / r
        h_hat = np.cross(rr, vv)
        h_hat = h_hat / norm(h_hat)
        t_hat = np.cross(h_hat, r_hat)
        T_rth = np.vstack((r_hat, t_hat, h_hat)).T

        # Transform SRP to inertial frame
        aa_srp = T_rth @ aa_srp_rth

        # Output results
        rr_dot = vv
        vv_dot = aa_2b + aa_srp
        return np.concat((rr_dot, vv_dot))