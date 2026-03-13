import numpy as np
from numpy.linalg import norm
from scipy.integrate import solve_ivp
import numba
import pygmo as pg

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
    
    def target(self, xx0, xxf_fun, tof_max, N=5, tof_guess=None, uus_guess=None, time_opt=False, is_mbh=True, verbose=True):
        solsail_udp = solsail.udp(self, N, xx0, xxf_fun, tof_max, time_opt)
        prob = pg.problem(solsail_udp)
        prob.c_tol = 1e-4

        sqp = pg.nlopt(solver="slsqp")
        sqp.xtol_rel = 1e-4
        sqp.xtol_abs = 1e-4
        if is_mbh:
            if time_opt:
                sqp.maxeval = 400
                sqp.maxtime = 120  # If this is reached, then something is really wrong
            else:
                sqp.maxeval = 200
                sqp.maxtime = 60

        local_algo = pg.algorithm(sqp)
        if is_mbh:
            local_algo.set_verbosity(0)
        
        algo = None
        if is_mbh:
            algo = pg.algorithm(pg.mbh(local_algo, stop=5, perturb=0.5))
        else:
            algo = local_algo

        if verbose:
            algo.set_verbosity(1)
        else:
            algo.set_verbosity(0)

        pop = None
        if tof_guess is None or uus_guess is None:
            pop = pg.population(prob, size=1)
        else:
            tof_guess_sc = tof_guess / solsail_udp.tof_unit
            uus_guess_vec = np.reshape(uus_guess, (2*N))
            x_guess = np.hstack((tof_guess_sc, uus_guess_vec))

            pop = pg.population(prob)
            pop.push_back(x=x_guess)
        
        pop = algo.evolve(pop)
        sol_f = pop.get_f()[0]  # Not really useful for anything because it's scaled
        sol_x = pop.get_x()[0]

        sol_tof = sol_x[0] * solsail_udp.tof_unit
        sol_uus = sol_x[1:].reshape((N, 2))

        return sol_tof, sol_uus, sol_f


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
    
    # Only used for Pygmo to optimize
    class udp:
        def __init__(self, ss, N, xx0, xxf_fun, tof_max, time_opt):
            self.ss = ss
            self.N = N
            self.xx0 = xx0
            self.xxf_fun = xxf_fun
            self.tof_max = tof_max
            self.time_opt = time_opt

            # Scaling parameters
            self.tof_unit = tof_max
            self.pos_unit = norm(xx0[0:3])
            self.vel_unit = norm(xx0[3:6])

        # x = [tof, phi0, th0, ph1, th1, ..., phiN-1, thN-1]
        def fitness(self, x):
            tof = x[0] * self.tof_unit  # tof is scaled
            uus = x[1:].reshape((self.N, 2))

            xxs, _ = self.ss.propagate(self.xx0, tof, uus)

            if xxs is None:
                if self.time_opt:
                    return [np.nan, np.nan]
                else:
                    return [np.nan]

            xxf = xxs[-1]

            xxf_targ = self.xxf_fun(tof)

            # Optimizes only position, or position and velocity
            constr = 0
            if xxf_targ.size == 3:
                constr = norm(xxf_targ - xxf[0:3]) / self.pos_unit
            elif xxf_targ.size == 6:
                constr = norm(xxf_targ[0:3] - xxf[0:3]) / self.pos_unit + \
                      norm(xxf_targ[3:6] - xxf[3:6]) / self.vel_unit

            if self.time_opt:
                return [x[0], constr]
            else:
                return [constr]
        
        def get_bounds(self):
            lb_tof = [1 / self.tof_unit]
            ub_tof = [self.tof_max / self.tof_unit]

            lb_uu = [-np.pi/2] * 2 * self.N
            ub_uu = [np.pi/2] * 2 * self.N

            lb = np.hstack((lb_tof, lb_uu))
            ub = np.hstack((ub_tof, ub_uu))

            return (lb, ub)
        
        # Number of equality constraints
        def get_nec(self):
            if self.time_opt:
                return 1
            else:
                return 0
        
        # Number of inequality constraints
        def get_nic(self):
            return 0

        # Supply numerical gradient
        def gradient(self, x):
            return pg.estimate_gradient(lambda x: self.fitness(x), x)