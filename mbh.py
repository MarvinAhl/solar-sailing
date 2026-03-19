import numpy as np
from scipy.optimize import minimize

# Implementation of Monotonic Basin Hopping (MBH) with scipy's minimize as local solver.
# Only unconstrained optimization possible for now.
# Inputs:
#   f: A callable, the objective function
#   bounds: List of tuples (lb, ub) for each element of x
#   x0: Initial guess, if none supplied, it's chosen randomly
#   stop: Algorithm terminated after this number of consecutive non-improvements
#   pert: Scalar or vector of shape of x. x perturbed by up to this much after each local iteration
#   method: String, can be any method that's allowed in scipy's minimized
#   local_ftol: Convergence tolerance for local solver
#   jac: Callable, string or None. Objective jacobian, if none supplied, 2-point finite differences used
#   verbose: Bool, whether to print status
# Output:
#   best_x: Decision vector of best solution
#   best_f: Objective value of best solution
def mbh(f, bounds, x0=None, stop=5, pert=1., method='SLSQP', local_ftol=1e-8, jac=None, verbose=True):
    lb = np.array([b[0] for b in bounds])
    ub = np.array([b[1] for b in bounds])

    if jac is None:
        jac = '2-point'

    N_no_change = 0
    N_f_evals = 0

    best_f = None
    best_x = None

    candidate_x = None
    if x0 is None:
        candidate_x = np.random.rand(len(bounds)) * (ub - lb) + lb
    else:
        candidate_x = np.array(x0)

    if verbose:
        print()
        print("===========================")
        print("  Monotonic Basin Hopping")
        print("fevals   f             stop")
        print("---------------------------")

    while N_no_change < stop:

        # Perform local optimization
        res = minimize(f, candidate_x, method=method, jac=jac, options={'ftol': local_ftol, 'disp': False}, bounds=bounds)

        # Update best solution
        if best_f is None or res.fun < best_f:
            best_x = res.x
            best_f = res.fun
            N_no_change = 0
        else:
            N_no_change += 1
        
        N_f_evals += res.nfev

        del res
        
        # Perturb best solution
        update = np.random.rand(len(bounds)) * 2 - 1
        update *= pert
        candidate_x = np.clip(best_x + update, lb, ub)

        # Plot a bit
        if verbose:
            print("{fevals:6n}   {f:.5e}   {stop:4n}".format(fevals=N_f_evals, f=best_f, stop=N_no_change))

    if verbose:
        print("MBH Terminated")
        print("===========================")

    return best_x, best_f