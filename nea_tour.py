import numpy as np
from numpy.linalg import norm
from scipy.integrate import solve_ivp
import numba
import spiceypy as spice
import pykep as pk
import time

from solsail import solsail
from asteroids import asteroids

spice.kclear()
spice.furnsh("kernels/metakernel.tm")

@numba.njit(numba.float64[:](numba.float64, numba.float64[:], numba.float64))
def kep_dyn(t, xx, mu):
    rr = xx[0:3]
    vv = xx[3:6]
    r = norm(rr)

    xx_dot = np.empty_like(xx)
    xx_dot[0:3] = vv
    xx_dot[3:6] = -mu * rr / r**3

    return xx_dot

def close_approach(t, xx, asts, ast_idx):
    ast_xx = asts.get_state(idx=ast_idx, et=t)

    drr = xx[0:3] - ast_xx[0:3]
    dvv = xx[3:6] - ast_xx[3:6]

    return np.dot(drr, dvv)

if __name__ == '__main__':
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

    # Start at Earth (v_inf = 0) in the beginning of 2030
    t0 = spice.str2et("2030-01-01T00:00:00")
    xx0, _ = spice.spkezr("EARTH", t0, "ECLIPJ2000", "NONE", "SSB")
    # Find asteroid encounters for no more than 5 years
    tf = spice.str2et("2035-01-01T00:00:00")

    # Propagation time for looking for encounters for a single arc is 6 months
    t_search = pk.DAY2SEC * 30 * 6
    # Number of grid points for each arc
    Nt = 50  # 20 was already enough but we wanna be on the safe side here

    # Number of piecewise constant control segments per solar sailing arc
    N_control = 3

    # How many nodes are retained in each level in beach search
    beam_width = 3

    verbose_optim = False

    # The beam search lists
    # Each node is a tuple of (current xx, current t, uu to get there, asteroid id, parent node, tree level)
    search_tree = []  # Fully explored nodes
    open_list = []  # Nodes to be expanded in next loop iteration
    earth_node = (xx0, t0, np.zeros((N_control, 2)), -1, -1, 0)
    open_list.append(earth_node)
    candidate_list = []  # Current node level, best L nodes become the next open_list

    tree_level = 0

    with open("log.txt", 'a') as log_file:
        print("Starting beam search ...", file=log_file)
        code_start_t = time.time()
        while len(open_list) > 0:
            tree_level += 1
            print("Starting tree level {}".format(tree_level), file=log_file, flush=True)
            
            for open_idx, node in enumerate(open_list):
                print("Expanding open list node {} / {} (tree level {})".format(open_idx + 1, len(open_list), tree_level), file=log_file, flush=True)

                xx0i = node[0]
                t0i = node[1]
                ast_i = node[3]
                tree_level_i = node[5]

                t1i = t0i + t_search
                tsi = np.linspace(t0i, t1i, Nt)
                dtsi = tsi - t0i

                # Propagate for one year with 100 time points 
                grid_sol_i = solve_ivp(lambda t, xx: kep_dyn(t, xx, mu_sun), (t0i, t1i), xx0i,
                                method='LSODA', t_eval=tsi, rtol=1e-6, atol=1e-6)
                xxsi = grid_sol_i.y.T

                # Propagate asteroids as well
                neas_xxsi = neas.get_state(idx=np.arange(neas.N), et=tsi)

                # Count number of asteroid close approaches from grid search

                # Maximum possible acceleration with solar sail
                ss_max_acc_i = 2 * I0/c * (au / norm(xx0i[0:3]))**2 * Am
                # Approximation for reachable distances at each time with solar sails
                drs_reach_i = ss_max_acc_i/2 * dtsi**2
                drs_reach_i /= 0.7  # Relax the distance requirement a bit, tests show it doesn't add much in practice though

                cas_idxs_i = []
                grid_cas_points_tot_i = 0  # Total number of collision risks
                for nea_idx, nea_xxs in enumerate(neas_xxsi):
                    if ast_i == nea_idx:
                        continue  # Don't go to same asteroid twice in a row
                    drs = norm(xxsi[:, 0:3] - nea_xxs[:, 0:3], axis=1)
                    reachable_et_idx = np.where(drs < drs_reach_i)

                    if np.any(reachable_et_idx):
                        cas_idxs_i.append(nea_idx)
                        grid_cas_points_tot_i += np.size(reachable_et_idx)

                N_ca_i = np.size(cas_idxs_i)  # Number of asteroids with close approaches
                print("Initial number of close approache indexes: {}".format(N_ca_i), file=log_file)
                print("(Avg grid triggers per idx: {})".format(grid_cas_points_tot_i / np.max([N_ca_i, 1])), file=log_file)

                # Propagate again but look for close approaches with the previously found candidates

                # Create close approach event functions
                ca_events_i = []
                for ca_idx in cas_idxs_i:
                    ca_event_i = lambda t, xx, ca_idx_t=ca_idx: close_approach(t, xx, neas, ca_idx_t)
                    ca_event_i.terminal = False
                    ca_event_i.direction = 1
                    ca_events_i.append(ca_event_i)

                # Propagate with close approach events
                ca_sol_i = solve_ivp(lambda t, xx: kep_dyn(t, xx, mu_sun), (t0i, t1i), xx0i,
                                method='LSODA', rtol=1e-6, atol=1e-6, events=tuple(ca_events_i))
                
                # Extract the close approach times and distances
                cas_i = []
                for i, ca_ts in enumerate(ca_sol_i.t_events):
                    if np.size(ca_ts) == 0:
                        continue  # No events detected for this
                    ca_t = ca_ts[0]  # Take the first close approach with each asteroid
                    ca_ast_id = cas_idxs_i[i]  # Get ID of close approach asteroid

                    # Get spacecraft and asteroid state at that approach time
                    ca_ss_xx = ca_sol_i.y_events[i][0]
                    ca_ast_xx = neas.get_state(ca_ast_id, ca_t)

                    # See if distance within reach (applying max acceleration heuristic)
                    ca_dist = norm(ca_ss_xx[0:3] - ca_ast_xx[0:3])
                    drs_reach = ss_max_acc_i/2 * (ca_t - t0)**2

                    # If reachable append to close approach list
                    if ca_dist < drs_reach:
                        cas_i.append((ca_ast_id, ca_t))

                print("Remaining number of close approache indexes: {}".format(np.shape(cas_i)[0]), file=log_file)

                # Sort close approaches by time
                cas_dtype = [('idx', int), ('et', float)]
                cas_arr_i = np.array(cas_i, dtype=cas_dtype)
                cas_arr_i = np.sort(cas_arr_i, order='et')

                # Find solar sailing transfers
                N_feas_trans = 0
                candidate_len_before = len(candidate_list)

                for ca in cas_arr_i:
                    ca_idx = ca['idx']
                    ca_et = ca['et']

                    et_margin = 5 * pk.DAY2SEC  # 5 days of margin for this filter

                    # Rule out all possible solutions that are after tf
                    if ca_et > tf + et_margin:
                        print("Optimization stop: Close approaches after final time", file=log_file)
                        break  # Candidates are sorted by time so none after will be inside

                    if len(candidate_list) >= beam_width and candidate_list[beam_width-1][1] + et_margin < ca_et:
                        print("Optimization stop: Solutions in candidate list likely better than current CAs", file=log_file)
                        break  # Candidates are sorted by time so none after will be inside

                    xxf_fun = lambda t: neas.get_state(idx=ca_idx, et=(t0i + t))[0:3]
                    sol_tofi, sol_uusi, sol_fi = ss.target(xx0i, xxf_fun, tof_max=(t1i-t0i), N=N_control,
                                                           verbose=verbose_optim)

                    # If optimization converged, i.e. solar sailing transfer possible
                    if sol_fi < 1e-4:
                        # Also tried reoptimizing here but time gain was never more than 2 days and just took too long

                        N_feas_trans += 1
                        # If within time frame, save result in candidate list
                        t2i = t0i + sol_tofi
                        if t2i < tf:
                            parent_id_i = len(search_tree) + open_idx  # Computes where parent will be inserted in search_tree
                            node_level_i = tree_level_i + 1

                            # Propagate solution to obtain final state
                            sol_xxsi, _ = ss.propagate(xx0i, sol_tofi, sol_uusi)
                            xx2i = sol_xxsi[-1]

                            candidate_node_i = (xx2i.copy(), t2i, sol_uusi.copy(), ca_idx.copy(), parent_id_i, node_level_i)
                            candidate_list.append(candidate_node_i)
                            print("Solution added to candidate list", file=log_file)
                        else:
                            print("Solution found but not added because after final time", file=log_file)
                    else:
                        print("No solution found", file=log_file)
                
                    # Sort for next 
                    candidate_list.sort(key=lambda node: node[1])

                print("Feasible transfers found: {}".format(N_feas_trans), file=log_file)
                print("Transfers added to candidate list: {}".format(len(candidate_list) - candidate_len_before), file=log_file)

            # Update all the lists

            # Append open list to search tree
            search_tree += open_list

            # Sort candidate list by time
            candidate_list.sort(key=lambda node: node[1])
            if len(candidate_list) > beam_width:
                candidate_list = candidate_list[:beam_width]
            
            # Make candidate list new open list
            open_list = candidate_list.copy()
            candidate_list = []

            search_tree_dtypes = [('state', '6f8'), ('et', float), ('control', '({},2)f8'.format(N_control)),
                                ('ast_id', int), ('parent_id', int), ('level', int)]
            search_tree_np = np.array(search_tree, dtype=search_tree_dtypes)
            open_list_np = np.array(open_list, dtype=search_tree_dtypes)
            np.save("search_tree.npy", search_tree_np)
            np.save("open_list.npy", open_list_np)

            print(file=log_file)  # Just adds a line break

        print("Beam search finished. Total computation time: {} s".format(time.time() - code_start_t), file=log_file)
        # Save results
        print("Saving results ...", file=log_file)

        search_tree_dtypes = [('state', '6f8'), ('et', float), ('control', '({},2)f8'.format(N_control)),
                            ('ast_id', int), ('parent_id', int), ('level', int)]
        search_tree_np = np.array(search_tree, dtype=search_tree_dtypes)
        np.save("search_tree.npy", search_tree_np)

        print("Saving successful", file=log_file)

    pass