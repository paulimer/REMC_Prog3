#!/usr/bin/env python3
import numpy as np
rng = np.random.default_rng()

from conformation import Conformation
from mcsearch import mc_search, BOLTZMANN

def remc(
        start_conformation, nb_replica, local_steps, step_limit, optimal_energy, t_min, t_max, search_neigh="no_pull"
):
    current_energy_list = np.zeros(nb_replica)
    array_t = np.linspace(t_min, t_max)
    replica_list = np.ndarray(nb_replica, dtype=Conformation)
    for i in range(len(replica_list)):
        replica_list[i] = start_conformation
    offset = 0
    global_steps = 0
    while(min(current_energy_list) > optimal_energy and global_steps < step_limit):
        for i in range(len(replica_list)):
            replica_list[i] = mc_search(current_conformation=replica_list[i],
                                        nb_steps=local_steps,
                                        search_neigh=search_neigh,
                                        temp=array_t[i])
            if replica_list[i].energy < current_energy_list[i]:
                current_energy_list[i] = replica_list[i].energy
        repl = offset
        while repl < nb_replica - 1:
            nex_repl = repl + 1
            delta_ener = (1/array_t[nex_repl]*BOLTZMANN - 1/array_t[repl]*BOLTZMANN)*\
                (current_energy_list[repl] - current_energy_list[nex_repl])
            if delta_ener <= 0:
                tmp_repl = replica_list[repl]
                replica_list[repl] = replica_list[nex_repl]
                replica_list[nex_repl] = tmp_repl
            else:
                q = rng.random()
                if q <= np.exp(delta_ener):
                    tmp_repl = replica_list[repl]
                    replica_list[repl] = replica_list[nex_repl]
                    replica_list[nex_repl] = tmp_repl
            repl += 2
        offset = 1 - offset
        global_steps += 1

    min_energy_index = np.argmin(current_energy_list)
    return replica_list[min_energy_index]
