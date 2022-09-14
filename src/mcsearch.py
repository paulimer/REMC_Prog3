#!/usr/bin/env python3
import numpy as np

from src.conformation import Conformation

BOLTZMANN = 0.0019872
rng = np.random.default_rng()

def mc_search(
        current_conformation=Conformation("AA"), nb_steps=1, search_neigh="no_pull", temp = 200
):
    """
    Performs the Monte Carlo search from a given conformation

    Parameters
    ----------
    current_conformation: Conformation
    The conformation at the start of the search
    nb_steps: int
    the number of steps to perform
    search_neigh: str
    the type of search to perform
    temp: int
    the temperature of the search

    Returns
    -------
    the Conformation at the end of the search
    """
    for _ in range(nb_steps):
        residue = rng.integers(low=0, high=len(current_conformation.sequence)-1)
        list_moves = current_conformation.get_possible_moves(residue, search_neigh=search_neigh)
        if len(list_moves) == 0:
            continue
        elif len(list_moves) == 1:
            chosen_move = list_moves[0]
        else:
            chosen_move = list_moves[int(rng.integers(0, len(list_moves) - 1))]
        # print(chosen_move)
        # print(chosen_move.conf)
        energy_delta = chosen_move.conf.energy - current_conformation.energy
        if energy_delta <= 0:
            current_conformation = chosen_move.conf
            current_conformation.evaluate_energy()
        else:
            threshold = rng.random()
            if threshold < np.exp(-energy_delta/(temp*BOLTZMANN)):
                current_conformation = chosen_move.conf
                current_conformation.evaluate_energy()
    return current_conformation
