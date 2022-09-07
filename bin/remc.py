#!/usr/bin/env python3

import argparse
import copy
import numpy as np

rng = np.random.default_rng()

HYDROPHOBIC_AMINOS = ["A", "V", "I", "L", "M", "F", "Y", "W"]


class AminoAcid:
    def __init__(self, position=np.array((0, 0)), one_letter_aa="A", index=0):
        self.position = position
        self.hp_type = self.get_type(one_letter_aa)
        self.index = index

    @staticmethod
    def get_type(one_letter_aa):
        if one_letter_aa in HYDROPHOBIC_AMINOS:
            return "H"
        return "P"


    def __str__(self):
        return(f"Cet acide aminé numéro {self.index} "
               f"de type {self.hp_type} a pour position "
               f"x = {self.position[0]}, y = {self.position[1]}.")


class Conformation:
    def __init__(self, sequence="", amino_list = None, lattice=None, energy=None):
        self.sequence = sequence

        if amino_list:
            self.amino_list = copy.deepcopy(amino_list)
        else:
            self.create_amino_list()
        if lattice:
            self.lattice = copy.deepcopy(lattice)
        else:
            self.assign_positions()
        if energy:
            self.energy = energy
        else:
            pass
            # self.evaluate_energy()

    def create_amino_list(self):
        self.amino_list = np.ndarray(len(self.sequence), dtype=AminoAcid)
        for i, aa in enumerate(self.sequence):
            self.amino_list[i] = AminoAcid(one_letter_aa=aa, index=i)

    def get_next_position(self, prev_aa_num):
        allowed_pos = self.get_free_pos((self.amino_list[prev_aa_num].position))
        if len(allowed_pos) == 0:
            raise ValueError
        return allowed_pos[int(rng.integers(len(allowed_pos), size=1))]

    def assign_positions(self):
        self.lattice = np.ndarray(shape=(len(self.sequence)*2, len(self.sequence)*2), dtype=AminoAcid)
        # starting at 1, amino acid 0 is in (len(self.sequence) - 1, len(self.sequence) -1)
        self.lattice[len(self.sequence) - 1, len(self.sequence) - 1] = self.amino_list[0]
        self.amino_list[0].position = np.array((len(self.sequence) - 1, len(self.sequence) - 1))
        for i in range(1, len(self.sequence)):
            try:
                new_position = self.get_next_position(i-1)
                self.lattice[new_position[0], new_position[1]] = self.amino_list[i]
                self.amino_list[i].position = np.array(new_position)
            except ValueError:
                print("No more room for sequence")
                return

    def get_free_pos(self, start_pos):
        res = []
        movements = np.array([(-1, 0), (0, -1), (1, 0), (0, 1)])
        for move in movements:
            cur_pos = start_pos + move
            if (cur_pos[0] >= 0 and cur_pos[0] < len(self.sequence)*2) and \
               (cur_pos[1] >= 0 and cur_pos[1] < len(self.sequence)*2) and \
               self.lattice[cur_pos[0], cur_pos[1]] is None:
                res.append((cur_pos[0], cur_pos[1]))
        return res


if __name__ == "__main__":
    conf = Conformation(sequence="AREAAR")
