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
            self.evaluate_energy()

    def create_amino_list(self):
        self.amino_list = np.ndarray(len(self.sequence), dtype=AminoAcid)
        for i, aa in enumerate(self.sequence):
            self.amino_list[i] = AminoAcid(one_letter_aa=aa, index=i)

    def get_next_position(self, prev_aa_num):
        allowed_pos = self.get_free_pos((self.amino_list[prev_aa_num].position))
        if len(allowed_pos) == 0:
            raise Exception
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
            except:
                # restart from beginning, this is not a limiting step
                print("No more room for sequence")
                self.assign_positions()
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

    def __str__(self):
        res = f"Cette conformation possède {len(self.sequence)} acides aminés\n"
        res += f"Cette conformation a pour énergie {self.energy}\n"
        res += " "
        res += "-" * len(self.lattice)
        res += "\n"
        for i in range(len(self.lattice)):
            res += "|"
            for j in range(len(self.lattice)):
                if self.lattice[i, j] is None:
                    res += (" ")
                elif self.lattice[i, j].hp_type == "H":
                    res += "H"
                elif self.lattice[i, j].hp_type == "P":
                    res += "P"
            res += "|\n"
        res += " "
        res += "-" * len(self.lattice)
        return res

    def evaluate_energy(self):
        self.energy = 0
        for i in range(len(self.amino_list)):
            for j in range(i + 1, len(self.amino_list)):
                if self.amino_list[i].hp_type == "H" and self.amino_list[j].hp_type == "H" \
                   and np.linalg.norm(self.amino_list[i].position - self.amino_list[j].position) <= 1.1:
                    self.energy += 1

    def get_possible_moves(self, aa_number):
        """ shows which moves are available"""
        res = []
        if aa_number == 0 or aa_number == len(self.sequence):
            if len(self.get_free_pos()) != 0:
                res += ["end"]
        else:
            # corner move
            prev_next_diff = self.amino_list[aa_number - 1].position - \
                      self.amino_list[aa_number + 1].position
            if sum(np.abs(prev_next_diff) == np.array((1, 1))) == 2:
                if self.amino_list[aa_number].position[0] == self.amino_list[aa_number -1].position[0]:
                    if self.lattice[self.amino_list[aa_number + 1].position[0], self.amino_list[aa_number - 1].position[1]]:
                        res += ["corner"]
                else:
                    if self.lattice[self.amino_list[aa_number - 1].position[0], self.amino_list[aa_number + 1].position[1]]:
                        res += ["corner"]

            # crankshaft move
            if np.linalg.norm(self.amino_list[aa_number - 2].position - self.amino_list[aa_number + 1].position) == 1:
                # need-to-be empty positions:
                in_front_cur = self.amino_list[aa_number].position + 2*(self.amino_list[aa_number - 1].position -  self.amino_list[aa_number].position)
                in_front_prev = self.amino_list[aa_number - 1].position + 2*(self.amino_list[aa_number - 2].position -  self.amino_list[aa_number - 1].position)
                if not (self.lattice[in_front_cur[0], in_front_cur[1]]) and not (self.lattice[in_front_prev[0], self.lattice[in_front_prev[1]]]):
                    res += ["crank"]
            if np.linalg.norm(self.amino_list[aa_number - 1].position - self.amino_list[aa_number + 2].position) == 1:
                # need-to-be empty positions:
                in_front_cur = self.amino_list[aa_number].position + 2*(self.amino_list[aa_number - 1].position -  self.amino_list[aa_number].position)
                in_front_next = self.amino_list[aa_number + 1].position + 2*(self.amino_list[aa_number + 2].position -  self.amino_list[aa_number + 1].position)
                if not (self.lattice[in_front_cur[0], in_front_cur[1]]) and not (self.lattice[in_front_prev[0], self.lattice[in_front_prev[1]]]):
                    res += ["crank"]
        return res

    def end_move(self, side):
        if side != 0 and side != len(self.sequence) - 1:
            raise ValueError
        new_position_list = self.get_free_pos(self.amino_list[side].position)
        new_position = new_position_list[rng.integers(len(new_position_list), size=1)]
        self.lattice[new_position[0], new_position[1]] = self.amino_list[side]
        self.lattice[self.amino_list[side].position[0], self.amino_list[side].position[1]] = None
        self.amino_list[side].position = np.array(new_position)


    def corner_move(self):


    # def pull_move(self):

if __name__ == "__main__":
    conf = Conformation(sequence="AREAAR")
    print(conf)
