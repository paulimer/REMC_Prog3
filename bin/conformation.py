#!/usr/bin/env python3

import argparse
import copy
import numpy as np

from moves import Move
from amino import AminoAcid

rng = np.random.default_rng()



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

    def __str__(self):
        res = f"Cette conformation possède {len(self.sequence)} acides aminés\n"
        res += f"Cette conformation a pour énergie {self.energy}\n"
        res += " "
        min_x = max(self.get_extr_coor(along="x", min=True) - 3, 0)
        max_x = min(self.get_extr_coor(along="x", min=False) + 3, 2*len(self.sequence))
        min_y = max(self.get_extr_coor(along="y", min=True) - 3, 0)
        max_y = min(self.get_extr_coor(along="y", min=False) + 3, 2*len(self.sequence))
        res += "-" * (max_y - min_y)
        res += "\n"
        for i in range(min_x, max_x):
            res += "|"
            for j in range(min_y, max_y):
                if self.lattice[i, j] is None:
                    res += (" ")
                elif self.lattice[i, j].hp_type == "H":
                    res += "H"
                elif self.lattice[i, j].hp_type == "P":
                    res += "P"
            res += "|\n"
        res += " "
        res += "-" * (max_y - min_y)
        return res

    def get_extr_coor(self, along="x", min=True):
        if min:
            extr = 2*self.size
            if along == "x":
                dim = 0
            else:
                dim = 1
            for amino in self.amino_list:
                if amino.position[dim] < extr:
                    extr = amino.position[dim]
            return extr
        else:
            extr = 0
            if along == "x":
                dim = 0
            else:
                dim = 1
            for amino in self.amino_list:
                if amino.position[dim] > extr:
                    extr = amino.position[dim]
            return extr

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

    def evaluate_energy(self):
        self.energy = 0
        for i in range(len(self.amino_list)):
            for j in range(i + 1, len(self.amino_list)):
                if self.amino_list[i].hp_type == "H" and self.amino_list[j].hp_type == "H" \
                   and np.linalg.norm(self.amino_list[i].position - self.amino_list[j].position) <= 1.1:
                    self.energy += 1

    def get_possible_moves(self, aa_number=0):
        """ shows which moves are available"""

        res = []
        if aa_number == 0 or aa_number == len(self.sequence) - 1:
            res += self.get_end_moves(aa_number=aa_number)
        else:
            # corner move
            prev_next_diff = self.amino_list[aa_number - 1].position - \
                      self.amino_list[aa_number + 1].position
            if sum(np.abs(prev_next_diff) == np.array((1, 1))) == 2:
                res += self.get_corner_move(aa_number=aa_number)

            # crankshaft move
            res += self.get_crankshaft_move(aa_number=aa_number)
        return res

    def get_end_moves(self, aa_number=0):
        res = []
        for pos in self.get_free_pos(self.amino_list[aa_number].position):
            next_move = Move(move_type="end",conf=self, number=aa_number,new_position=pos)
            next_move.end_move()
            res += [next_move]
        return res

    def get_corner_move(self, aa_number=0):
        res = []
        if self.amino_list[aa_number].position[0] == self.amino_list[aa_number -1].position[0]:
            if not self.lattice[self.amino_list[aa_number + 1].position[0], self.amino_list[aa_number - 1].position[1]]:
                next_move = Move(move_type="corner", conf=self, number=aa_number, new_position=\
                                 np.array((self.amino_list[aa_number + 1].position[0], self.amino_list[aa_number - 1].position[1])))
                next_move.corner_move()
                res += [next_move]
        else:
            if not self.lattice[self.amino_list[aa_number - 1].position[0], self.amino_list[aa_number + 1].position[1]]:
                next_move = Move(move_type="corner", conf=self, number=aa_number, new_position=\
                                 np.array((self.amino_list[aa_number - 1].position[0], self.amino_list[aa_number + 1].position[1])))
                next_move.corner_move()
                res += [next_move]
        return res

    def get_crankshaft_move(self, aa_number=0):
        res = []
        if aa_number >= 2 and aa_number <= len(self.sequence) - 2 \
           and np.linalg.norm(self.amino_list[aa_number - 2].position - self.amino_list[aa_number + 1].position) == 1:
            # need-to-be empty positions:
            in_front_cur = self.amino_list[aa_number].position + 2*(self.amino_list[aa_number - 1].position -  self.amino_list[aa_number].position)
            in_front_prev = self.amino_list[aa_number - 1].position + 2*(self.amino_list[aa_number - 2].position -  self.amino_list[aa_number - 1].position)
            if not (self.lattice[in_front_cur[0], in_front_cur[1]]) and not (self.lattice[in_front_prev[0], in_front_prev[1]]):
                next_move = Move(move_type="crank", conf=self, number=aa_number, new_position=in_front_cur)
                next_move.crankshaft_move(neightbour_pos=in_front_prev, front=False)
                res += [next_move]
        if aa_number >= 1 and aa_number <= len(self.sequence) - 3 \
           and np.linalg.norm(self.amino_list[aa_number - 1].position - self.amino_list[aa_number + 2].position) == 1:
            # need-to-be empty positions:
            in_front_cur = self.amino_list[aa_number].position + 2*(self.amino_list[aa_number - 1].position -  self.amino_list[aa_number].position)
            in_front_next = self.amino_list[aa_number + 1].position + 2*(self.amino_list[aa_number + 2].position -  self.amino_list[aa_number + 1].position)
            if not (self.lattice[in_front_cur[0], in_front_cur[1]]) and not (self.lattice[in_front_next[0], in_front_next[1]]):
                next_move = Move(move_type="crank", conf=self, number=aa_number, new_position=in_front_cur)
                next_move.crankshaft_move(neightbour_pos=in_front_next, front=True)
                res += [next_move]

        return res



if __name__ == "__main__":
    conf = Conformation(sequence="AREAAR")
    print(conf)
