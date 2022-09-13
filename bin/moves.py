#!/usr/bin/env python3

import copy
import numpy as np

class Conformation():
    def __init__(self, sequence="", amino_list = None, lattice=None, energy=None):
        pass
    pass

class Move():
    def __init__(self, move_type=None, conf=Conformation(sequence="AA"), number=0, new_position=np.array((0, 0))):
        self.move_type = move_type
        self.conf = copy.deepcopy(conf)
        self.new_position = new_position
        self.old_position = copy.deepcopy(self.conf.amino_list[number].position)
        self.number = number

    def __str__(self):
        return f"Move type : {self.move_type}, New position : {self.new_position}, Old position : {self.old_position}"

    def end_move(self):
        if self.move_type != "end":
            raise Exception("Error in move assignment")

        self.conf.lattice[self.new_position[0], self.new_position[1]] = self.conf.amino_list[self.number]
        self.conf.lattice[self.conf.amino_list[self.number].position[0], self.conf.amino_list[self.number].position[1]] = None
        self.conf.amino_list[self.number].position = np.array(self.new_position)

    def corner_move(self):
        if self.move_type != "corner":
            raise Exception("Error in move assignment")

        self.conf.lattice[self.new_position[0], self.new_position[1]] = self.conf.amino_list[self.number]
        self.conf.lattice[self.conf.amino_list[self.number].position[0], self.conf.amino_list[self.number].position[1]] = None
        self.conf.amino_list[self.number].position = np.array(self.new_position)

    def crankshaft_move(self, neightbour_pos=np.array((0, 0)), front=True):
        if self.move_type != "crank":
            raise Exception("Error in move assignment")
        #current aa
        self.conf.lattice[self.new_position[0], self.new_position[1]] = self.conf.amino_list[self.number]
        self.conf.lattice[self.conf.amino_list[self.number].position[0], self.conf.amino_list[self.number].position[1]] = None
        self.conf.amino_list[self.number].position = np.array(self.new_position)
        #neighbour aa
        if front:
            direction = 1
        else:
            direction = -1
        self.conf.lattice[neightbour_pos[0], neightbour_pos[1]] = self.conf.amino_list[self.number + direction]
        self.conf.lattice[self.conf.amino_list[self.number + direction].position[0], self.conf.amino_list[self.number + direction].position[1]] = None
        self.conf.amino_list[self.number + direction].position = np.array(self.new_position)
