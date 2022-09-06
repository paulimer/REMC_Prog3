#!/usr/bin/env python3

import copy
import numpy as np
import argparse

class amino_acid_2d:
    # _amino_number = 0
    def __init__(self, x=0, y=0, hp_type = "H"):
        self.x = x
        self.y = y
        self.hp_type = hp_type
        # self._amino_number += 1


    def __str__(self):
        return(f"Cet acide aminé numéro {self._amino_number} "
               f"de type {self.hp_type} est à la position "
               f"x = {self.x}, y = {self.y}")

class conformation:
    def __init__(self, amino_list = [], energy = 0):
        self.amino_list = copy.deepcopy(amino_list)
        self.energy = energy

    # def __getslice__() <- peut-être utile à un moment
