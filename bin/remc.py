#!/usr/bin/env python3

import copy
import numpy as np
import argparse

HYDROPHOBIC_AMINOS = ["A", "V", "I", "L", "M", "F", "Y", "W"]

class amino_acid:
    def __init__(self, number=0, one_letter_aa = "A"):
        self.number = number
        self.hp_type = self.get_type(one_letter_aa)

    @staticmethod
    def get_type(one_letter_aa):
        if one_letter_aa in HYDROPHOBIC_AMINOS:
            return "H"
        return "P"


    def __str__(self):
        return(f"Cet acide aminé "
               f"de type {self.hp_type} a pour numéro "
               f"{self.number}")

class conformation:
    def __init__(self, sequence = [], amino_list = None, amino_positions = None, energy = None):
        self.sequence = sequence
        if amino_list:
            self.amino_list = copy.deepcopy(amino_list)
        else:
            self.create_amino_list()

        if amino_positions:
            self.amino_positions = amino_positions
        else:
            self.assign_positions()

        if energy:
            self.energy = energy
        else:
            self.evaluate_energy()

    def create_amino_list(self):
        res = []
        for i, one_letter_aa in self.sequence:
            res.append(amino_acid(i, one_letter_aa))
        return res

    def get_next_position(self, i):
        allowed_movement = []
        # à remplir avec les positions vides
        # puis tirer à la fin
        # Gérer l'exception de pas de mouvement possible
        # probablement un try catch
        if self.amino_positions[i, :] +


    def assign_positions(self):
        res = np.zeros((len(self.sequence), 2))
        # starting at 1, amino acid 0 is in (0, 0)
        for i in range(1, len(self.sequence)):
            res[:, i] = self.get_next_position(i-1)
