#!/usr/bin/env python3

import numpy as np
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
