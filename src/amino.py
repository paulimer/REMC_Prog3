#!/usr/bin/env python3

import numpy as np
HYDROPHOBIC_AMINOS = ["A", "V", "I", "L", "M", "F", "Y", "W"]

class AminoAcid:
    """
    Class representing an amino acid

    Attributes
    ----------
    position: np.array (shape : 2)
    the position on the lattice of the AminoAcid
    hp_type: str
    the type of the AminoAcid according to the HP model
    index: int
    the index of the AminoAcid in its sequence

    Methods
    -------
    get_type(one_letter_aa)
    gets the HP type of the amino acid based on its residue one letter name
    """
    def __init__(self, position=np.array((0, 0)), one_letter_aa="A", index=0):
        self.position = position
        self.hp_type = self.get_type(one_letter_aa)
        self.index = index

    @staticmethod
    def get_type(one_letter_aa):
        """gets the HP type of the amino acid based on its residue one letter name"""
        if one_letter_aa in HYDROPHOBIC_AMINOS:
            return "H"
        return "P"

    def __str__(self):
        return(f"Cet acide aminé numéro {self.index} "
               f"de type {self.hp_type} a pour position "
               f"x = {self.position[0]}, y = {self.position[1]}.")
