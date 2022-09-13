#!/usr/bin/env python3

import argparse
import numpy as np

from conformation import Conformation
from mcsearch import mc_search
from remc import remc

if __name__ == "__main__":
        # parser = argparse.ArgumentParser(description="Runs a Replica Exchange Monte Carlo algorithm to fold a protein, following the HP model")
        # parser.add_argument(
        #     "--file",
        #     nargs="?",
        #     help="Read the protein sequence from this file"
        # )
        # parser.add_argument(
        #     "--sequence",
        #     nargs="?",
        #     help="The sequence of the protein as a string"
        # )
        # parser.add_argument(
        #     "steps",
        #     nargs="?",
        #     type=int,
        #     help="The number of local steps before an attempt to exchange labels"
        # )
        # conf = Conformation(sequence="MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAKSVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHERCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNSSCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACPGRDRRTEEENLRKKGEPHHELPPGSTKRALPN NTSSSPQPKKKPLDGEYFTLQIRGRERFEMFRELNEALELKDAQAGKEPGGSRAHSSHLKSKKGQSTSRHKKLMFKTEGPDSD")
        conf_start = Conformation(sequence="LAED"*4)
        print(conf_start)
        conf_opt = mc_search(conf_start, nb_steps=20)
        # conf_opt = remc(conf_start, nb_replica=5, local_steps=50, step_limit=50, optimal_energy=-30,t_min=160, t_max=220)
        # print("Final")
        print(conf_opt)
