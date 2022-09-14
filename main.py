#!/usr/bin/env python3

import argparse
import sys

from src.conformation import Conformation
from src.remc import remc
from src.sequence_ff import get_sequence_from_file

if __name__ == "__main__":
        parser = argparse.ArgumentParser(description="Runs a Replica Exchange Monte Carlo algorithm to fold a protein, following the HP model")
        sequence_entry = parser.add_mutually_exclusive_group(required=True)
        sequence_entry.add_argument(
                "--file",
                nargs="?",
                help="Read the protein sequence from this file"
        )
        sequence_entry.add_argument(
                "--sequence",
                nargs="?",
                help="The sequence of the protein as a string"
        )
        parser.add_argument(
                "nb_replicas",
                nargs="?",
                type=int,
                help="The number of replicas (default 5)",
                default=5
        )
        parser.add_argument(
                "steps",
                nargs="?",
                type=int,
                help="The maximal number of replica exchange (default 20)",
                default=20
        )
        parser.add_argument(
                "local_steps",
                nargs="?",
                type=int,
                help="The number of local steps before an attempt to exchange labels (default 50)",
                default=50
        )
        parser.add_argument(
                "t_min",
                nargs="?",
                type=int,
                help="The minimum temperature of the replicas(default 160)",
                default=160
        )
        parser.add_argument(
                "t_max",
                nargs="?",
                type=int,
                help="The maximum temperature of the replicas (default 220)",
                default=220
        )
        parser.add_argument(
                "--energy",
                nargs="?",
                type=int,
                help="The energy to achieve"
        )
        args = parser.parse_args()

        if args.file:
                try:
                        sequence = get_sequence_from_file(args.file)
                except FileNotFoundError:
                        sys.exit(f"File {args.file} not found")
        else:
                sequence = args.sequence

        start_conf = Conformation(sequence=sequence)
        print("The starting conformation : ")
        print(start_conf)
        if args.energy:
                end_conf = remc(
                        start_conformation=start_conf,
                        nb_replica=5,
                        local_steps=args.local_steps,
                        step_limit=args.steps,
                        t_min=args.t_min,
                        t_max=args.t_max,
                        search_neigh="no_pull",
                        optimal_energy=args.energy
                )
        else:
                end_conf = remc(
                        start_conformation=start_conf,
                        nb_replica=5,
                        local_steps=args.local_steps,
                        step_limit=args.steps,
                        t_min=args.t_min,
                        t_max=args.t_max,
                        search_neigh="no_pull"
                )
        print("The final conformation : ")
        print(end_conf)



        # conf = Conformation(sequence="MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAKSVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHERCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNSSCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACPGRDRRTEEENLRKKGEPHHELPPGSTKRALPN NTSSSPQPKKKPLDGEYFTLQIRGRERFEMFRELNEALELKDAQAGKEPGGSRAHSSHLKSKKGQSTSRHKKLMFKTEGPDSD")
        # conf_start = Conformation(sequence="LAED"*4)
        # print(conf_start)
        # conf_opt = mc_search(conf_start, nb_steps=20)
        # conf_opt = remc(conf_start, nb_replica=5, local_steps=50, step_limit=50, optimal_energy=-30,t_min=160, t_max=220)
        # print("Final")
        # print(conf_opt)
