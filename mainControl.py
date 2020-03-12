# ===========================================================================
# Author:       Max ZIMMER
# Project:      multi-commodity-flows-over-time 2020
# File:         mainControl.py
# Description:  Script to be called from Terminal
# ===========================================================================
from multiFlowClass import MultiFlow
from argparse import ArgumentParser
import pickle


def load_graph(path):
    """
    Load graph instance from '.cg' file
    :param path: Path to load graph from
    :return: returns networkx graph instance
    """
    with open(path, 'rb') as f:
        network = pickle.load(f)
    return network

def read_files(networkFile, inflowFile):
    """
    Reads the files and initiates MultiFlow instance
    :param networkFile: networkx graph
    :param inflowFile: File containing commodities
    :return: MultiFlow object
    """
    network = load_graph(networkFile)
    multiFlow = MultiFlow(network)
    with open(inflowFile, 'r') as fRead:
        for line in fRead:
            line = line.strip()
            rate, startTime, endTime, path = line.split(" ")
            path = path.split(",")
            multiFlow.add_commodity(path, startTime, endTime, rate)

    return MultiFlow



# Parser configuration
parser = ArgumentParser(description='Compute multi-commodity Flows over Time.')
parser.add_argument("networkFile", help="Spillback network in .cg-Format.")
parser.add_argument("inflowFile", help="Inflow .txt File containing the commodities.")

args = parser.parse_args()

# Positional arguments
networkFile = args.networkFile
inflowFile = args.inflowFile

# Main
multiFlow = read_files(networkFile, inflowFile)
