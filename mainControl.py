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
    # TODO: Check requirements here (e.g. is \tau_e > 0 f.a. e \in E?)
    mf = MultiFlow(network)
    with open(inflowFile, 'r') as fRead:
        firstLine = True
        for line in fRead:
            if firstLine:
                firstLine = False
            else:
                line = line.strip()
                rate, interval, path = line.split()
                startTime, endTime = interval.split(",")
                path = tuple(path.split(","))
                mf.add_commodity(path, float(startTime), float(endTime), float(rate))
    return mf



# Parser configuration
parser = ArgumentParser(description='Compute multi-commodity Flows over Time.')
parser.add_argument("networkFile", help="Spillback network in .cg-Format.")
parser.add_argument("inflowFile", help="Inflow .txt File containing the commodities.")

args = parser.parse_args()

# Positional arguments
networkFile = args.networkFile
inflowFile = args.inflowFile

# Main
mf = read_files(networkFile, inflowFile)
mf.compute()
