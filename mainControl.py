# ===========================================================================
# Author:       Max ZIMMER
# Project:      multi-commodity-flows-over-time 2020
# File:         mainControl.py
# Description:  Script to be called from Terminal
# ===========================================================================
from multiFlowClass import MultiFlow
from argparse import ArgumentParser
import pickle
import sys
import os
import timeit


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

    rc = mf.validate_input()
    if rc != 0:
        # Return code is error message
        sys.exit(rc)

    return mf


# Parser configuration
parser = ArgumentParser(description='Compute multi-commodity Flows over Time.')
parser.add_argument("networkFile", help="Network in .cg-Format.")
parser.add_argument("inflowFile", help="Inflow .txt File containing the commodities.")

args = parser.parse_args()

# Positional arguments
networkFile = args.networkFile
inflowFile = args.inflowFile

# Main
mf = read_files(networkFile, inflowFile)
print("Starting computation.\n")
mf.compute()
baseName = os.path.basename(networkFile)
baseName = baseName[:baseName.rfind(".")]
print("Generating output.\n")
startTime = timeit.default_timer()
mf.generate_output(os.path.dirname(os.path.abspath(networkFile)), baseName)
endTime = timeit.default_timer()
print("Finished generating output. Total elapsed: {0:.2f} seconds".format(endTime - startTime))
