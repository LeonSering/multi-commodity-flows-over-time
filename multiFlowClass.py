# ===========================================================================
# Author:       Max ZIMMER
# Project:      multi-commodity-flows-over-time 2020
# File:         multiFlowClass.py
# Description:  Class managing the multi commodity flow
# ===========================================================================

from heapq import *

class MultiFlow:
    """Manages the multi commodity flow"""

    def __init__(self, network):
        """
        :param network: Networkx Digraph instance
        """
        self.network = network.copy()
        self.pathCommodityDict = dict() # Format: path:(startTime, endTime, rate) //type(path) == tuple
        self.timeCommodityDict = dict() # Format: (startTime, endTime): (path, rate)

    def add_commodity(self, path, startTime, endTime, rate):
        """Adds commodity to the dictionaries"""
        self.pathCommodityDict[path] = (startTime, endTime, rate)
        self.timeCommodityDict[(startTime, endTime)] = (path, rate)


    def compute(self):
        # Get the minimal start time, this is the time up to which we know everything
        initialTime = min(self.timeCommodityDict.keys())[0]
        # Init priority heap
        self.priority = heapify([(initialTime, node) for node in self.network.nodes])

