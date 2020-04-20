# ===========================================================================
# Author:       Max ZIMMER
# Project:      multi-commodity-flows-over-time 2020
# File:         multiFlowClass.py
# Description:  Class managing the multi commodity flow
# ===========================================================================

import heapq
from collections import OrderedDict

class MultiFlow:
    """Manages the multi commodity flow"""

    def __init__(self, network):
        """
        :param network: Networkx Digraph instance
        """
        self.network = network.copy()
        self.pathCommodityDict = dict() # Format: path:(startTime, endTime, rate) //type(path) == tuple
        self.timeCommodityDict = dict() # Format: (startTime, endTime): (path, rate)
        self.priority = None
        self.inflow = {edge:OrderedDict() for edge in self.network.edges}  # f^+_e
        self.outflow = {edge:OrderedDict() for edge in self.network.edges} # f^-_e
        self.bIn = {edge:OrderedDict() for edge in self.network.edges} # b^+_e
        self.bOut = {edge:OrderedDict() for edge in self.network.edges} # b^-_e
        self.spillback = {node:OrderedDict() for node in self.network.nodes}   # c_v

    def add_commodity(self, path, startTime, endTime, rate):
        """Adds commodity to the dictionaries"""
        self.pathCommodityDict[path] = (startTime, endTime, rate)
        self.timeCommodityDict[(startTime, endTime)] = (path, rate)

    def init_values(self, t):
        """Init f, c, b up to initialTime t"""
        minInf = -float('inf')
        for v in self.network.nodes:
            self.spillback[v][(minInf, t)] = 1
        for e in self.network.edges:
            v, w = e
            tau = self.network[v][w]['transitTime']
            self.inflow[e][(minInf, t)] = 0
            self.outflow[e][(minInf, t)] = 0
            self.bIn[e][(minInf, t + tau)] = self.network[v][w]['inCapacity'] # not full
            self.bOut[e][(minInf, t + tau)] = 0

    def compute(self):
        """The priority heap maintained works as follows:
            - entries of the form (time, node)
            - sorted by min(time) // the lower time, the higher priority
            - if (t, v) is in the heap, that means we know outflow rates of all incoming edges into v and inflow rates
            of all outgoing edges of v up to time t.
        """

        # Get the minimal start time, this is the time up to which we know everything, as nothing happens before that
        initialTime = min(self.timeCommodityDict.keys())[0]

        # Init priority heap (No need to heapify, as this is sorted)
        self.priority = [(initialTime, node) for node in self.network.nodes]

        # Init f+, f-, c_v, b+_e, b-_e
        self.init_values(initialTime)

        # Access first element of heap
        theta, v = heapq.heappop(self.priority)

