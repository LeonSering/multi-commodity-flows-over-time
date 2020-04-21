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
        self.commodityInflow = dict()
        self.commodityOutflow = dict()
        self.bIn = {edge:OrderedDict() for edge in self.network.edges} # b^+_e
        self.bOut = {edge:OrderedDict() for edge in self.network.edges} # b^-_e
        self.spillback = {node:OrderedDict() for node in self.network.nodes}   # c_v

    def add_commodity(self, path, startTime, endTime, rate):
        """Adds commodity to the dictionaries"""
        self.pathCommodityDict[path] = (startTime, endTime, rate)
        self.timeCommodityDict[(startTime, endTime)] = (path, rate)

    def edge_on_path(self, path, edge):
        """Check whether edge lies on path"""
        v, w = edge
        for j in range(len(path)-1):
            if (path[j] == v and path[j+1] == w):
                return True
        return False

    def validate_input(self):
        """Checks validity of input and returns error message if necessary"""

        def get_error_message(errorCode):
            errorDescription = {
                1: "Transittime must be greater than zero for all edges.",
                2: "Storage of initial edges (of comodities) must be infinite.",
                3: "Incapacities of initial edges (of comodities) have to be greater than the maximum cumulative rate.",
            }

            return errorDescription[errorCode]

        for edge in self.network.edges:
            v, w = edge
            if self.network[v][w]['transitTime'] <= 0:
                return get_error_message(1)

        mergedFlows = {}

        for path in self.pathCommodityDict:
            v, w = path[0], path[1]
            if self.network[v][w]['storage'] < float('inf'):
                return get_error_message(2)
            mergedFlows[(v, w)] = []

        sortedIntervalDict = OrderedDict(sorted(self.timeCommodityDict.items()))
        for interval in sortedIntervalDict:
            t_0, t_1 = interval
            path, rate = sortedIntervalDict[interval]
            v, w = path[0], path[1]
            e = (v, w)
            if len(mergedFlows[e]) == 0:
                mergedFlows[e].append((t_0, t_1, rate))
            else:
                staticL = list(mergedFlows[e])
                idx = 0
                idxShift = 0
                while idx < len(staticL):
                    t_l, t_u, r = staticL[idx]
                    if t_0 == t_u:  # Edge case
                        idx += 1
                    elif t_l <= t_0 < t_u and t_1 <= t_u:
                        # (t_0, t_1) completely contained in previous interval -> easy
                        lowSplit, highSplit = (t_l < t_0), (t_1 < t_u)
                        newL = []
                        if lowSplit:
                            newL.append((t_l, t_0, r))
                        newL.append((t_0, t_1, r + rate))
                        if highSplit:
                            newL.append((t_1, t_u, r))
                        mergedFlows[e][idx + idxShift:idx + idxShift + 1] = newL
                        idxShift += len(newL)-1
                        t_0 = t_1
                        break   # Nothing else to do here
                    elif t_l <= t_0 < t_u and t_1 > t_u:
                        lowSplit = (t_l < t_0)
                        newL = []
                        if lowSplit:
                            newL.append((t_l, t_0, r))
                        newL.append((t_0, t_u, r + rate))
                        mergedFlows[e][idx + idxShift:idx + idxShift + 1] = newL
                        idxShift += len(newL) - 1
                        t_0 = t_u   # Adjust interval for next iterations
                    else:
                        idx += 1
                if t_0 < t_1:
                    # Add to last case
                    mergedFlows[e].append((t_0, t_1, rate))
        for e in mergedFlows:
            v, w = e
            m = max(mergedFlows[e], key=lambda entry: entry[2])[2]
            if m > self.network[v][w]['inCapacity']:
                return get_error_message(3)

        return 0

    def init_values(self, t):
        """Init f, c, b up to initialTime t"""
        minInf = -float('inf')
        maxInf = float('inf')
        for v in self.network.nodes:
            self.spillback[v][(minInf, t)] = 1
        for e in self.network.edges:
            v, w = e
            tau = self.network[v][w]['transitTime']
            self.inflow[e][(minInf, t)] = 0
            self.outflow[e][(minInf, t)] = 0
            self.bIn[e][(minInf, t + tau)] = self.network[v][w]['inCapacity'] # not full
            self.bOut[e][(minInf, t + tau)] = 0

        for path in self.pathCommodityDict:
            self.commodityInflow[path] = {edge:OrderedDict() for edge in self.network.edges}
            self.commodityOutflow[path] = {edge:OrderedDict() for edge in self.network.edges}
            for e in self.network.edges:
                if self.edge_on_path(path, e):
                    self.commodityInflow[path][e][(minInf, t)] = 0
                    self.commodityOutflow[path][e][(minInf, t)] = 0
                else:
                    # Edge never used by this commodity
                    self.commodityInflow[path][e][(minInf, maxInf)] = 0
                    self.commodityOutflow[path][e][(minInf, maxInf)] = 0

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

        # LOOP
        # Access first element of heap
        theta, v = heapq.heappop(self.priority)

        # STEP 1: Compute push rates into node

        # TODO: DONT FORGET TO USE TIE!!

