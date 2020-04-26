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

    def dictInSortAdd(self, OD, newValues):
        """
        Update OrderedDict OD by valList suitable for ADDING a flow over time
        :param OD: OrderedDict where keys are of the form (startTime, endTime) and values are flow rates s.t.
        for (s, t): r before (s', t'):r' we have s <= t <= s' <= t'
        :param newValues: List of values (s,t,r) that need to be added to previous dict, i.e. by adding rates at times
        """
        if type(newValues) is tuple:
            newValues = [newValues]
        newValues = sorted(newValues)
        keyValList = [(time[0], time[1], OD[time]) for time in OD]
        for t_0, t_1, rate in newValues:
            if len(keyValList) == 0:
                # Just add
                keyValList.append((t_0, t_1, rate))
            else:
                staticL = list(keyValList)
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
                        keyValList[idx + idxShift:idx + idxShift + 1] = newL
                        idxShift += len(newL)-1
                        t_0 = t_1
                        break   # Nothing else to do here
                    elif t_l <= t_0 < t_u < t_1:
                        lowSplit = (t_l < t_0)
                        newL = []
                        if lowSplit:
                            newL.append((t_l, t_0, r))
                        newL.append((t_0, t_u, r + rate))
                        keyValList[idx + idxShift:idx + idxShift + 1] = newL
                        idxShift += len(newL) - 1
                        t_0 = t_u   # Adjust interval for next iterations
                    else:
                        idx += 1
                if t_0 < t_1:
                    # Add to last case
                    keyValList.append((t_0, t_1, rate))

        OD.clear()
        OD.update(OrderedDict([((triplet[0], triplet[1]), triplet[2]) for triplet in keyValList]))

    def dictInSort(self, OD, newValues):
        """
        Update OrderedDict OD by valList suitable for inserting a flow over time
        :param OD: OrderedDict where keys are of the form (startTime, endTime) s.t.
        for (s, t): r before (s', t'):r' we have s <= t <= s' <= t'
        :param newValues: List of values (s,t,val) that need to be added to previous dict, i.e. by inserting val at times
        """
        if type(newValues) is tuple:
            newValues = [newValues]
        newValues = sorted(newValues)
        keyValList = [(time[0], time[1], OD[time]) for time in OD]
        for t_0, t_1, rate in newValues:
            if len(keyValList) == 0:
                # Just add
                keyValList.append((t_0, t_1, rate))
            else:
                staticL = list(keyValList)
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
                        keyValList[idx + idxShift:idx + idxShift + 1] = newL
                        idxShift += len(newL)-1
                        t_0 = t_1
                        break   # Nothing else to do here
                    elif t_l <= t_0 < t_u < t_1:
                        lowSplit = (t_l < t_0)
                        newL = []
                        if lowSplit:
                            newL.append((t_l, t_0, r))
                        newL.append((t_0, t_u, r + rate))
                        keyValList[idx + idxShift:idx + idxShift + 1] = newL
                        idxShift += len(newL) - 1
                        t_0 = t_u   # Adjust interval for next iterations
                    else:
                        idx += 1
                if t_0 < t_1:
                    # Add to last case
                    keyValList.append((t_0, t_1, rate))

        OD.clear()
        OD.update(OrderedDict([((triplet[0], triplet[1]), triplet[2]) for triplet in keyValList]))

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
                2: "Storage of initial edges (of commodities) have to be infinite.",
                3: "Incapacities of initial edges (of commodities) have to be infinite.",
            }

            return errorDescription[errorCode]

        for edge in self.network.edges:
            v, w = edge
            if self.network[v][w]['transitTime'] <= 0:
                return get_error_message(1)

        #mergedFlows = {}

        for path in self.pathCommodityDict:
            v, w = path[0], path[1]
            if self.network[v][w]['storage'] < float('inf'):
                return get_error_message(2)
            if self.network[v][w]['inCapacity'] < float('inf'):
                return get_error_message(3)
            #mergedFlows[(v, w)] = []
        """
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
                    elif t_l <= t_0 < t_u < t_1:
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
        """
        return 0

    def init_values(self, t, startingEdges, isolatedNodes):
        """Init f, c, b up to initialTime t using the fact that we have known startingEdges and isolatedNodes"""
        minInf = -float('inf')
        maxInf = float('inf')
        for v in self.network.nodes:
            if v not in isolatedNodes:
                self.spillback[v][(minInf, t)] = 1
            else:
                # Nodes with no incoming edges have spillback factor always equal to 1
                self.spillback[v][(minInf, maxInf)] = 1

        # Init dictionaries
        for path in self.pathCommodityDict:
            self.commodityInflow[path] = {edge:OrderedDict() for edge in self.network.edges}
            self.commodityOutflow[path] = {edge:OrderedDict() for edge in self.network.edges}

        for e in self.network.edges:
            v, w = e
            tau = self.network[v][w]['transitTime']
            #self.inflow[e][(minInf, t)] = 0
            #self.outflow[e][(minInf, t)] = 0
            if e not in startingEdges:
                self.bIn[e][(minInf, t + tau)] = self.network[v][w]['inCapacity'] # Not full in this time frame
            else:
                # Starting edges have infinite storage, are hence never full
                self.bIn[e][(minInf, maxInf)] = self.network[v][w]['inCapacity']

            self.bOut[e][(minInf, t + tau)] = 0

            for path in self.pathCommodityDict:
                if v == path[0] and w == path[1]:
                    # This is the initial edge, so we can safely map the flow
                    startTime, endTime, rate = self.pathCommodityDict[path]
                    self.commodityInflow[path][e][(minInf, startTime)] = 0
                    self.commodityInflow[path][e][(startTime, endTime)] = rate
                    self.commodityInflow[path][e][(endTime, maxInf)] = 0

                    self.commodityOutflow[path][e][(minInf, t + tau)] = 0   # Here we can't know more than this
                elif self.edge_on_path(path, e):
                    # Edge is on path, but due to elif not the initial one
                    # Get the sum of transitTimes up to that edge
                    sTT, idx = 0, 0
                    while idx < len(path)-1:
                        if path[idx] == v:
                            break
                        else:
                            sTT += self.network[path[idx]][path[idx+1]]['transitTime']
                            idx += 1
                    self.commodityInflow[path][e][(minInf, t + sTT)] = 0
                    self.commodityOutflow[path][e][(minInf, t + sTT + tau)] = 0
                else:
                    # Edge never used by this commodity
                    self.commodityInflow[path][e][(minInf, maxInf)] = 0
                    self.commodityOutflow[path][e][(minInf, maxInf)] = 0

    def compute_alpha(self, theta, v, in_edges):
        """Computes alpha = min_{e in delta^-(v)} alpha_e such that commodity inflow constant"""
        alpha = float('inf')
        for e in in_edges:
            # Compute alpha to get extension size
            u, _ = e    # e = (u,v)
            tau = self.network[u][v]['transitTime']
            t = theta - tau
            print("time ", t)

            # Find maximal alpha such that commodity flows stay constant
            for path in self.pathCommodityDict:
                partAlpha = 0
                lastRate = None
                found = False
                for interval in self.commodityInflow[path][e]:
                    t_l, t_u = interval
                    if t_l <= t < t_u:
                        lastRate = self.commodityInflow[path][e][interval]
                        partAlpha += t_u - t
                        found = True
                    elif found:
                        if lastRate == self.commodityInflow[path][e][interval]:
                            # We can extend further
                            partAlpha += t_u - t_l
                        else:
                            break
                alpha = min(alpha, partAlpha)   # TODO: This works as extension only in the case of no spillback!
        return alpha

    def extend_bOut(self, theta, alpha, in_edges):
        thetaFixed, alphaFixed = theta, alpha
        for e in in_edges:
            u, v = e
            tau = self.network[u][v]['transitTime']
            theta = thetaFixed
            alpha = alphaFixed
            if theta < tau:
                # Adjust extension interval, as bOut(t) = 0 for time < tau
                # Note: No need to update bOut, we have that information already (see initialization)
                alpha = alpha - (tau-theta)
                theta = tau
            t = theta - tau
            print(theta, alpha)
            for path in self.pathCommodityDict:
                print(self.commodityInflow[path][e])
                print(self.commodityOutflow[path][e])



    def compute(self):
        """The priority heap maintained works as follows:
            - entries of the form (time, node)
            - sorted by min(time) // the lower time, the higher priority
            - if (t, v) is in the heap, that means we know outflow rates of all incoming edges into v and inflow rates
            of all outgoing edges of v up to time t.
        """

        # Get the minimal start time, this is the time up to which we know everything, as nothing happens before that
        initialTime = min(self.timeCommodityDict.keys())[0]
        startingEdges = [(path[0], path[1]) for path in self.pathCommodityDict]
        isoNodes = [v for v in self.network.nodes if self.network.in_degree(v) == 0]

        # Init f+, f-, c_v, b+_e, b-_e
        self.init_values(initialTime, startingEdges, isoNodes)

        # Init priority heap (No need to heapify, as this is sorted) (Note: isolated nodes not included)
        self.priority = [(initialTime, node) for node in self.network.nodes if node not in isoNodes]
        print(self.priority)

        # LOOP
        # Access first element of heap
        theta, v = heapq.heappop(self.priority)
        print("v: ", v, " theta: ", theta)
        # STEP 1: Compute alpha extension size
        in_edges = list(self.network.in_edges(v))
        alpha = self.compute_alpha(theta, v, in_edges)
        print(alpha)

        # STEP 2: Compute push rates into node
        self.extend_bOut(theta, alpha, in_edges)


        # We compute bOut of that edge  # TODO: Do we have to take care of the case theta < tau (i.e. bOut = 0)?



        # TODO: DONT FORGET TO USE TIE!!

