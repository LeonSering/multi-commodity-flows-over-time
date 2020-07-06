# ===========================================================================
# Author:       Max ZIMMER
# Project:      multi-commodity-flows-over-time 2020
# File:         multiFlowClass.py
# Description:  Class managing the multi commodity flow
# ===========================================================================

import heapq
from collections import OrderedDict
from utilitiesClass import Utilities
import os
import timeit


class MultiFlow:
    """Manages the multi commodity flow"""

    def __init__(self, network):
        """
        :param network: Networkx Digraph instance
        """
        self.network = network.copy()
        self.pathCommodityDict = dict()  # Format: path:(startTime, endTime, rate) //type(path) == tuple
        self.timeCommodityDict = dict()  # Format: (startTime, endTime): (path, rate)
        self.priority = None
        self.commodityInflow = dict()
        self.commodityOutflow = dict()
        self.hasBeenActive = dict()
        self.bIn = {edge: OrderedDict() for edge in self.network.edges}  # b^+_e
        self.bOut = {edge: OrderedDict() for edge in self.network.edges}  # b^-_e

    def add_commodity(self, path, startTime, endTime, rate):
        """Adds commodity to the dictionaries"""
        self.pathCommodityDict[path] = (startTime, endTime, rate)
        self.timeCommodityDict[(startTime, endTime)] = (path, rate)

    def edge_on_path(self, path, edge):
        """Check whether edge lies on path"""
        v, w = edge
        for j in range(len(path) - 1):
            if path[j] == v and path[j + 1] == w:
                return True
        return False

    def ensure_general_model(self):
        """Ensures that the model is according to Skutella-Koch"""
        for e in self.network.edges:
            v, w = e
            self.network[v][w]['storage'] = float('inf')
            self.network[v][w]['inCapacity'] = float('inf')

    def validate_input(self):
        """Checks validity of input and returns error message if necessary"""
        self.ensure_general_model()

        def get_error_message(errorCode):
            errorDescription = {
                1: "Transittime must be greater than zero for all edges.",
            }

            return errorDescription[errorCode]

        for edge in self.network.edges:
            v, w = edge
            if self.network[v][w]['transitTime'] <= 0:
                return get_error_message(1)

        return 0

    def init_values(self, t, startingEdges):
        """Init f, c, b up to initialTime t using the fact that we have known startingEdges and isolatedNodes"""
        minInf = -float('inf')
        maxInf = float('inf')

        # Init dictionaries
        for path in self.pathCommodityDict:
            self.commodityInflow[path] = {edge: OrderedDict() for edge in self.network.edges}
            self.commodityOutflow[path] = {edge: OrderedDict() for edge in self.network.edges}
            self.hasBeenActive[path] = {edge: False for edge in self.network.edges}

        for e in self.network.edges:
            v, w = e
            tau = self.network[v][w]['transitTime']
            if e not in startingEdges:
                self.bIn[e][(minInf, t + tau)] = self.network[v][w]['inCapacity']  # Not full in this time frame
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

                    self.commodityOutflow[path][e][(minInf, t + tau)] = 0  # Here we can't know more than this
                elif self.edge_on_path(path, e):
                    # Edge is on path, but due to elif not the initial one
                    # Get the sum of transitTimes up to that edge
                    sTT, idx = 0, 0
                    while idx < len(path) - 1:
                        if path[idx] == v:
                            break
                        else:
                            sTT += self.network[path[idx]][path[idx + 1]]['transitTime']
                            idx += 1
                    self.commodityInflow[path][e][(minInf, t + sTT)] = 0
                    self.commodityOutflow[path][e][(minInf, t + sTT + tau)] = 0
                else:
                    # Edge never used by this commodity
                    self.commodityInflow[path][e][(minInf, maxInf)] = 0
                    self.commodityOutflow[path][e][(minInf, maxInf)] = 0

    def compute_alpha(self, theta, v, in_edges):
        """Computes alpha = min_{e in delta^-(v)} alpha_e such that commodity inflow constant and queue does not vanish"""
        alpha = float('inf')
        for e in in_edges:
            # Compute alpha to get extension size
            u, _ = e  # e = (u,v)
            tau = self.network[u][v]['transitTime']
            t = theta - tau

            # Find maximal alpha such that commodity flows stay constant
            for path in self.pathCommodityDict:
                if not self.edge_on_path(path, e):
                    continue
                partAlpha = 0
                lastRate = None
                found = False
                for interval, rate in self.commodityInflow[path][e].items():
                    t_l, t_u = interval
                    if Utilities.is_leq_tol(t_l, t) and Utilities.is_greater_tol(t_u, t):
                        lastRate = rate
                        partAlpha += t_u - t
                        found = True
                    elif found:
                        if lastRate == rate:
                            # We can extend further
                            partAlpha += t_u - t_l
                        else:
                            break

                alpha = min(alpha, partAlpha)

            # Find maximal alpha such that queues do not vanish
            qSize = self.queue_size(e, theta)
            if Utilities.is_greater_tol(qSize, 0.0):
                output = self.network[u][v]['outCapacity']  # Constant as long as there is queue
                input = self.inflow_rate(e, theta - tau)  # Note that this is constant in this alpha phase
                vanishRate = float(output - input)
                if Utilities.is_greater_tol(vanishRate, 0.0):
                    # The queue decreases.
                    # Otherwise the queue increases and hence we do not have to worry about it vanishing
                    alpha = min(alpha, qSize / vanishRate)

            phi = self.inverse_travel_time(e, theta)
            inflow_e = self.inflow_rate(e, phi)
            if Utilities.is_greater_tol(inflow_e, 0.0):
                for path in self.pathCommodityDict:
                    if self.edge_on_path(path, e):
                        inflow_ratio = float(self.inflow_rate(e, phi, commodityPath=path)) / inflow_e
                        outflow_com = inflow_ratio * self.network[u][v]['outCapacity']
                        if Utilities.is_greater_tol(outflow_com, 0.0):
                            timeToVanish = self.queue_size(e, theta, commodityPath=path) / outflow_com
                            if Utilities.is_greater_tol(timeToVanish, 0.0):
                                alpha = min(alpha, timeToVanish)
        return alpha

    def time_to_inflow_change(self, e, t, commodityPath=None):
        """Returns the maximum known alpha such that the inflow of e is constant on [t, t + alpha["""
        r = float('inf')
        paths = self.commodityInflow if not commodityPath else [commodityPath]
        for path in paths:
            if self.edge_on_path(path, e):
                for interval, flowRate in reversed(self.commodityInflow[path][e].items()):
                    t_l, t_u = interval
                    if Utilities.is_between_tol(t_l, t, t_u):
                        r = min(r, t_u - t)
                        break
        return r

    def cum_inflow(self, v, w, t, commodityPath=None):
        """
        :param v: tail of edge
        :param w: head of edge
        :param t: time
        :param commodityPath: If specified, compute for one commodity only
        :return: F_(v,w)^+(t)
        """
        if Utilities.is_leq_tol(t, 0):
            return 0
        e = (v, w)
        s = 0
        paths = self.commodityInflow if not commodityPath else [commodityPath]
        for path in paths:
            for interval, inflowVal in self.commodityInflow[path][e].items():
                t_l, t_u = interval
                if t_l == -float('inf') or t_u == float('inf'):
                    # In or outflow rate must be zero anyway
                    continue
                if Utilities.is_between_tol(t_l, t, t_u):
                    # This is the interval in which t lies
                    s += (t - t_l) * inflowVal
                elif t > t_u:
                    s += (t_u - t_l) * inflowVal
                elif t <= t_l:
                    break
        return s

    def cum_outflow(self, v, w, t, commodityPath=None):
        """
        :param v: tail of edge
        :param w: head of edge
        :param t: time
        :param commodityPath: If specified, compute for one commodity only
        :return: F_(v,w)^-(t)
        """
        if Utilities.is_leq_tol(t, 0):
            return 0
        e = (v, w)
        s = 0
        paths = self.commodityOutflow if not commodityPath else [commodityPath]
        for path in paths:
            for interval, outflowVal in self.commodityOutflow[path][e].items():
                t_l, t_u = interval
                if t_l == -float('inf') or t_u == float('inf'):
                    # In or outflow rate must be zero anyway
                    continue
                if Utilities.is_between_tol(t_l, t, t_u):
                    # This is the interval in which t lies
                    s += (t - t_l) * outflowVal
                elif t > t_u:
                    s += (t_u - t_l) * outflowVal
                elif t <= t_l:
                    break
        return s

    def queue_size(self, e, t, commodityPath=None):
        """Returns z_e(t). If commodyPath given, get just the values for that commodity"""
        v, w = e
        tau = self.network[v][w]['transitTime']
        qSize = float(self.cum_inflow(v, w, t - tau, commodityPath=commodityPath) - self.cum_outflow(v, w, t,
                                                                                                     commodityPath=commodityPath))
        assert (Utilities.is_geq_tol(qSize, 0.0))
        return qSize

    def inflow_rate(self, e, t, commodityPath=None):
        """Returns f^+_e(t)"""
        r = 0
        paths = self.commodityInflow if not commodityPath else [commodityPath]
        for path in paths:
            for interval, flowRate in reversed(self.commodityInflow[path][e].items()):
                t_l, t_u = interval
                if Utilities.is_between_tol(t_l, t, t_u):
                    # if t_l <= t < t_u:
                    r += flowRate
                    break
        return r

    def path_travel_time(self, path, t):
        """Returns the travel time along the entire path starting at time t"""
        if abs(t) == float('inf'):
            return t

        headArrivalTime = t
        for i in range(len(path) - 1):
            v, w = path[i], path[i + 1]
            e = (v, w)
            headArrivalTime = self.travel_time(e, headArrivalTime)
        return headArrivalTime - t

    def travel_time(self, e, t):
        """Returns T_e(t) = t + tau_e + q_e(t)"""
        if abs(t) == float('inf'):
            return t
        v, w = e
        tau = self.network[v][w]['transitTime']
        return t + tau + float(self.queue_size(e, t + tau)) / self.network[v][w]['outCapacity']

    def inverse_travel_time(self, e, theta):
        """Finds phi s.t. T_e(phi) = theta"""
        # Check whether we dont have a queue by chance
        v, w = e
        tau = self.network[v][w]['transitTime']
        if Utilities.is_eq_tol(self.travel_time(e, theta - tau), theta, tol=1e-7):
            return theta - tau
        # Find smallest time interval [a,b] s.t. in and outflow constant on [a,b] and theta \in [T_e(a), T_e(b)]
        # This allows us to linearly interpolate
        fc = lambda t: self.travel_time(e, t)
        a, b = -float('inf'), float('inf')
        T_a, T_b = fc(a), fc(b)
        for path in self.commodityInflow:
            if not self.edge_on_path(path, e):
                continue
            for idx, L in enumerate(
                    [reversed(self.commodityInflow[path][e].keys()), reversed(self.commodityOutflow[path][e].keys())]):
                for t_l, t_u in L:
                    if idx == 1:
                        # This is the outflow path, breakpoints hence differ
                        t_l, t_u = t_l - tau, t_u - tau
                    for bound in [t_l, t_u]:
                        if Utilities.is_between_tol(a, bound, b):
                            # We can make the interval smaller
                            z = fc(bound)
                            if Utilities.is_between_tol(T_a, theta, z):
                                b = bound
                                T_b = z
                            elif Utilities.is_between_tol(z, theta, T_b):
                                a = bound
                                T_a = z

        x = a + (b - a) * ((theta - T_a) / (T_b - T_a))
        assert (Utilities.is_eq_tol(theta, fc(x), tol=1e-8))
        return x

    def inverse_travel_time_old(self, e, theta):
        """Finds phi s.t. T_e(phi) = theta"""
        # Check whether we dont have a queue by chance
        v, w = e
        tau = self.network[v][w]['transitTime']
        if Utilities.is_eq_tol(self.travel_time(e, theta - tau), theta, tol=1e-6):
            return theta - tau

        for path in self.commodityInflow:
            if not self.edge_on_path(path, e):
                continue
            for interval, inflowVal in reversed(self.commodityInflow[path][e].items()):
                t_l, t_u = interval
                T_l, T_u = self.travel_time(e, t_l), self.travel_time(e, t_u)
                if T_l <= theta <= T_u:
                    if Utilities.is_eq_tol(T_l, theta, tol=1e-6):
                        phi = t_l
                    elif Utilities.is_eq_tol(T_u, theta, tol=1e-6):
                        phi = t_u
                    else:
                        # Must lie in this interval -> binary search
                        phi = Utilities.binary_search((t_l, t_u), lambda t: self.travel_time(e, t), theta)
                    return phi
                elif theta > T_u:
                    # Maybe check other path
                    break
                elif theta < T_l:
                    continue

    def generate_output(self, filePath, baseName):
        """Outputs the following:
        - Path travel times
        - Total cumulative inflow
        - Cumulative inflow per commodity
        """
        pathTTFile = os.path.join(filePath, baseName + "-" + "path_travel_times.txt")
        commodityCumFile = os.path.join(filePath, baseName + "-" + "cumulative_inflow_commodity_")
        totalCumFile = os.path.join(filePath, baseName + "-" + "total_cumulative_inflow.txt")

        # Cumulative inflow per commodity
        for idx, path in enumerate(self.pathCommodityDict):
            with open(commodityCumFile + str(idx) + ".txt", "w") as file:
                file.write("arc cumulative_inflow_function\n")
                for i in range(len(path) - 1):
                    v, w = path[i], path[i + 1]
                    e = (v, w)
                    s = str(v) + "," + str(w) + " "
                    breakPoints = [0]
                    for interval, rate in self.commodityInflow[path][e].items():
                        t_l, t_u = interval
                        if t_l >= 0:
                            breakPoints.append(t_l)
                        breakPoints.append(t_u)
                    breakPoints = sorted(list(set(breakPoints)))
                    tupleList = [(x, self.cum_inflow(v, w, x, commodityPath=path)) for x in breakPoints]
                    prettyList = Utilities.cleanup_output_list(tupleList)

                    s = s + ",".join([str(pair) for pair in prettyList]) + "\n"
                    file.write(s)

        # Total cumulative inflow
        with open(totalCumFile, "w") as file:
            file.write("arc total_cumulative_inflow_function\n")
            for e in self.network.edges:
                v, w = e
                s = str(v) + "," + str(w) + " "
                breakPoints = [0]
                for path in self.pathCommodityDict:
                    if not self.edge_on_path(path, e):
                        continue
                    for interval, rate in self.commodityInflow[path][e].items():
                        t_l, t_u = interval
                        if t_l >= 0:
                            breakPoints.append(t_l)
                        breakPoints.append(t_u)
                breakPoints = sorted(list(set(breakPoints)))
                tupleList = [(x, self.cum_inflow(v, w, x)) for x in breakPoints]
                prettyList = Utilities.cleanup_output_list(tupleList)

                s = s + ",".join([str(pair) for pair in prettyList]) + "\n"
                file.write(s)

        # Path travel times
        with open(pathTTFile, "w") as file:
            file.write("path path_travel_time\n")
            for path in self.pathCommodityDict:
                breakPoints = [0.0]
                for i in reversed(range(len(path) - 1)):
                    v, w = path[i], path[i + 1]
                    e = (v, w)
                    for xPath in self.pathCommodityDict:
                        if self.edge_on_path(xPath, e):
                            for interval in self.commodityInflow[xPath][e].keys():
                                t_l, t_u = interval
                                breakPoints = breakPoints + [t_l, t_u]
                            for interval in self.commodityOutflow[xPath][e].keys():
                                t_l, t_u = interval
                                breakPoints = breakPoints + [t_l, t_u]
                    breakPoints = sorted(list(set(breakPoints)))
                    breakPoints = [bp for bp in breakPoints if 0.0 <= bp < float('inf')]

                    breakPoints = breakPoints + [self.inverse_travel_time(e, bp) for bp in breakPoints]
                    breakPoints = sorted(list(set([bp for bp in breakPoints if 0.0 <= bp < float('inf')])))
                breakPoints = Utilities.get_unique_tol(L=breakPoints, tol=1e-5)

                s = ",".join([str(node) for node in path]) + " "
                tupleList = [(x, self.path_travel_time(path, x)) for x in breakPoints]
                prettyList = Utilities.cleanup_output_list(tupleList)
                prettyList = prettyList[:-1]
                prettyList.append((float('inf'), prettyList[-1][1]))
                if len(prettyList) > 2 and Utilities.is_eq_tol(prettyList[-1][1], prettyList[-3][1]):
                    prettyList.pop(-3)

                s = s + ",".join([str(pair) for pair in prettyList]) + "\n"
                file.write(s)

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
        self.init_values(initialTime, startingEdges)

        # Init priority heap (Note: isolated nodes not included)
        # Structure: Each element of the heap is a 3 tuple (theta, topologicalDistance, nodeName)
        # topologicalDistance (0, 1, 2, ...) is a tie breaker, lower values preferred
        # topologicalDistance: Topological ordering in reverse order (being far away top. is preferred)
        topologicalDistance = dict((node, 0) for node in self.network.nodes if node not in isoNodes)
        self.priority = [(initialTime, topologicalDistance[node], node) for node in self.network.nodes if
                         node not in isoNodes]
        heapq.heapify(self.priority)

        idx = 1
        startTime = timeit.default_timer()
        while self.priority:
            # Access first element of heap
            theta, topDist, v = heapq.heappop(self.priority)
            # STEP 1: Compute alpha extension size
            in_edges = list(self.network.in_edges(v))
            alpha = self.compute_alpha(theta, v, in_edges)

            # STEP 2: Extend outflow of incoming edges
            for e in in_edges:
                u, v = e
                tau = self.network[u][v]['transitTime']

                y = self.inflow_rate(e, theta - tau) if Utilities.is_eq_tol(self.queue_size(e, theta), 0.0) else float(
                    'inf')
                outflow_e = min(y, self.network[u][v]['outCapacity'])
                if Utilities.is_eq_tol(0, outflow_e, tol=1e-6):
                    for path in self.commodityOutflow:
                        upper = float('inf') if self.hasBeenActive[path][e] else theta + alpha
                        Utilities.dictInSort(self.commodityOutflow[path][e], (theta, upper, 0.0))
                else:
                    # We need to find phi s.t. T_e(phi) = theta
                    phi = self.inverse_travel_time(e, theta)
                    for path in self.commodityOutflow:
                        if not self.edge_on_path(path, e):
                            continue
                        inflow_e = self.inflow_rate(e, phi)
                        if Utilities.is_eq_tol(0, inflow_e, tol=1e-4):
                            upper = float('inf') if self.hasBeenActive[path][e] else theta + alpha
                            Utilities.dictInSort(self.commodityOutflow[path][e], (theta, upper, 0.0))
                        else:
                            inflow_ratio = float(self.inflow_rate(e, phi, commodityPath=path)) / inflow_e
                            outflow_ext = inflow_ratio * outflow_e
                            upper = float('inf') if self.hasBeenActive[path][e] and Utilities.is_eq_tol(outflow_ext,
                                                                                                        0.0) else theta + alpha
                            Utilities.dictInSort(self.commodityOutflow[path][e], (theta, upper, outflow_ext))
                            if Utilities.is_greater_tol(outflow_ext, 0.0): self.hasBeenActive[path][e] = True

                # STEP 3: Extend inflow rates along commodity paths
                for path in self.commodityOutflow:
                    edges_on_path = [(path[i], path[i + 1]) for i in range(len(path) - 1)]
                    try:
                        idx_on_path = edges_on_path.index(e)
                    except ValueError:
                        continue
                    if idx_on_path == len(edges_on_path) - 1:
                        # Last edge of commodity, nothing to do here
                        continue
                    else:
                        e_next = edges_on_path[idx_on_path + 1]
                        lastOutflowEntry = next(reversed(self.commodityOutflow[path][e].items()))
                        (n_l, n_u), r = lastOutflowEntry
                        Utilities.dictInSort(self.commodityInflow[path][e_next], (n_l, n_u, r))

                # STEP 4: Update bOut
                queue_in_future = self.queue_size(e, theta + alpha)
                y_ext = self.inflow_rate(e, theta - tau) if Utilities.is_eq_tol(queue_in_future, 0.0, 1e-6) else float(
                    'inf')
                bOut_ext = min(y_ext, self.network[u][v]['outCapacity'])
                Utilities.dictInSort(self.bOut[e], (theta, theta + alpha, bOut_ext))
            # STEP 5: Update priority queue
            theta_new = theta + alpha
            if theta_new < float('inf'):
                heapq.heappush(self.priority, (theta_new, topDist, v))

            if idx % 20 == 1 or len(self.priority) == 0:
                print("Iteration {0:d}: Node {1} | Theta {2:.2f} | Alpha {3:.2f} | {4:d} nodes in queue".format(idx, v,
                                                                                                                theta,
                                                                                                                alpha,
                                                                                                                len(
                                                                                                                    self.priority)))
            idx += 1
        endTime = timeit.default_timer()
        print("Done after {0:d} iterations. Elapsed time: {1:.2f} seconds.".format(idx - 1, endTime - startTime))
