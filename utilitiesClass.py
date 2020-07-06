# ===========================================================================
# Author:       Max ZIMMER
# Project:      multi-commodity-flows-over-time 2020
# File:         utilitiesClass.py
# Description:  Class containing utility functions
# ===========================================================================

import os
import time
from collections import OrderedDict

# ======================================================================================================================


TOL = 1e-6  # Tolerance


class Utilities:
    def __init__(self):
        pass

    @staticmethod
    def create_dir(path):
        """
        Create directory if not existing
        :param path: path to the directory
        """
        if not os.path.isdir(path):
            os.mkdir(path)

    @staticmethod
    def get_time():
        """Get formatted time"""
        return time.strftime("%d_%m_%Y-%H_%M_%S", time.localtime())

    @staticmethod
    def get_time_for_log():
        """Get formatted time"""
        return time.strftime("%H:%M:%S", time.localtime())

    @staticmethod
    def is_eq_tol(a, b, tol=TOL):
        """Is equal to with tolerance"""
        return abs(a - b) <= tol

    @staticmethod
    def is_not_eq_tol(a, b, tol=TOL):
        """Is not equal to with tolerance"""
        return abs(a - b) > tol

    @staticmethod
    def is_geq_tol(a, b, tol=TOL):
        """Is greater-or-equal than with tolerance"""
        return a - b + tol >= 0

    @staticmethod
    def is_leq_tol(a, b, tol=TOL):
        """Is less-or-equal than with tolerance"""
        return Utilities.is_geq_tol(b, a, tol)

    @staticmethod
    def is_greater_tol(a, b, tol=TOL):
        """Is 'a' greater than 'b' with tolerance"""
        return a - tol - b > 0

    @staticmethod
    def is_between_tol(a, x, b):
        """Returns true if a <= x <= b w.r.t tol"""
        return Utilities.is_leq_tol(a, x) and Utilities.is_leq_tol(x, b)

    @staticmethod
    def binary_search(interval, fc, y, tol=1e-6):
        # Try to find x \in interval s.t. fc(x) = y
        # Called recursively
        t_l, t_u = interval
        if t_l == -float('inf'):
            # We need to change the interval in a suitable manner
            t_l = 0
        if t_u == float('inf'):
            M = 10.0
            while fc(M) <= y:
                M *= 2
            t_u = M

        mid = float(t_u + t_l) / 2
        y_mid = fc(mid)
        if Utilities.is_eq_tol(y, y_mid, tol=tol):
            return mid  # Solution found
        elif t_u - t_l <= 1e-7:
            # The interval is sufficiently small to interpolate linearly
            T_l, T_u = fc(t_l), fc(t_u)
            x = t_l + (t_u - t_l) * ((y - T_l) / (T_u - T_l))
            return x
        elif y_mid < y:
            return Utilities.binary_search((mid, t_u), fc, y)
        elif y_mid > y:
            return Utilities.binary_search((t_l, mid), fc, y)

    @staticmethod
    def dictInSortAdd(OD, newValues):
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
                        idxShift += len(newL) - 1
                        t_0 = t_1
                        break  # Nothing else to do here
                    elif t_l <= t_0 < t_u < t_1:
                        lowSplit = (t_l < t_0)
                        newL = []
                        if lowSplit:
                            newL.append((t_l, t_0, r))
                        newL.append((t_0, t_u, r + rate))
                        keyValList[idx + idxShift:idx + idxShift + 1] = newL
                        idxShift += len(newL) - 1
                        t_0 = t_u  # Adjust interval for next iterations
                    else:
                        idx += 1
                if t_0 < t_1:
                    # Add to last case
                    keyValList.append((t_0, t_1, rate))

        OD.clear()
        OD.update(OrderedDict([((triplet[0], triplet[1]), triplet[2]) for triplet in keyValList]))

    @staticmethod
    def dictInSort(OD, newValues):
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
                t_l, t_u, r = keyValList[-1]
                if t_u <= t_0:
                    if t_u == t_0 and r == rate:
                        keyValList[-1] = (t_l, t_1, rate)
                    else:
                        keyValList.append((t_0, t_1, rate))
                elif t_0 <= t_u < t_1:
                    if r == rate:
                        keyValList[-1] = (t_l, t_1, rate)
                    else:
                        keyValList.append((t_u, t_1, rate))

        OD.clear()
        OD.update(OrderedDict([((triplet[0], triplet[1]), triplet[2]) for triplet in keyValList]))

    @staticmethod
    def get_unique_tol(L, tol=TOL):
        """Cleanup sorted list L from duplicates (wrt tolerance tol), keeping the order"""
        final = []
        while L:
            first = L.pop(0)
            final.append(first)
            for item in list(L):
                if Utilities.is_eq_tol(first, item, tol=tol):
                    L.remove(item)
                else:
                    # This is strictly bigger and hence all elements coming after that
                    break
        return final

    @staticmethod
    def cleanup_output_list(tupleList):
        finalList = []
        first_idx = 0
        while first_idx <= len(tupleList) - 1:
            x_0, y_0 = tupleList[first_idx]
            last_idx = first_idx
            while last_idx <= len(tupleList) - 1 and Utilities.is_eq_tol(y_0, tupleList[last_idx][1], tol=1e-4):
                last_idx += 1
            finalList.append((x_0, y_0))
            if last_idx - 1 > first_idx:
                x, y = tupleList[last_idx - 1]
                finalList.append((x, y))
            first_idx = last_idx

        # Find linear segments
        tupleList = finalList
        finalList = []
        slope = lambda a, b: float(b[1] - a[1]) / (b[0] - a[0])

        if len(tupleList) <= 2:
            first_idx = len(tupleList)
            finalList = tupleList
        else:
            first_idx = 0
        while first_idx <= len(tupleList) - 1:
            if first_idx == len(tupleList) - 1:
                finalList.append(tupleList[first_idx])
                break
            elif first_idx == len(tupleList) - 2:
                finalList.append(tupleList[first_idx + 1])
                break
            # first_idx <= len(tupleList) - 3
            x_0, y_0 = tupleList[first_idx]
            if first_idx == 0:
                finalList.append((x_0, y_0))
            last_idx = first_idx + 1
            x_1, y_1 = tupleList[last_idx]
            m_l = slope((x_0, y_0), (x_1, y_1))

            while last_idx <= len(tupleList) - 1 \
                    and Utilities.is_eq_tol(m_l, slope((x_0, y_0), tupleList[last_idx]), tol=1e-5) \
                    and tupleList[last_idx][0] < float('inf'):
                last_idx += 1

            x, y = tupleList[last_idx - 1]
            finalList.append((x, y))
            first_idx = last_idx - 1

        prettyList = []
        for x, y in finalList:
            if float('inf') > x == int(x):
                x = int(x)
            # elif x < float('inf'):
            #    x = float(str("%.4f" % x))
            if int(y) == y:
                y = int(y)
            # else:
            #    y = float(str("%.4f" % y))
            prettyList.append((x, y))
        return prettyList
