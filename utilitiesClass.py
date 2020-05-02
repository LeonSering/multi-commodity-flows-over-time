# ===========================================================================
# Author:       Max ZIMMER
# Project:      multi-commodity-flows-over-time 2020
# File:         utilitiesClass.py
# Description:  Class containing utility functions
# ===========================================================================

import bisect
import os
import time

import matplotlib
import networkx as nx
import numpy as np
from collections import OrderedDict

# ======================================================================================================================


TOL = 1e-3  # Tolerance


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
    def is_between_tol(a, x, b,tol=TOL):
        """Returns true if a <= x <= b w.r.t tol"""
        return (Utilities.is_leq_tol(a, x) and Utilities.is_leq_tol(x, b))

    @staticmethod
    def get_insertion_point_left(L, el):
        """Get insertion point of element el in list L"""
        return bisect.bisect_left(L, el)

    @staticmethod
    def binary_search(interval, fc, y):
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

        mid = float(t_u+t_l)/2
        y_mid = fc(mid)
        if Utilities.is_eq_tol(y, y_mid, tol=1e-9):
            return mid  # Solution found
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
                t_l, t_u, _ = keyValList[-1]
                if t_u <= t_0:
                    keyValList.append((t_0, t_1, rate))
                elif t_0 <= t_u < t_1:
                    keyValList.append((t_u, t_1, rate))

        OD.clear()
        OD.update(OrderedDict([((triplet[0], triplet[1]), triplet[2]) for triplet in keyValList]))

    @staticmethod
    def get_edge_label_rotation(axes, src, dst, pos):
        """
        Modified from networkx plotting edge label function
        """
        x1, y1 = src
        x2, y2 = dst
        angle = np.arctan2(y2 - y1, x2 - x1) / (2.0 * np.pi) * 360  # degrees
        # make label orientation "right-side-up"
        if angle > 90:
            angle -= 180
        if angle < - 90:
            angle += 180
        # transform data coordinate angle to screen coordinate angle
        xy = np.array(pos)
        trans_angle = axes.transData.transform_angles(np.array((angle,)),
                                                      xy.reshape((1, 2)))[0]
        return trans_angle

    @staticmethod
    def get_shortest_path_network(network, labels=None):
        """
        Computes shortest path network
        :param network: the graph on which to compute the shortest path network
        :param labels: label functions. If None, then use edge transitTimes
        :return: networkx digraph: shortest path network of network given labels
        """

        if not labels:
            # Use transit-times as edge weight
            labels = nx.single_source_dijkstra_path_length(G=network, source='s',
                                                           weight='transitTime')  # Compute node distance from source

        # Create shortest path network containing _all_ shortest paths
        # shortestPathEdges = [(edge[0], edge[1]) for edge in network.edges() if labels[edge[0]] + network[edge[0]][edge[1]]['transitTime'] <= labels[edge[1]]]

        shortestPathEdges = [(edge[0], edge[1]) for edge in network.edges() if Utilities.is_geq_tol(labels[edge[1]],
                                                                                                    labels[edge[0]] +
                                                                                                    network[edge[0]][
                                                                                                        edge[1]][
                                                                                                        'transitTime'])]
        shortestPathNetwork = nx.DiGraph()
        shortestPathNetwork.add_nodes_from(network)
        shortestPathNetwork.add_edges_from(shortestPathEdges)

        for edge in shortestPathEdges:
            v, w = edge[0], edge[1]
            shortestPathNetwork[v][w]['outCapacity'] = network[v][w]['outCapacity']
            shortestPathNetwork[v][w]['transitTime'] = network[v][w]['transitTime']
            '''
            try:
                shortestPathNetwork[v][w]['inflowBound'] = network[v][w]['inflowBound'].items()[-1]
            except KeyError:
                continue
            '''

        for w in shortestPathNetwork:
            shortestPathNetwork.nodes[w]['dist'] = labels[w]
            shortestPathNetwork.nodes[w]['label'] = network.nodes[w]['label']
            shortestPathNetwork.nodes[w]['position'] = network.nodes[w]['position']

        return shortestPathNetwork

    @staticmethod
    def compute_min_attr_of_network(network, attr='outCapacity'):
        """
        Computes the minimal attr of all edges
        :param network: network with edge attribute attr
        :param attr: string-name of the attribute accessible through network[v][w][attr]
        :return: min(network[v][w][attr] | (v,w) in network.edges)
        """
        minAttr = float('inf')
        for edge in network.edges():
            v, w = edge[0], edge[1]
            minAttr = min([minAttr, network[v][w][attr]])
        return minAttr

    @staticmethod
    def join_intersect_dicts(dict1, dict2=None, *more_dicts):
        """Intersect two or more dict"""
        if not dict2:
            return dict1
        if len(more_dicts) == 0:
            return {key: (dict1[key], dict2[key]) for key in dict1 if key in dict2}
        else:
            d = {}
            for key in dict1:
                l = [dict1[key]]
                if key in dict2:
                    l.append(dict2[key])
                    for dAdditional in more_dicts:
                        if key not in dAdditional:
                            continue
                        l.append(dAdditional[key])
                    d[key] = tuple(l)
            return d

    @staticmethod
    def add_and_round_up(x, n):
        """Add n and round up to 10s"""
        x += n
        return x if x % 10 == 0 else x + 10 - x % 10

    @staticmethod
    def round_up(x):
        """Round up to 10s"""
        return x if x % 10 == 0 else x + 10 - x % 10

    #   All following functions are modifications of networkx functions to draw different edge styles
    '''Modified function networkx.draw_networkx_edges'''

    @staticmethod
    def draw_edges_with_boxes(G, pos,
                              edgelist=None,
                              width=1.0,
                              edge_color='k',
                              style='solid',
                              alpha=1.0,
                              edge_cmap=None,
                              edge_vmin=None,
                              edge_vmax=None,
                              ax=None,
                              arrows=True,
                              label=None,
                              **kwds):
        try:
            import matplotlib
            import matplotlib.pyplot as plt
            import matplotlib.cbook as cb
            from matplotlib.colors import colorConverter, Colormap
            from matplotlib.collections import LineCollection
            import numpy
        except ImportError:
            raise ImportError("Matplotlib required for draw()")
        except RuntimeError:
            print("Matplotlib unable to open display")
            raise

        if ax is None:
            ax = plt.gca()

        if edgelist is None:
            edgelist = list(G.edges())

        if not edgelist or len(edgelist) == 0:  # no edges!
            return None

        # set edge positions

        box_pos = numpy.asarray([(pos[e[0]], pos[e[1]]) for e in edgelist])
        p = 0.25
        edge_pos = []
        for edge in edgelist:
            src, dst = np.array(pos[edge[0]]), np.array(pos[edge[1]])
            s = dst - src
            # src = src + p * s  # Box at beginning
            # dst = src + (1-p) * s   # Box at the end
            dst = src  # No edge at all
            edge_pos.append((src, dst))
        edge_pos = numpy.asarray(edge_pos)

        if not cb.iterable(width):
            lw = (width,)
        else:
            lw = width

        if not cb.is_scalar_or_string(edge_color) \
                and cb.iterable(edge_color) \
                and len(edge_color) == len(edge_pos):
            if numpy.alltrue([cb.is_scalar_or_string(c)
                              for c in edge_color]):
                # (should check ALL elements)
                # list of color letters such as ['k','r','k',...]
                edge_colors = tuple([colorConverter.to_rgba(c, alpha)
                                     for c in edge_color])
            elif numpy.alltrue([not cb.is_scalar_or_string(c)
                                for c in edge_color]):
                # If color specs are given as (rgb) or (rgba) tuples, we're OK
                if numpy.alltrue([cb.iterable(c) and len(c) in (3, 4)
                                  for c in edge_color]):
                    edge_colors = tuple(edge_color)
                else:
                    # numbers (which are going to be mapped with a colormap)
                    edge_colors = None
            else:
                raise ValueError('edge_color must consist of either color names or numbers')
        else:
            if cb.is_scalar_or_string(edge_color) or len(edge_color) == 1:
                edge_colors = (colorConverter.to_rgba(edge_color, alpha),)
            else:
                raise ValueError(
                    'edge_color must be a single color or list of exactly m colors where m is the number or edges')
        edge_collection = LineCollection(edge_pos,
                                         colors=edge_colors,
                                         linewidths=width,
                                         antialiaseds=(1,),
                                         linestyle=style,
                                         transOffset=ax.transData,
                                         )

        edge_collection.set_zorder(1)  # edges go behind nodes
        edge_collection.set_label(label)
        ax.add_collection(edge_collection)

        # Note: there was a bug in mpl regarding the handling of alpha values for
        # each line in a LineCollection.  It was fixed in matplotlib in r7184 and
        # r7189 (June 6 2009).  We should then not set the alpha value globally,
        # since the user can instead provide per-edge alphas now.  Only set it
        # globally if provided as a scalar.
        if cb.is_numlike(alpha):
            edge_collection.set_alpha(alpha)

        if edge_colors is None:
            if edge_cmap is not None:
                assert (isinstance(edge_cmap, Colormap))
            edge_collection.set_array(numpy.asarray(edge_color))
            edge_collection.set_cmap(edge_cmap)
            if edge_vmin is not None or edge_vmax is not None:
                edge_collection.set_clim(edge_vmin, edge_vmax)
            else:
                edge_collection.autoscale()

        box_collection = Utilities.get_boxes(edge_colors=edge_colors, edge_pos=box_pos)
        box_collection.set_zorder(1)  # edges go behind nodes
        box_collection.set_label(label)
        ax.add_collection(box_collection)

        arrow_collection = Utilities.get_arrows_on_edges(edge_colors=edge_colors, edge_pos=box_pos)
        arrow_collection.set_zorder(2)
        if arrows:
            ax.add_collection(arrow_collection)

        return edge_collection, box_collection, arrow_collection

    '''Modified function networkx.draw_networkx_edges'''

    @staticmethod
    def draw_animation_edges(G, pos,
                             edgelist=None,
                             width=1.0,
                             edge_color='k',
                             style='solid',
                             alpha=1.0,
                             edge_cmap=None,
                             edge_vmin=None,
                             edge_vmax=None,
                             ax=None,
                             arrows=True,
                             label=None,
                             **kwds):
        try:
            import matplotlib
            import matplotlib.pyplot as plt
            import matplotlib.cbook as cb
            from matplotlib.colors import colorConverter, Colormap
            from matplotlib.collections import LineCollection
            import numpy
        except ImportError:
            raise ImportError("Matplotlib required for draw()")
        except RuntimeError:
            print("Matplotlib unable to open display")
            raise

        if ax is None:
            ax = plt.gca()

        if edgelist is None:
            edgelist = list(G.edges())

        if not edgelist or len(edgelist) == 0:  # no edges!
            return None

        # set edge positions

        box_pos = numpy.asarray([(pos[e[0]], pos[e[1]]) for e in edgelist])
        p = 0.25
        edge_pos = []
        for edge in edgelist:
            src, dst = np.array(pos[edge[0]]), np.array(pos[edge[1]])
            s = dst - src
            # src = src + p * s  # Box at beginning
            # dst = src + (1-p) * s   # Box at the end
            dst = src  # No edge at all
            edge_pos.append((src, dst))
        edge_pos = numpy.asarray(edge_pos)

        if not cb.iterable(width):
            lw = (width,)
        else:
            lw = width

        if not cb.is_scalar_or_string(edge_color) \
                and cb.iterable(edge_color) \
                and len(edge_color) == len(edge_pos):
            if numpy.alltrue([cb.is_scalar_or_string(c)
                              for c in edge_color]):
                # (should check ALL elements)
                # list of color letters such as ['k','r','k',...]
                edge_colors = tuple([colorConverter.to_rgba(c, alpha)
                                     for c in edge_color])
            elif numpy.alltrue([not cb.is_scalar_or_string(c)
                                for c in edge_color]):
                # If color specs are given as (rgb) or (rgba) tuples, we're OK
                if numpy.alltrue([cb.iterable(c) and len(c) in (3, 4)
                                  for c in edge_color]):
                    edge_colors = tuple(edge_color)
                else:
                    # numbers (which are going to be mapped with a colormap)
                    edge_colors = None
            else:
                raise ValueError('edge_color must consist of either color names or numbers')
        else:
            if cb.is_scalar_or_string(edge_color) or len(edge_color) == 1:
                edge_colors = (colorConverter.to_rgba(edge_color, alpha),)
            else:
                raise ValueError(
                    'edge_color must be a single color or list of exactly m colors where m is the number or edges')
        '''
        modEdgeColors = list(edge_colors)
        modEdgeColors = tuple(modEdgeColors + [colorConverter.to_rgba('w', alpha)
                                     for c in edge_color])
        #print(modEdgeColors)
        edge_collection = LineCollection(np.asarray(list(edge_pos)*2),
                                         colors=modEdgeColors,
                                         linewidths=[6]*len(list(edge_colors))+[4]*len(list(edge_colors)),
                                         antialiaseds=(1,),
                                         linestyle=style,
                                         transOffset=ax.transData,
                                         )
        '''
        edge_collection = LineCollection(edge_pos,
                                         colors=edge_colors,
                                         linewidths=6,
                                         antialiaseds=(1,),
                                         linestyle=style,
                                         transOffset=ax.transData,
                                         )

        edge_collection.set_zorder(1)  # edges go behind nodes
        edge_collection.set_label(label)
        ax.add_collection(edge_collection)

        tube_collection = LineCollection(edge_pos,
                                         colors=tuple([colorConverter.to_rgba('lightgrey', alpha)
                                                       for c in edge_color]),
                                         linewidths=4,
                                         antialiaseds=(1,),
                                         linestyle=style,
                                         transOffset=ax.transData,
                                         )

        tube_collection.set_zorder(1)  # edges go behind nodes
        tube_collection.set_label(label)
        ax.add_collection(tube_collection)
        # Note: there was a bug in mpl regarding the handling of alpha values for
        # each line in a LineCollection.  It was fixed in matplotlib in r7184 and
        # r7189 (June 6 2009).  We should then not set the alpha value globally,
        # since the user can instead provide per-edge alphas now.  Only set it
        # globally if provided as a scalar.
        if cb.is_numlike(alpha):
            edge_collection.set_alpha(alpha)

        if edge_colors is None:
            if edge_cmap is not None:
                assert (isinstance(edge_cmap, Colormap))
            edge_collection.set_array(numpy.asarray(edge_color))
            edge_collection.set_cmap(edge_cmap)
            if edge_vmin is not None or edge_vmax is not None:
                edge_collection.set_clim(edge_vmin, edge_vmax)
            else:
                edge_collection.autoscale()


        box_collection = Utilities.get_boxes(edge_colors=edge_colors, edge_pos=box_pos)
        box_collection.set_zorder(1)  # edges go behind nodes
        box_collection.set_label(label)
        ax.add_collection(box_collection)

        arrow_collection = Utilities.get_arrows_on_edges(edge_colors=edge_colors, edge_pos=box_pos)
        arrow_collection.set_zorder(0)

        if arrows:
            # Visualize them only if wanted
            ax.add_collection(arrow_collection)

        return edge_collection, box_collection, tube_collection, arrow_collection

    '''Modified function networkx.draw_networkx_edges'''

    @staticmethod
    def get_arrows_on_edges(edge_colors=None, edge_pos=None, width=1.0):
        import matplotlib.pyplot as plt
        import matplotlib.cbook as cb
        from matplotlib.colors import colorConverter
        ax = plt.gca()

        if not cb.iterable(width):
            lw = (width,)
        else:
            lw = width

        if edge_colors is None:
            edge_colors = tuple([colorConverter.to_rgba(c)
                                 for c in 'k'])

        arrows = []
        arrow_colors = edge_colors
        p = 0.3
        for src, dst in edge_pos:
            src = np.array(src)
            dst = np.array(dst)
            d = np.sqrt(np.sum(((dst - src) * p) ** 2))
            s = dst - src
            start = src

            if d == 0:  # source and target at same position
                continue
            #angle = np.rad2deg(np.arctan2(s[1], s[0]))
            #t = matplotlib.transforms.Affine2D().rotate_deg_around(start[0], start[1], angle)

            end = start + p*s
            arrow = FancyArrowPatch(start, end,
                                    arrowstyle='simple',
                                    shrinkA=2,
                                    shrinkB=2,
                                    mutation_scale=10,
                                    linewidth=10,
                                    connectionstyle=None,
                                    zorder=1)  # arrows go behind nodes
            arrows.append(arrow)

            # ax.add_patch(rec)
            # ax.plot([src[0], dst[0]], [src[1], dst[1]], lw=2, solid_capstyle="butt", zorder=0, color='r')

        arrow_collection = PatchCollection(arrows,
                                           linewidths=[ww for ww in lw],
                                           edgecolors=arrow_colors,
                                           facecolors='none',
                                           antialiaseds=(1,),
                                           transOffset=ax.transData, )

        return arrow_collection

    @staticmethod
    def get_boxes(edge_colors=None, edge_pos=None, width=1.0):
        import matplotlib.pyplot as plt
        import matplotlib.cbook as cb
        from matplotlib.colors import colorConverter
        ax = plt.gca()

        if not cb.iterable(width):
            lw = (width,)
        else:
            lw = width

        if edge_colors is None:
            edge_colors = tuple([colorConverter.to_rgba(c)
                                 for c in 'k'])

        rectangles = []
        box_colors = edge_colors
        # p = 0.25  # 1/4 of edge should be the box
        p = 1
        radius = 7
        for src, dst in edge_pos:
            src = np.array(src)
            dst = np.array(dst)
            d = np.sqrt(np.sum(((dst - src) * p) ** 2))
            s = dst - src
            # box_location = src  # Box at Beginning
            # box_location = src + (1-p)*s  # Box at End
            box_location = src  # Entire edge is box

            if d == 0:  # source and target at same position
                continue
            angle = np.rad2deg(np.arctan2(s[1], s[0]))
            delta = np.array([0, radius])
            t = matplotlib.transforms.Affine2D().rotate_deg_around(box_location[0], box_location[1], angle)
            rec = Rectangle(box_location - delta, width=d, height=radius * 2,
                            transform=t)
            rectangles.append(rec)

        box_collection = PatchCollection(rectangles,
                                           linewidths=[ww for ww in lw],
                                           edgecolors=box_colors,
                                           facecolors='none',
                                           antialiaseds=(1,),
                                           transOffset=ax.transData, )
        # alpha=0.5)
        return box_collection

    '''Modified function networkx.draw_networkx_nodes'''

    @staticmethod
    def draw_nodes(G,
                   pos,
                   nodelist=None,
                   node_size=300,
                   node_color='r',
                   node_shape='o',
                   alpha=1.0,
                   cmap=None,
                   vmin=None,
                   vmax=None,
                   ax=None,
                   linewidths=None,
                   label=None,
                   **kwds):

        try:
            import matplotlib.pyplot as plt
            import numpy
        except ImportError:
            raise ImportError("Matplotlib required for draw()")
        except RuntimeError:
            print("Matplotlib unable to open display")
            raise

        if ax is None:
            ax = plt.gca()

        if nodelist is None:
            nodelist = list(G.nodes())

        if not nodelist or len(nodelist) == 0:  # empty nodelist, no drawing
            return None

        try:
            xy = numpy.asarray([pos[v] for v in nodelist])
        except KeyError as e:
            raise nx.NetworkXError('Node %s has no position.' % e)
        except ValueError:
            raise nx.NetworkXError('Bad value in node positions.')

        radius = 7

        circles = []
        for v in nodelist:
            circ = Circle(pos[v], radius=radius, fill=True)
            circles.append(circ)

        node_collection = PatchCollection(circles, facecolors='r', edgecolors='k', linewidths=2, alpha=0.8)
        node_collection.set_zorder(3)
        ax.add_collection(node_collection)
        return node_collection
