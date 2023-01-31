from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from random import random


__all__ = ['quadmesh_herringbone_edges']


def quadmesh_edge_herringbone_neighbors(quadmesh, edge, orientation):
    """Find the other edges in the unit module of an initial edge for a herringbone pattern based on this quad mesh.
    The quad mesh must be double two-colorable, mainly with off-boundary singularities with valencies proportional to 4 (4, 8...).
    More on quad mesh coloring in Oval et al., 2021, Two-Colour Topology Finding of Quad-Mesh Patterns: https://www.sciencedirect.com/science/article/pii/S0010448521000415.

    Inputs
    ------
    quadmesh : Quad mesh
        A quad mesh with the valid topology.
    edge : tuple
        The identifier of the initial edge.
    orientation : bool
        The orientation of the herringbone pattern. 0 and 1 yield leftward and rightward around the initial edge compared to mesh normal, respectively.

    Returns
    -------
    edges : list
        List of edge identifiers that complete the unit module for a herringbone pattern. 

    """

    edge2type = {}
    for u, v in [edge, tuple(reversed(edge))]:
        if quadmesh.halfedge[u][v] is not None:

            if orientation == 0:
                t, u = quadmesh.halfedge_before(u, v)
                if quadmesh.halfedge[u][t] is not None:
                    edge2type[quadmesh.halfedge_after(u, t)] = 0
                v, w = quadmesh.halfedge_after(u, v)
                if quadmesh.halfedge[t][w] is not None:
                    edge2type[quadmesh.halfedge_after(t, w)] = 1

            elif orientation == 1:
                v, w = quadmesh.halfedge_after(u, v)
                if quadmesh.halfedge[w][v] is not None:
                    edge2type[quadmesh.halfedge_before(w, v)] = 1
                t, u = quadmesh.halfedge_before(u, v)
                if quadmesh.halfedge[t][w] is not None:
                    edge2type[quadmesh.halfedge_before(t, w)] = 0

    return edge2type


def quadmesh_herringbone_edges(quadmesh, edge0, orientation):
    """Find the edges in a quad mesh starting from an initial edge that form a herringbone pattern.
    The quad mesh must be double two-colorable, mainly with off-boundary singularities with valencies proportional to 4 (4, 8...).
    More on quad mesh coloring in Oval et al., 2021, Two-Colour Topology Finding of Quad-Mesh Patterns: https://www.sciencedirect.com/science/article/pii/S0010448521000415.

    Inputs
    ------
    quadmesh : Quad mesh
        A quad mesh with the valid topology.
    edge : tuple
        The identifier of the initial edge.
        orientation : bool
        The orientation of the herringbone pattern. 0 and 1 yield leftward and rightward around the initial edge compared to mesh normal, respectively.

    Returns
    -------
    edges : list
        List of edge identifiers that form a staggered pattern. 

    """

    edge2type = {(min(edge0), max(edge0)): orientation}
    to_visit, next_to_visit = [edge0], []

    while len(to_visit) != 0:
        for edge in to_visit:
            edge = (min(edge), max(edge))
            for edge, edge_type in quadmesh_edge_herringbone_neighbors(quadmesh, edge, edge2type[edge]).items():
                edge = (min(edge), max(edge))
                if edge not in edge2type:
                    edge2type[edge] = edge_type
                    next_to_visit.append(edge)
        to_visit, next_to_visit = next_to_visit, []

    return edge2type.keys()


def quadmesh_herringbone_unstructured_random_edges(quadmesh, nb_hb, candidate_edges=None):
    """Select random edges among a set of candiate edges for a random unstructured herringbone pattern.

    Inputs
    ------
    quadmesh : Quad mesh
        A quad mesh.
    nb_hb : int
        The number of edges to select.
    candidate_edges : list, optional
        A list of edges as tuples of vertex keys that forms the set of candidate edges to randomly select.

    Returns
    -------
    selected_edges : list
        List of edge identifiers that form a random unstructured herringbone pattern. 

    """

    if candidate_edges is None:
        candidate_edges = list(quadmesh.edges())

    selected_edges = []
    selected_faces = {fkey: False for fkey in quadmesh.faces()}
    count = quadmesh.number_of_edges()

    while nb_hb and count:

        count -= 1
        x = random()
        i = int(x * quadmesh.number_of_edges())
        edge = candidate_edges[i]
        fkeys = quadmesh.edge_faces(*edge)

        if any([selected_faces[fkey] for fkey in fkeys if fkey is not None]):
            continue

        selected_faces.update({fkey: True for fkey in fkeys})
        selected_edges.append(edge)
        nb_hb -= 1

    return selected_edges


def quadmesh_dome_herringbone_unstructured_equilength_edges(quadmesh, arc_length, length_offset=2):
    """Select edges in a dome to form an unstructured herringbone pattern with constant arc length of the herringbone arches.

    Inputs
    ------
    quadmesh : Quad mesh
        A quad mesh.
    arc_length : float
        The target length for the herringbone arches.
    length_offset : int, optional
        The offset between consecutive courses as the fraction of the arc length.

    Returns
    -------
    edges : list
        List of edge identifiers that form an unstructured herringbone pattern with quasi constant herringbone arch length. 

    """

    edges = []

    v0, u0 = mesh.edges_on_boundary()[0]
    v0, w0 = mesh.halfedge_after(u0, v0)
    p0 = mesh.collect_polyedge(v0, w0)

    for k, (u, v) in enumerate(pairwise(p0[:-1])):
        v, w = mesh.halfedge_after(u, v)
        p = mesh.collect_polyedge(v, w)
        m = len(p)
        length = sum([mesh.edge_length(a, b) for a, b in pairwise(p)])
        n = int(length / arch_length)

        for i in range(n):
            idx = int((i + (k % length_offset) / length_offset) / n * m)
            edges.append(p[idx: idx + 2])

    return edges
