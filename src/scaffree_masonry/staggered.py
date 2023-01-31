from __future__ import absolute_import
from __future__ import division
from __future__ import print_function


__all__ = ['quadmesh_staggered_edges']


def quadmesh_edge_staggered_neighbors(quadmesh, edge, both_sides=True):
    """Find the other edges in the unit module of an initial edge for a staggered pattern based on this quad mesh.
    The quad mesh must be two-colorable, mainly with off-boundary singularities with valencies proportional to 2 (2, 4, 6, 8...).
    More on quad mesh coloring in Oval et al., 2021, Two-Colour Topology Finding of Quad-Mesh Patterns: https://www.sciencedirect.com/science/article/pii/S0010448521000415.

    Inputs
    ------
    quadmesh : Quad mesh
        A quad mesh with the valid topology.
    edge : tuple
        The identifier of the initial edge.
    both_sides : bool
        Whether the neighbors are collected in both sides of the halfedge or not.

    Returns
    -------
    edges : list
        List of edge identifiers that complete the unit module for a staggered pattern. 

    """

    edges = set()
    halfedges = [edge]
    if both_sides:
        halfedges.append(tuple(reversed(edge)))
    for u, v in halfedges:
        if quadmesh.halfedge[u][v] is not None:
            v, w = quadmesh.halfedge_after(u, v)
            if quadmesh.halfedge[w][v] is not None:
                edges.add(quadmesh.halfedge_before(w, v))
            t, u = quadmesh.halfedge_before(u, v)
            if quadmesh.halfedge[u][t] is not None:
                edges.add(quadmesh.halfedge_after(u, t))
    return edges


def quadmesh_staggered_edges(quadmesh, edge0, both_sides=True, flip=False, offset=False):
    """Find the edges in a quad mesh starting from an initial edge that form a staggered pattern.
    The quad mesh must be two-colorable, mainly with off-boundary singularities with valencies proportional to 2 (2, 4, 6, 8...).
    More on quad mesh coloring in Oval et al., 2021, Two-Colour Topology Finding of Quad-Mesh Patterns: https://www.sciencedirect.com/science/article/pii/S0010448521000415.

    Inputs
    ------
    quadmesh : Quad mesh
        A quad mesh with the valid topology.
    edge : tuple
        The identifier of the initial edge.
    both_sides : bool
        Whether the neighbors are collected in both sides of the halfedge or not.
    flip : bool, optional
        Whether to flip the staggered pattern.
    offset : bool, optional
        Whether to offset the staggered pattern.

    Returns
    -------
    edges : list
        List of edge identifiers that form a staggered pattern. 

    """

    if quadmesh.halfedge_face(*edge0) is None:
        edge0 = tuple(reversed(edge0))

    if flip:
        edge0 = quadmesh.halfedge_after(*edge0)

    if offset:
        edge0 = quadmesh.halfedge_after(*quadmesh.halfedge_after(*edge0))

    edges = set([(min(edge0), max(edge0))])
    to_visit, next_to_visit = [edge0], []

    while len(to_visit) != 0:
        for edge in to_visit:
            for edge in quadmesh_edge_staggered_neighbors(quadmesh, edge, both_sides):
                edge = (min(edge), max(edge))
                if edge not in edges:
                    edges.add(edge)
                    next_to_visit.append(edge)
        to_visit, next_to_visit = next_to_visit, []

    return edges
