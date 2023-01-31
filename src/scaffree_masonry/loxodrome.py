from __future__ import absolute_import
from __future__ import division
from __future__ import print_function


__all__ = ['quadmesh_edge_loxodromes']


def quadmesh_edge_front_right_edge(quadmesh, edge):
    """Get the edge at the front right of an input edge in a quad mesh.
    Front means in the direction of the face of the halfedge.
    Right means in the direction of the halfedge.

    Inputs
    ------
    quadmesh : Quad mesh
        A quad mesh.
    edge : tuple
        The identifier of the initial edge.

    Returns
    -------
    tuple
        The identifier of the front right edge. 

    """

    if quadmesh.halfedge_face(*edge) is not None:
        v, w = quadmesh.halfedge_after(*edge)
        if quadmesh.halfedge_face(w, v) is not None:
            return tuple(reversed(quadmesh.halfedge_before(w, v)))


def quadmesh_edge_front_left_edge(quadmesh, edge):
    """Get the edge at the front left of an input edge in a quad mesh.
    Front means in the direction of the face of the halfedge.
    Left means in the opposite direction of the halfedge.

    Inputs
    ------
    quadmesh : Quad mesh
        A quad mesh.
    edge : tuple
        The identifier of the initial edge.

    Returns
    -------
    tuple
        The identifier of the front left edge. 

    """

    if quadmesh.halfedge_face(*edge) is not None:
        t, u = quadmesh.halfedge_before(*edge)
        if quadmesh.halfedge_face(u, t) is not None:
            return tuple(reversed(quadmesh.halfedge_after(u, t)))


def quadmesh_edge_back_right_edge(quadmesh, edge):
    """Get the edge at the back right of an input edge in a quad mesh.
    Back means in the opposite direction of the face of the halfedge.
    Right means in the direction of the halfedge.

    Inputs
    ------
    quadmesh : Quad mesh
        A quad mesh.
    edge : tuple
        The identifier of the initial edge.

    Returns
    -------
    tuple
        The identifier of the back right edge. 

    """

    edge = tuple(reversed(edge))
    edge = quadmesh_edge_front_left_edge(quadmesh, edge)
    if edge is not None:
        edge = tuple(reversed(edge))
    return edge


def quadmesh_edge_back_left_edge(quadmesh, edge):
    """Get the edge at the back left of an input edge in a quad mesh.
    Back means in the opposite direction of the face of the halfedge.
    Left means in the opposite direction of the halfedge.

    Inputs
    ------
    quadmesh : Quad mesh
        A quad mesh.
    edge : tuple
        The identifier of the initial edge.

    Returns
    -------
    tuple
        The identifier of the back left edge. 

    """

    edge = tuple(reversed(edge))
    edge =  quadmesh_edge_front_right_edge(quadmesh, edge)
    if edge is not None:
        edge = tuple(reversed(edge))
    return edge


def quadmesh_edge_loxodrome_direction(quadmesh, edge0, func, max_edges=None):
    """Get the edges in the diagonal/loxodrome direction of an input edge in a quad mesh until the boundary or a maximum number.
    Use quadmesh_edge_front/back_right/left_edge functions.
    Front means in the direction of the face of the halfedge.
    Back means in the opposite direction of the face of the halfedge.
    Right means in the direction of the halfedge.
    Left means in the opposite direction of the halfedge.

    Inputs
    ------
    quadmesh : Quad mesh
        A quad mesh.
    edge0 : tuple
        The identifier of the initial edge.
    func : function
        One of the quadmesh_edge_front/back_right/left_edge functions.
    max_edges : int, opt
        The maximum number of edges to return. Infinite if None.

    Returns
    -------
    lox : list
        The list of edges in the specific diagonal/loxodrome direction. 

    """

    if max_edges is None:
        max_edges = float('inf')

    lox = [edge0]
    while max_edges:
        max_edges -= 1
        edge = func(quadmesh, lox[-1])
        if edge is None:
            break
        lox.append(edge)
        if quadmesh.is_edge_on_boundary(*edge):
            break
    return lox


def quadmesh_edge_loxodromes(quadmesh, edge0, front_right=True, front_left=True, back_right=True, back_left=True):
    """Get the edges in the diagonal/loxodrome directions of an input edge in a quad mesh until the boundary.
    Front means in the direction of the face of the halfedge.
    Back means in the opposite direction of the face of the halfedge.
    Right means in the direction of the halfedge.
    Left means in the opposite direction of the halfedge.

    Inputs
    ------
    quadmesh : Quad mesh
        A quad mesh.
    edge0 : tuple
        The identifier of the initial edge.
    front_right : bool
        Whether to collect the loxdrome in the front right direction. 
    front_left : bool
        Whether to collect the loxdrome in the front left direction. 
    back_right : bool
        Whether to collect the loxdrome in the back right direction. 
    back_left : bool
        Whether to collect the loxdrome in the back left direction. 

    Returns
    -------
    lox : list
        The list of edges in the specific diagonal/loxodrome direction. 

    """

    loxs = [edge0]
    
    if front_right:
        loxs += quadmesh_edge_loxodrome_direction(
            quadmesh, edge0, quadmesh_edge_front_right_edge)[1:]
    if front_left:
        loxs += quadmesh_edge_loxodrome_direction(
            quadmesh, edge0, quadmesh_edge_front_left_edge)[1:]
    if back_right:
        loxs += quadmesh_edge_loxodrome_direction(
            quadmesh, edge0, quadmesh_edge_back_right_edge)[1:]
    if back_left:
        loxs += quadmesh_edge_loxodrome_direction(
            quadmesh, edge0, quadmesh_edge_back_left_edge)[1:]
    return loxs
