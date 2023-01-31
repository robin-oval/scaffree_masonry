from __future__ import absolute_import
from __future__ import division
from __future__ import print_function


__all__ = [
    'mesh_merge_faces_by_edge',
    'mesh_merge_faces_by_edges'
]


def mesh_merge_faces_by_edge(mesh, edge):
    """Merge faces that are adjacent to each other in a mesh.

    Parameters
    ----------
    mesh : Mesh
            The mesh data structure.
    edge : tuple
            The identifier of the edge between the faces to merge.

    Returns
    -------
    int
            The identifier of the newly merge face.

    """

    u, v = edge
    f0, f1 = mesh.halfedge[u][v], mesh.halfedge[v][u]
    fv0, fv1 = mesh.face_vertices(f0), mesh.face_vertices(f1)
    i0 = fv0.index(v)
    fv0 = fv0[i0:] + fv0[:i0]
    i1 = fv1.index(u)
    fv1 = fv1[i1:] + fv1[:i1]
    fv = fv0 + fv1[1:-1]
    mesh.delete_face(f0)
    mesh.delete_face(f1)
    return mesh.add_face(fv)


def mesh_merge_faces_by_edges(mesh, edges):
    """Merge multiple faces that are adjacent to each other in a mesh.

    Parameters
    ----------
    mesh : Mesh
            The mesh data structure.
    edges : list
            The list of edge identifiers between the faces to merge.
            The edges can not share faces.

    Returns
    -------
    list
            The list of identifiers of the newly merge faces.

    """

    return [mesh_merge_faces_by_edge(mesh, edge) for edge in edges]

def parabola_arc_length(h, b):
    """Arc length of parabola f(x) = 4 * h / b ** 2 * x ** 2 between - b /2 and + b / 2."""

    from math import asinh

    def func(x):
        return 1 / 2 * x * (4 * x ** 2 + 1) ** 0.5 + 1 / 4 * asinh(2 * x)

    return b ** 2 / 2 / h * (func(2 * h / b) - func(0))
