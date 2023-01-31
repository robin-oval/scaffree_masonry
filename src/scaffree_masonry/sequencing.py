from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from compas.datastructures import Mesh, Network

from compas_quad.datastructures import QuadMesh

from compas.topology import breadth_first_ordering

from compas.utilities import pairwise


__all__ = ['Sequencer']


class Sequencer(QuadMesh):

    def __init__(self):
        super(Sequencer, self).__init__()
        self.default_face_attributes.update({'step': None, 'substep': None})

    @classmethod
    def from_faces(cls, mesh, initial_faces, min_nbrs=1, face2group={}):
        """Get sequence steps form initial faces. Store data in 'step' face attribute.

        Parameters
        ----------
        mesh : Mesh
                The mesh data structure.
        initial_faces : list
                The list of mesh face keys where the sequence starts.
        min_nbrs : int
                The minimum number of neighbors upstream in the sequence needed to add the face to the sequence.
        face2group : dict
                Dictionary of face keys pointing to their group of larger faces.

        Returns
        -------
        seq : Sequencer
                A Sequencer object.

        """

        def faces_neighbors(mesh, fkeys):
            fkeys = set(fkeys)
            nbrs = []
            for fkey in fkeys:
                for nbr in mesh.face_neighbors(fkey):
                    if nbr not in fkeys:
                        nbrs.append(nbr)
            return set(nbrs)

        seq = cls.from_vertices_and_faces(*mesh.to_vertices_and_faces())
        step = -1
        to_visit = initial_faces
        next_to_visit = []

        for fkey0 in to_visit:
            for fkey in face2group[fkey0]:
                seq.face_attribute(fkey, 'step', -1)

        while len(to_visit) != 0:
            step += 1

            for fkey0 in to_visit:
                seq.faces_attribute('step', step, face2group[fkey0])

            for fkey0 in to_visit:
                for nbr0 in faces_neighbors(seq, face2group[fkey0]):
                    if seq.face_attribute(nbr0, 'step') is None:
                        supp_nbr2 = []
                        for nbr2 in faces_neighbors(seq, face2group[nbr0]):
                            if seq.face_attribute(nbr2, 'step') is not None and seq.face_attribute(nbr2, 'step') != -1:
                                supp_nbr2.append(nbr2)
                        if len(set(supp_nbr2)) >= min_nbrs:
                            next_to_visit += face2group[nbr0]
                            seq.faces_attribute('step', -1, face2group[nbr0])

            to_visit = set(next_to_visit)
            next_to_visit = []

        return seq

    def number_of_steps(self):
        return len(set([self.face_attribute(fkey, 'step') for fkey in self.faces()]))

    def subsequencing_topo(self, face2group={}):
        """Get sequence of substeps per step. Store data in 'substep' face attribute.

        """

        # create graph for adjacency between faces within same step
        graph = Network()
        for fkey in self.faces():
            graph.add_node(key=fkey, attr_dict={key:value for key, value in zip('xyz', self.face_centroid(fkey))})
        for u, v in self.edges():
            if not self.is_edge_on_boundary(u, v):
                f0, f1 = self.halfedge[u][v], self.halfedge[v][u]
                if self.face_attribute(f0, 'step') == self.face_attribute(f1, 'step'):
                    graph.add_edge(f0, f1)
        adjacency = graph.adjacency

        # subsequence starting from leaves of graph
        leaves = set(graph.leaves())
        while len(leaves) != 0:
            leaf = leaves.pop()
            keys = breadth_first_ordering(adjacency, leaf)
            substep = -1
            for key in keys:
                if self.face_attribute(key, 'substep') is None:
                    substep += 1
                    for key2 in face2group[key]:                
                        self.face_attribute(key2, 'substep', substep)
                        if key in leaves:
                            leaves.remove(key)
    
    def subsequencing_geom(self, face2group={}):
        """Get sequence of substeps per step. Store data in 'substep' face attribute.

        """

        # create graph for adjacency between faces within same step
        step2faces = {step: [] for step in range(self.number_of_steps())}
        for fkey in self.faces():
            step2faces[self.face_attribute(fkey, 'step')].append(fkey)

        for step, faces in step2faces.items():
            substep = -1
            face2height = {fkey: self.face_centroid(fkey)[2] for fkey in faces}
            for fkey in sorted(face2height.keys(), key=face2height.get):
                if self.face_attribute(fkey, 'substep') is None:
                    substep += 1
                    for fkey2 in face2group[fkey]:                
                        self.face_attribute(fkey2, 'substep', substep)

    def assembly_at_step_k(self, k):
        return [fkey for fkey in self.faces() if self.face_attribute(fkey, 'step') == k]

    def assembly_until_step_k(self, k):
        return [fkey for fkey in self.faces() if self.face_attribute(fkey, 'step') <= k]

    def vertices_and_edges_at_step_k(self, k):
        vertices = {vkey: self.vertex_coordinates(vkey) for vkey in self.vertices() if self.vertex_attribute(vkey, 'step') <= k}
        edges = [(u, v) for u, v in self.edges() if self.vertex_attribute(u, 'step') <= k and self.vertex_attribute(v, 'step') <= k]
        return vertices, edges