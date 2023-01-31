from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from math import acos

from compas.geometry import distance_point_point, add_vectors, subtract_vectors, normalize_vector, scale_vector, dot_vectors, cross_vectors, centroid_points, angle_vectors_signed

from compas_quad.datastructures import Mesh, QuadMesh

from scaffree_masonry.sequencing import Sequencer

from scaffree_masonry.staggered import quadmesh_staggered_edges
from scaffree_masonry.herringbone import quadmesh_herringbone_edges
from scaffree_masonry.loxodrome import quadmesh_edge_loxodromes

from compas.topology import adjacency_from_edges, connected_components, dijkstra_distances

from compas_quad.coloring import quad_mesh_polyedge_2_coloring

from compas.utilities import pairwise


__all__ = ['MasonryPattern']


class MasonryPattern(QuadMesh):

    def __init__(self):
        super(MasonryPattern, self).__init__()
        self.default_face_attributes.update({'face2block': None})
        self.default_face_attributes.update({'is_supported': False})
        self.attributes['block2faces'] = {}
        self.attributes['block2type'] = {}
        self.attributes['block2volume'] = {}
        self.attributes['block2step'] = {}
        self.attributes['block2substep'] = {}

    def dummy_blocks(self, faces=None):
        if faces is None:
            faces = self.faces()
        bkey = self.max_block_key() + 1
        bkey0 = bkey
        for fkey in faces:
            if self.face_block(fkey) is None:
                self.face_block(fkey, bkey)
                self.block_faces(bkey, [fkey])
                self.block_type(bkey, 'dummy')
                bkey += 1
        print('{} new dummy blocks'.format(bkey - bkey0))

    def staggered_blocks(self, edge0, both_sides=True, flip=False, offset=False, faces_to_skip=None):
        if faces_to_skip is None:
            faces_to_skip = set()
        bkey = self.max_block_key() + 1
        bkey0 = bkey
        edges = quadmesh_staggered_edges(self, edge0, both_sides, flip, offset)
        for edge in edges:
            if any([fkey in faces_to_skip for fkey in self.edge_faces(*edge)]):
                continue
            if not self.is_edge_on_boundary(*edge):
                fkeys = self.edge_faces(*edge)
                self.store_block_data(bkey, fkeys, 'hb')
                bkey += 1
        print('{} new staggered blocks'.format(bkey - bkey0))
        return edges

    def herringbone_blocks(self, edge0, orientation=0):
        bkey = self.max_block_key() + 1
        bkey0 = bkey
        edges = quadmesh_herringbone_edges(self, edge0, orientation)
        for edge in edges:
            if not self.is_edge_on_boundary(*edge):
                fkeys = self.edge_faces(*edge)
                self.store_block_data(bkey, fkeys, 'hb')
                bkey += 1
        print('{} new staggered blocks'.format(bkey - bkey0))
        return edges

    def loxodrome_blocks(self, edges0, front_right=True, front_left=True, back_right=True, back_left=True):
        bkey = self.max_block_key() + 1
        bkey0 = bkey
        all_edges = set()
        for edge0 in edges0:
            edges = [(min(u, v), max(u, v)) for u, v in quadmesh_edge_loxodromes(self, edge0, front_right, front_left, back_right, back_left)]
            for edge in edges:
                if edge not in all_edges:
                    if not self.is_edge_on_boundary(*edge):
                        fkeys = self.edge_faces(*edge)
                        self.store_block_data(bkey, fkeys, 'lox')
                        bkey += 1
            all_edges.update(edges)
        print('{} new loxdrome blocks'.format(bkey - bkey0))
        return all_edges

    def create_blocks(self, edges, faces_to_skip=None, name=None):
        if faces_to_skip is None:
            faces_to_skip = set()
        bkey = self.max_block_key() + 1
        bkey0 = bkey
        for edge in edges:
            if any([fkey in faces_to_skip for fkey in self.edge_faces(*edge)]):
                continue
            if not self.is_edge_on_boundary(*edge):
                fkeys = self.edge_faces(*edge)
                self.store_block_data(bkey, fkeys, name)
                bkey += 1
        print('{} new blocks'.format(bkey - bkey0))
        return edges

    def loxodrome_patch(self, polyedge0, nb_lox, directions=0):
        """

        Inputs
        ------
        polyedge0 : list
            The polyedge as a list of successive vertices pointing towards the quad mesh patch.
        nb_lox : int
            The number of loxodromes starting from the initial polyedge.
        directions : int
            The direction of the loxdromes. 0: left and right; 1: right; 2: left.

        """

        back_right, back_left = False, False
        if directions == 0:
            front_right, front_left = True, True
        elif directions == 1:
            front_right, front_left = True, False
        elif directions == 2:
            front_right, front_left = False, True

        n = len(polyedge0)

        edges0 = []
        spacing = int((n - nb_lox + 1) / nb_lox)
        for k in range(nb_lox):
            i = int((spacing - 1) / 2 + k * (spacing + 1))
            edges0.append(tuple(polyedge0[i:i + 2]))
        
        edges1, edges2 = [], []
        if polyedge0[0] != polyedge0[-1]:
            
            strip_first = self.collect_strip(polyedge0[0], polyedge0[1], both_sides=False)
            m = len(strip_first)
            for k in range(int(m / n * nb_lox)):
                i = int((spacing - 1) / 2 + k * (spacing + 1))
                edges1.append(strip_first[i])

            strip_last = self.collect_strip(polyedge0[-2], polyedge0[-1], both_sides=False)
            m = len(strip_last)
            for k in range(int(m / n * nb_lox)):
                i = int((spacing - 1) / 2 + k * (spacing + 1))
                edges2.append(strip_last[i])

        lox_blocks0 = self.loxodrome_blocks(edges0, front_right, front_left, back_right, back_left)
        lox_blocks1 = self.loxodrome_blocks(edges1, front_right, front_left, back_right, back_left)
        lox_blocks2 = self.loxodrome_blocks(edges2, front_right, front_left, back_right, back_left)
        
        lox_blocks0.update(lox_blocks1)
        lox_blocks0.update(lox_blocks2)
        lox_blocks = set()
        for u, v in lox_blocks0:
                lox_blocks.add((min(u, v), max(u, v)))
        
        return lox_blocks

    ### MASONRY VOLUMES ###

    def block_volumes(self, block2thickness, triangulate=False, block2ecc=None):
        """Produce block volumes as closed meshes. Store data in blocks.

        Parameters
        ----------
        block2thickness : dict
            Dictionary with block keys pointing to block thickness values.
        triangulate : bool
            Whether to triangulate the block meshes.
        ecc : float
            Eccentricity distance of block compared to medial mesh. Positive values create an eccentricity in the direction of the mesh normal, and vice versa.

        """

        if block2ecc is None:
            block2ecc = {bkey: 0.0 for bkey in self.blocks()}

        block2volume = {}
        to_visit = set(list(self.faces()))
        vkey2normal = {vkey: self.vertex_normal(vkey) for vkey in self.vertices()}
        
        while to_visit:
        
            fkey = to_visit.pop()
            bkey = self.face_block(fkey)
            thickness = block2thickness[bkey]
            ecc = block2ecc[bkey]
            bfkeys = self.block_faces(bkey)
        
            # border
            if len(bfkeys) == 1:
                border = self.face_vertices(fkey)
            else:
                to_visit.difference_update(bfkeys)
                fkey1, fkey2 = bfkeys
                u, v = self.face_adjacency_halfedge(fkey1, fkey2)
                border = [v, self.face_vertex_descendant(fkey1, v), self.face_vertex_descendant(
                    fkey1, v, n=2), u, self.face_vertex_descendant(fkey2, u), self.face_vertex_descendant(fkey2, u, n=2)]
        
            # offset (thickness + eccentricity)

            # vertices
            top, bottom = [], []
            for vkey in border:
                normal = vkey2normal[vkey]
                xyz0 = self.vertex_coordinates(vkey)
                xyz = add_vectors(xyz0, scale_vector(normal, thickness / 2 + ecc))
                top.append(xyz)
                xyz = add_vectors(xyz0, scale_vector(normal, - thickness / 2 + ecc))
                bottom.append(xyz)
            block_vertices = bottom + top
            
            # faces
            block_faces = []
            n = int(len(block_vertices) / 2)
        
            # extrados and intrados faces
            bottom = list(range(0, n))
            top = list(range(n, 2 * n))
            if triangulate:
                block_faces += [[v, u, bottom[0]]
                                for u, v in pairwise(bottom[1:])]
                block_faces += [[top[0], u, v] for u, v in pairwise(top[1:])]
            else:
                block_faces += [list(reversed(bottom)), top]
        
            # lateral faces
            for k in range(n):
                face = [bottom[k], bottom[(k + 1) %
                                          n], top[(k + 1) % n], top[k]]
                if triangulate:
                    block_faces += [[face[0], u, v]
                                    for u, v in pairwise(face[1:])]
                else:
                    block_faces.append(face)
            
            bkey = self.face_block(fkey)
            bmesh = Mesh.from_vertices_and_faces(block_vertices, block_faces)
            self.block_volume(bkey, bmesh)
            block2volume[bkey] = bmesh
            
        return block2volume

    ### MASONRY ASSEMBLY ###

    def identify_supports(self, supported_vertices):
        for fkey in self.faces():
            supp_vkeys = [vkey for vkey in self.face_vertices(
                fkey) if vkey in supported_vertices]
            if len(supp_vkeys) >= 2:
                for bfkey in self.block_faces(self.face_block(fkey)):
                    self.face_attribute(bfkey, 'is_supported', True)

    def assembly(self, include_substeps=False):
        """Compute assembly sequence propagating from supports. Store data in 'block2step'."""

        supported_faces = [fkey for fkey in self.faces(
        ) if self.face_attribute(fkey, 'is_supported')]
        face2group = {fkey: self.face_block_faces(
            fkey) for fkey in self.faces()}

        seq = Sequencer.from_faces(
            self, supported_faces, min_nbrs=1, face2group=face2group)

        for bkey, fkeys in self.blocks(data=True):
            self.block_step(bkey, seq.face_attribute(fkeys[0], 'step'))

        if include_substeps:
            seq.subsequencing_geom(face2group=face2group)
            for bkey, fkeys in self.blocks(data=True):
                self.block_substep(
                    bkey, seq.face_attribute(fkeys[0], 'substep'))

    def assembly_polyedges_twocolor(self, supported_pkeys):
        
        edge2pkey = {(u, v): pkey for pkey, polyedge in self.polyedges(data=True) for u0, v0 in pairwise(polyedge) for u, v in ((u0, v0), (v0, u0))}

        pkey0 = supported_pkeys.pop()
        pkey2color = quad_mesh_polyedge_2_coloring(self)
        color0 = pkey2color[pkey0]
        pkeys0 = set([pkey for pkey in self.polyedges() if pkey2color[pkey] == color0])

        edges = set()
        for fkey in self.faces():
            if self.face_degree(fkey) == 4:
                u, v, w, x = list(self.face_vertices(fkey))
                for a, b, c, d in [(u, v, w, x), (v, w, x, u)]:
                    if edge2pkey[(a, b)] in pkeys0:
                        pab, pcd = edge2pkey[(a, b)], edge2pkey[(c, d)]
                        p0, p1 = sorted([pab, pcd])
                        if (p0, p1) not in edges:
                            edges.add((p0, p1))
                        break

        convert_pkeys = {pkey: pkey if pkey not in supported_pkeys else pkey0 for pkey in self.polyedges()}
        edges = set(tuple(sorted([convert_pkeys[p0], convert_pkeys[p1]])) for p0, p1 in edges)

        adjacency = adjacency_from_edges(edges)
        weight = {edge: 1.0 for u, v in edges for edge in [(u, v), (v, u)]}
        dist_pkey2pkey0 = dijkstra_distances(adjacency, weight, pkey0)
        for pkey in pkeys0:
            if pkey not in dist_pkey2pkey0:
                dist_pkey2pkey0[pkey] = 999
        
        pkey2step = {pkey: dist_pkey2pkey0[pkey] if pkey not in supported_pkeys else 0 for pkey in pkeys0}
        max_step = int(max(pkey2step.values()))

        fkey2pkey = {fkey: [edge2pkey[edge] for edge in self.face_halfedges(fkey)] for fkey in self.faces()}

        step2bkeys = {i: [] for i in range(max_step + 1)}
        for bkey, fkeys in self.blocks(data=True):
            pkeys = [pkey for fkey in fkeys for pkey in fkey2pkey[fkey]]
            step = int(min([pkey2step[pkey] for pkey in pkeys if pkey in pkey2step]))
            self.block_step(bkey, step)
            step2bkeys[step].append(bkey)

    def assembly_type_dome(self, include_substeps=False):

        edge2pkey = {(u, v): pkey for pkey, polyedge in self.polyedges(data=True) for u0, v0 in pairwise(polyedge) for u, v in ((u0, v0), (v0, u0))}
        
        cpkey2centroid = {pkey: centroid_points([self.vertex_coordinates(vkey) for vkey in polyedge]) for pkey, polyedge in self.polyedges(data=True) if polyedge[0] == polyedge[-1]}

        sorted_cpkey = sorted(cpkey2centroid.items(), key=lambda item: item[1][2])
        cpkey2step = {pkey: i for i, (pkey, centroid) in enumerate(sorted_cpkey)}

        fkey2pkey = {fkey: [edge2pkey[edge] for edge in self.face_halfedges(fkey)] for fkey in self.faces()}

        is_bkey_herringbone = {}
        step2bkeys = {}
        for bkey, fkeys in self.blocks(data=True):
            pkeys = [pkey for fkey in fkeys for pkey in fkey2pkey[fkey]]
            is_bkey_herringbone[bkey] = len(pkeys) > 2
            step = int(min([cpkey2step[pkey] for pkey in pkeys if pkey in cpkey2step]))
            self.block_step(bkey, step)
            if step in step2bkeys:
                step2bkeys[step].append(bkey)
            else:
                step2bkeys[step] = [bkey]

        if include_substeps:
            o = self.centroid()
            for step, bkeys in step2bkeys.items():
                fkeys = set([fkey for bkey in bkeys for fkey in self.block_faces(bkey)])
                
                edges = []
                for fkey in fkeys:
                    for fkey2 in self.face_neighbors(fkey):
                        if fkey2 in fkeys:
                            edges.append((fkey, fkey2))
                adjacency = adjacency_from_edges(edges)
                for fkey in fkeys:
                    if fkey not in adjacency:
                        adjacency[fkey] = []
                face_groups = connected_components(adjacency)
                
                gidx2angle = {}
                for i, fkeys in enumerate(face_groups):
                    g = centroid_points([self.face_centroid(fkey) for fkey in fkeys])
                    ox = [1.0, 0.0, 0.0]
                    og = subtract_vectors(g, o)
                    og[2] = 0.0
                    theta = angle_vectors_signed(ox, og, [0.0 ,0.0, 1.0])
                    gidx2angle[i] = theta
                
                fkey2substep = {}
                sorted_gidx = sorted(gidx2angle.items(), key=lambda item: item[1])
                for substep, (i, angle) in enumerate(sorted_gidx):
                    for fkey in face_groups[i]:
                        fkey2substep[fkey] = substep
                
                for fkey, substep in fkey2substep.items():
                    bkey = self.face_block(fkey)
                    self.block_substep(bkey, fkey2substep[fkey])
                
        for bkey in self.blocks():
            if bkey not in self.attributes['block2substep']:
                self.block_substep(bkey, -1)

    def assembly_type_barrel(self, supported_vertices, pkey0, include_substeps=False):

        edge2pkey = {(u, v): pkey for pkey, polyedge in self.polyedges(data=True) for u0, v0 in pairwise(polyedge) for u, v in ((u0, v0), (v0, u0))}

        pkey2color = quad_mesh_polyedge_2_coloring(self)
        color0 = pkey2color[pkey0]
        pkeys0 = set([pkey for pkey in self.polyedges() if pkey2color[pkey] == color0])

        edges = set()
        for fkey in self.faces():
            u, v, w, x = list(self.face_vertices(fkey))
            for a, b, c, d in [(u, v, w, x), (v, w, x, u)]:
                if edge2pkey[(a, b)] in pkeys0:
                    pab, pcd = edge2pkey[(a, b)], edge2pkey[(c, d)]
                    p0, p1 = sorted([pab, pcd])
                    if (p0, p1) not in edges:
                        edges.add((p0, p1))
                    break
        
        adjacency = adjacency_from_edges(edges)
        weight = {edge: 1.0 for u, v in edges for edge in [(u, v), (v, u)]}
        dist_pkey2pkey0 = dijkstra_distances(adjacency, weight, pkey0)
        
        pkey2step = {pkey: dist_pkey2pkey0[pkey] if pkey != pkey0 else 0 for pkey in pkeys0}
        max_step = int(max(pkey2step.values()))

        fkey2pkey = {fkey: [edge2pkey[edge] for edge in self.face_halfedges(fkey)] for fkey in self.faces()}

        step2bkeys = {i: [] for i in range(max_step + 1)}
        for bkey, fkeys in self.blocks(data=True):
            pkeys = [pkey for fkey in fkeys for pkey in fkey2pkey[fkey]]
            step = int(min([pkey2step[pkey] for pkey in pkeys if pkey in pkey2step]))
            self.block_step(bkey, step)
            step2bkeys[step].append(bkey)

        if include_substeps:
            for step, bkeys in step2bkeys.items():
                fkeys = set([fkey for bkey in bkeys for fkey in self.block_faces(bkey)])
                
                edges = []
                for fkey in fkeys:
                    for fkey2 in self.face_neighbors(fkey):
                        if fkey2 in fkeys:
                            edges.append((fkey, fkey2))
                adjacency = adjacency_from_edges(edges)
                for fkey in fkeys:
                    if fkey not in adjacency:
                        adjacency[fkey] = []
                face_groups = connected_components(adjacency)
                if step == 0:
                    print(face_groups, len(face_groups[0]))

                included_faces = set([face for faces in face_groups for face in faces])
                # print(len(fkeys), len(included_faces))
                for fkey in fkeys:
                    if fkey not in included_faces:
                        print('!')
                        face_groups.append(fkey)

                
                gidx2xz = {}
                for i, fkeys in enumerate(face_groups):
                    # x = min([self.vertex_coordinates(vkey)[0] for fkey in fkeys for vkey in self.face_vertices(fkey)])
                    # z = min([self.vertex_coordinates(vkey)[2] for fkey in fkeys for vkey in self.face_vertices(fkey)])
                    x, y, z = centroid_points([self.vertex_coordinates(vkey) for fkey in fkeys for vkey in self.face_vertices(fkey)])
                    gidx2xz[i] = (x, z)
                
                fkey2substep = {}
                sorted_gidx = sorted(gidx2xz.items(), key=lambda item: (item[1][1], item[1][0]))
                for substep, (i, xz) in enumerate(sorted_gidx):
                    for fkey in face_groups[i]:
                        fkey2substep[fkey] = substep
                if step == 0:
                    print(fkey2substep, len(fkey2substep))
                for fkey, substep in fkey2substep.items():
                    bkey = self.face_block(fkey)
                    self.block_substep(bkey, fkey2substep[fkey])
                
        for bkey in self.blocks():
            if bkey not in self.attributes['block2substep']:
                fkeys = self.block_faces(bkey)
                print('no substep for bkey with fkeys', fkeys, self.block_step(bkey))
                self.block_substep(bkey, -1)

    ### UTILITIES ###

    def max_block_key(self):
        list_blocks = list(self.blocks())
        if len(list_blocks) == 0:
            return -1
        else:
            return max(list_blocks)

    def number_of_blocks(self):
        return len(list(self.blocks()))
            
    def number_of_assembly_steps(self):
        return len(set([self.block_step(bkey) for bkey in self.blocks()]))

    ### iterators ###

    def blocks(self, data=False):
        for bkey in self.attributes['block2faces']:
            if data:
                yield bkey, self.block_faces(bkey)
            else:
                yield bkey

    def supported_faces(self):
        for fkey in self.faces():
            if self.face_attribute(fkey, 'is_supported'):
                yield fkey

    def supported_blocks(self):
        for bkey in self.blocks():
            if self.face_attribute(self.block_faces(bkey)[0], 'is_supported'):
                yield bkey

    ### getters/setters ###

    def face_block(self, fkey, value=None):
        if value is not None:
            self.face_attribute(fkey, 'face2block', value)
        return self.face_attribute(fkey, 'face2block')

    def block_faces(self, bkey, value=None):
        if value is not None:
            self.attributes['block2faces'][bkey] = value
        return self.attributes['block2faces'][bkey]

    def block_type(self, bkey, value=None):
        if value is not None:
            self.attributes['block2type'][bkey] = value
        return self.attributes['block2type'][bkey]

    def block_volume(self, bkey, value=None):
        if value is not None:
            self.attributes['block2volume'][bkey] = value
        return self.attributes['block2volume'][bkey]

    def block_step(self, bkey, value=None):
        if value is not None:
            self.attributes['block2step'][bkey] = value
        return self.attributes['block2step'][bkey]

    def block_substep(self, bkey, value=None):
        if value is not None:
            self.attributes['block2substep'][bkey] = value
        return self.attributes['block2substep'][bkey]

    def store_block_data(self, block_key, block_faces, block_type):
        self.block_faces(block_key, block_faces)
        self.block_type(block_key, block_type)
        for face in block_faces:
            self.face_block(face, block_key)

    def delete_block_data(self, bkey):
        for fkey in self.block_faces(bkey):
            self.face_attribute(fkey, 'face2block', None)
        del self.attributes['block2faces'][bkey]
        del self.attributes['block2type'][bkey]

    def number_of_herringbone_blocks(self):
        return len([bkey for bkey in self.blocks() if self.block_type(bkey) == 'hb'])

    def number_of_platebande_blocks(self):
        return len([bkey for bkey in self.blocks() if self.block_type(bkey) == 'pb'])

    def face_block_faces(self, fkey):
        return self.block_faces(self.face_block(fkey))

    ### analysis ###

    def block_neighbors(self, bkey):

        blocks = []
        for fkey in self.block_faces(bkey):
            for fkey_nbr in self.face_neighbors(fkey):
                bkey_nbr = self.face_block(fkey_nbr)
                if bkey_nbr != bkey:
                    blocks.append(bkey_nbr)
        return set(blocks)

    def block_adjacency_edge(self, bkey0, bkey1):

        for fkey0 in self.blocks_faces(bkey0):
            for fkey1 in self.block_faces(bkey1):
                if self.face_adjacency_halfedge(fkey0, fkey1) is not None:
                    return self.face_adjacency_halfedge(fkey0)

    def edge_normal(self, edge, area_weight=False):

        u, v = edge
        nu, nv = self.vertex_normal(u), self.vertex_normal(v)
        if not area_weight:
            return scale_vector(add_vectors(nu, nv), 0.5)
        else:
            au, av = self.vertex_area(u), self.vertex_area(v)
            a = au + av
            return add_vectors(scale_vector(nu, au / a), scale_vector(nv, av / a))

    def interface_angles(self):

        z0 = [0.0, 0.0, 1.0]

        interface_angles = {}

        for u, v in self.edges():
            if not self.is_edge_on_boundary(u, v):
                f0, f1 = self.halfedge[u][v], self.halfedge[v][u]
                b0 = self.face_block(f0)
                b1 = self.face_block(f1)
                if b0 != b1:

                    step0, step1 = self.block_step(b0), self.block_step(b1)
                    if step0 == step1:
                        continue
                    if step0 > step1:
                        u, v = v, u
                        f0, f1 = f1, f0
                        b0, b1 = b1, b0
                        step0, step1 = step1, step0

                    fc0, fc1 = self.face_centroid(f0), self.face_centroid(f1)
                    n01 = normalize_vector(subtract_vectors(fc1, fc0))
                    uv = self.edge_vector(u, v)
                    nuv = self.edge_normal((u, v), area_weight=True)
                    t = cross_vectors(uv, nuv)
                    if dot_vectors(t, n01) < 0:
                        t = scale_vector(t, -1)

                    interface_angles[(u, v)] = acos(dot_vectors(z0, t))

        return interface_angles
