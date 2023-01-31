import os

from numpy import mean, std

from compas_quad.datastructures import CoarseQuadMesh

from scaffree_masonry.datastructure import MasonryPattern

from compas.datastructures import mesh_weld, mesh_disconnected_faces

from compas_quad.coloring import quad_mesh_polyedge_2_coloring

from scaffree_masonry.loxodrome import quadmesh_edge_loxodrome_direction, quadmesh_edge_front_right_edge, quadmesh_edge_front_left_edge, quadmesh_edge_back_right_edge, quadmesh_edge_back_left_edge
from scaffree_masonry.staggered import quadmesh_staggered_edges

from compas.geometry import distance_point_point, dot_vectors
from compas.geometry import oriented_bounding_box_numpy

from compas_view2.app import App

from compas.utilities import pairwise


brick_width = .12 # [m]
lox_thickness, brick_thickness = .03, .03 # [m]

n_lox = 2 # [-] per spine branch

dual = True

triangulate_blocks = False

export_json = True

view = True
view_blocks = True
view_courses = False

HERE = os.path.dirname(__file__)
FILENAME = "tripod_3d_bf_a.json"
PATH = os.path.join(HERE, FILENAME)

mesh = CoarseQuadMesh.from_json(PATH)
mesh = mesh_weld(mesh)
mesh.collect_strips()
mesh.set_strips_density_target(brick_width)
mesh.densification()
v, f = mesh.dense_mesh().to_vertices_and_faces()
mesh = MasonryPattern.from_vertices_and_faces(v, f)
if dot_vectors([0.0, 0.0, 1.0], mesh.face_normal(mesh.get_any_face())) < 0:
    mesh.flip_cycles()
if dual:
    mesh = mesh.dual()
mesh.collect_polyedges()
print(mesh)

supported_pkeys = set()
for pkey, polyedge in mesh.polyedges(data=True):
    if all([abs(mesh.vertex_coordinates(vkey)[2]) < 0.1 for vkey in polyedge]):
        supported_pkeys.add(pkey)
print('supported_pkeys', supported_pkeys)

supported_vkeys = set([vkey for pkey in supported_pkeys for vkey in mesh.polyedge_vertices(pkey)])

spine_pkeys = set()
spine_polyedges = []
if not dual:
    for pkey, polyedge in mesh.polyedges(data=True):
        start, end  = polyedge[0], polyedge[-1]
        if (start in supported_vkeys and mesh.vertex_degree(end) == 6) or (end in supported_vkeys and mesh.vertex_degree(start) == 6):
            spine_polyedges.append(polyedge)
            spine_pkeys.add(pkey)
else:
    for pkey, polyedge in mesh.polyedges(data=True):
        for polyedge in [polyedge, list(reversed(polyedge))]:
            if polyedge[0] in supported_vkeys:
                for u, v in pairwise(polyedge):
                    if mesh.halfedge[u][v] is not None and len(mesh.face_vertices(mesh.halfedge[u][v])) == 6:
                        spine_polyedges.append(polyedge)
                        spine_pkeys.add(pkey)

# tessellation

all_lox_edges = []
if not dual:
    for polyedge in spine_polyedges[:1]:
        for polyedge in [polyedge, list(reversed(polyedge))][:1]:
            lox_edges += mesh.loxodrome_patch(polyedge, n_lox, directions=0)
else:
    for polyedge in spine_polyedges:
        edge = polyedge[:2]
        edge = mesh.halfedge_after(*mesh.halfedge_after(*edge))
        polyedge = mesh.collect_polyedge(*edge)
        n = len(polyedge)
        toggle = 0
        for i, edge in enumerate(pairwise(polyedge)):
            if mesh.face_degree(mesh.halfedge_face(*list(reversed(edge)))) == 6:
                toggle = i
                continue
            if toggle:
                if (i - toggle) % 12 == 10:
                    lox_edges = quadmesh_edge_loxodrome_direction(mesh, edge, quadmesh_edge_front_left_edge)[1:]
                    mesh.create_blocks(lox_edges, name='lox')
                    all_lox_edges.append(lox_edges)
        toggle = 0
        for i, edge in enumerate(pairwise(list(reversed(polyedge)))):
            if mesh.face_degree(mesh.halfedge_face(*edge)) == 6:
                toggle = i
                continue
            if toggle:
                if (i - toggle) % 12 == 10:
                    lox_edges = quadmesh_edge_loxodrome_direction(mesh, edge, quadmesh_edge_back_left_edge)[1:]
                    mesh.create_blocks(lox_edges, name='lox')
                    all_lox_edges.append(lox_edges)
faces_to_skip = set([fkey for lox_edges in all_lox_edges for edge in lox_edges for fkey in mesh.edge_faces(*edge)])
lox_bkeys = set([mesh.face_block(fkey) for fkey in faces_to_skip if fkey is not None])


sing = [fkey for fkey in mesh.faces() if mesh.face_degree(fkey) == 6][0]
spine_fkeys = set([sing])
for edge in mesh.face_halfedges(sing):
    edge = list(reversed(edge))
    edges = mesh.collect_strip(*edge, both_sides=False)
    if edges[-1][-1] in supported_vkeys:
        spine_fkeys.update([mesh.halfedge_face(*edge) for edge in edges[:-1]])

mesh_disc = mesh.copy()
for fkey in spine_fkeys:
    mesh_disc.delete_face(fkey)

for polyedge in spine_polyedges:
    stagg_edges = quadmesh_staggered_edges(mesh_disc, polyedge[:2], flip=True, offset=True)
    mesh.create_blocks(stagg_edges, faces_to_skip=faces_to_skip, name='stagg')

mesh.dummy_blocks()

# voussoirs

block2thickness = {bkey: lox_thickness if mesh.block_type(bkey) == 'lox' else brick_thickness for bkey in mesh.blocks()}
block2ecc = {bkey: - block2thickness[bkey] / 2 for bkey in mesh.blocks()}
block2volume = mesh.block_volumes(block2thickness, triangulate=triangulate_blocks, block2ecc=block2ecc)

# assembly

mesh.assembly_polyedges_twocolor(spine_pkeys)
max_step = int(max([mesh.block_step(bkey) for bkey in mesh.blocks()]))
print('max step', max_step)

# results

block_lengths, block_widths, block_thicknesses = [], [], []
for bkey in mesh.blocks():
    if mesh.block_type(bkey) == 'lox':
        points = [mesh.vertex_coordinates(vkey) for fkey in mesh.block_faces(bkey) for vkey in mesh.face_vertices(fkey)]
        bbox = oriented_bounding_box_numpy(points)
        a = distance_point_point(bbox[1], bbox[0])
        b = distance_point_point(bbox[3], bbox[0])
        c = distance_point_point(bbox[4], bbox[0])
        t, w, l = list(sorted([a, b, c]))
        block_lengths.append(l)
        block_widths.append(w)
        block_thicknesses.append(t)

pkey2color = quad_mesh_polyedge_2_coloring(mesh)
color0 = pkey2color[next(iter(spine_pkeys))]
pkeys0 = set([pkey for pkey in mesh.polyedges() if pkey2color[pkey] == color0])

edges1 = [edge for pkey in mesh.polyedges() for edge in mesh.polyedge_edges(pkey) if pkey not in pkeys0]
course_widths = [mesh.edge_length(*edge) for edge in edges1]

for name in ['block_lengths', 'block_widths', 'block_thicknesses', 'course_widths']:
    values = eval(name)
    print("herringbone block {} min {} m, max {} m, mean {} m and standard deviation {} m".format(name, round(min(values), 3), round(max(values), 3), round(mean(values), 3), round(std(values), 3)))


# export

if export_json:
    PATH = os.path.join(HERE, 'tmp.json')
    mesh.to_json(PATH)
    print('Exported design {} to {}'.format(mesh, PATH))

# view

if view:
    viewer = App()

    # viewer.add(mesh)

    if view_blocks:

        for bkey in mesh.blocks():
            step = mesh.block_step(bkey)
            if step <= 1.0 * max_step:
                try:
                    if mesh.block_faces(bkey)[0] in spine_fkeys:
                        facecolor = (0.5, 0.5, 0.5)
                    elif mesh.block_type(bkey) == 'lox':
                        facecolor = (1.0, 0.0, 1.0)
                    else:
                        facecolor = (0.8 * step / max_step, 0.8 * step / max_step, 1.0)
                    viewer.add(mesh.block_volume(bkey), facecolor=facecolor)
                except:
                    print('pbm view block', bkey)
                    pass

    if view_courses:

        from compas.geometry import Line, Polyline

        edges0 = set([tuple(sorted(edge)) for pkey in pkeys0 for edge in mesh.polyedge_edges(pkey)])
        
        # for bkey in mesh.blocks():
        #     if mesh.block_type(bkey) == 'lox':
        #         f0, f1 = mesh.block_faces(bkey)
        #         edge_mid = mesh.face_adjacency_halfedge(f0, f1)
        #         edge_front = mesh.face_opposite_edge(*edge_mid)
        #         edge_back = mesh.face_opposite_edge(*tuple(reversed(edge_mid)))
        #         for edge in [edge_mid, edge_front, edge_back]:
        #             edge = tuple(sorted(edge))
        #             if edge in edges0:
        #                 edges0.remove(edge)
        #         viewer.add(Line(mesh.edge_midpoint(*edge_front), mesh.edge_midpoint(*edge_back)), linecolor=(1.0, 0.0, 1.0))
        
        for lox_edges in all_lox_edges:
            viewer.add(Polyline([mesh.edge_midpoint(*edge) for edge in lox_edges]), linecolor=(1.0, 0.0, 1.0))
            for edge_mid in lox_edges:
                edge_front = mesh.face_opposite_edge(*edge_mid)
                edge_back = mesh.face_opposite_edge(*tuple(reversed(edge_mid)))
                for edge in [edge_mid, edge_front, edge_back]:
                    if edge is not None:
                        edge = tuple(sorted(edge))
                        if edge in edges0:
                            edges0.remove(edge)

        for edge in edges0:
            step = min([mesh.block_step(mesh.face_block(fkey)) for fkey in mesh.edge_faces(*edge) if fkey is not None])
            linecolor = (0.8 * step / max_step, 0.8 * step / max_step, 1.0)
            viewer.add(Line(*mesh.edge_coordinates(*edge)), linecolor=linecolor)

    viewer.show()
