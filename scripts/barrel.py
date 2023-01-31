import os

from compas.geometry import discrete_coons_patch
from compas.geometry import dot_vectors

from scaffree_masonry.utilities import parabola_arc_length

from compas_quad.datastructures import QuadMesh

from scaffree_masonry.datastructure import MasonryPattern

from scaffree_masonry.interoperability import masonry_pattern_to_vrml

from compas.utilities import pairwise

from compas_view2.app import App

# parameters

span, height, length, offset = 10.0, 4.0, 10.0, 1.0
lox_thick, thick = 0.25, 0.25
brick_length, brick_width = 0.125, 0.0625

n_lox = 6
right, left = True, True
toggle_orientation, toggle_offset = 1, 0

HERE = os.path.dirname(__file__)
export_3dec = False
export_json = True
filename = 'tmp'

view = True
view_percentage = 0.5

# subdivision

n_edges_span_target = int(parabola_arc_length(height, span) / brick_length)
spacing = int((n_edges_span_target - n_lox + 1) / n_lox)
spacing = 2 * int(spacing / 2) + 1
d1 = spacing * n_lox + n_lox - 1
d2 = int(parabola_arc_length(offset, length) / brick_width)
print('Mesh subdivision with {} arches divided by {}'.format(d2, d1))

# geometry

ab, bc, dc, ad = [], [], [], []
for i in range(d1 + 1):
    x = (i / d1 - 0.5) * span
    y = - length / 2
    z = - 4 * height / span ** 2 * (x - span / 2) * (x + span / 2)
    ab.append([x, y, z])
    y = length / 2
    dc.append([x, y, z])
for j in range(d2 + 1):
    y = (j / d2 - 0.5) * length
    z = 0.0
    x = - 4 * offset / length ** 2 * (y - length / 2) * (y + length / 2) - span /2
    ad.append([x, y, z])
    x *= -1
    bc.append([x, y, z])
cd = list(reversed(dc))
da = list(reversed(ad))
vertices, faces = discrete_coons_patch(ab, bc, dc, ad)
mesh = MasonryPattern.from_vertices_and_faces(vertices, faces)
if dot_vectors(mesh.face_normal(mesh.get_any_face()), [0.0, 0.0, 1.0]) < 0:
    mesh.flip_cycles()
mesh.collect_polyedges()
edgeperpkey = {}
for pkey, polyedge in mesh.polyedges(data=True):
    uv0 = polyedge[:2]
    edgeperpkey[(min(uv0), max(uv0))] = pkey

# supports

support_polyedges = []
for pkey, polyedge in mesh.polyedges(data=True):
    if mesh.is_edge_on_boundary(*polyedge[:2]):
        if abs(mesh.vertex_coordinates(polyedge[1])[2]) < 1e-6:
            if mesh.halfedge[polyedge[0]][polyedge[1]] is None:
                polyedge = list(reversed(polyedge))
            support_polyedges.append(polyedge)
support_polyedge = support_polyedges[0]

u = support_polyedge[int(d2/2)]
v = [vkey for vkey in mesh.vertex_neighbors(u) if not mesh.is_vertex_on_boundary(vkey)][0]
spine = mesh.collect_polyedge(u, v)
pkey0 = None
for edge in [spine[:2], spine[-2:]]:
    u, v = min(edge), max(edge)
    if (u, v) in edgeperpkey:
        pkey0 = edgeperpkey[(u, v)]
        break
edge0 = spine[:2]

# tessellation

lox_edges0 = mesh.loxodrome_patch(spine, n_lox, directions=0)
lox_edges1 = mesh.loxodrome_patch(list(reversed(spine)), n_lox, directions=0)
lox_edges = set()
for u, v in lox_edges0:
    lox_edges.add((min(u, v), max(u, v)))
for u, v in lox_edges1:
    lox_edges.add((min(u, v), max(u, v)))
faces_to_skip = set([fkey for edge in lox_edges for fkey in mesh.edge_faces(*edge)])
lox_bkeys = set([mesh.face_block(fkey) for fkey in faces_to_skip if fkey is not None])
mesh.staggered_blocks(edge0, flip=True, offset=False, faces_to_skip=faces_to_skip)
mesh.dummy_blocks()

# voussoirs

block2thickness = {bkey: lox_thick if bkey in lox_bkeys else thick for bkey in mesh.blocks()}
block2ecc = {bkey: block2thickness[bkey] / 2 for bkey in mesh.blocks()}
block2volume = mesh.block_volumes(block2thickness, triangulate=True, block2ecc=block2ecc)

# assembly

supported_vertices = [vkey for polyedge in support_polyedges for vkey in polyedge]
mesh.assembly_type_barrel(supported_vertices, pkey0, include_substeps=True)
max_step = int(max([mesh.block_step(bkey) for bkey in mesh.blocks()]))
step2substeps = {step: [] for step in range(max_step + 1)}
for bkey in mesh.blocks():
    step, substep = mesh.block_step(bkey), mesh.block_substep(bkey)
    step2substeps[step].append(substep)
step2maxsubstep = {}
for step, substeps in step2substeps.items():
    step2maxsubstep[step] = max(substeps)
print(step2maxsubstep)

# export

if export_json:
    PATH = os.path.join(HERE, 'barrel_{}.json'.format(filename))
    mesh.to_json(PATH)
    print('Exported design {} to {}'.format(mesh, PATH))

if export_3dec:
                
    step_export = 10
    nb_steps = mesh.number_of_assembly_steps()
    step_splits = [step for step in range(nb_steps) if step == 0 or step == nb_steps - 1 or step % step_export == 0]    

    data = {}

    for step0, step1 in pairwise(step_splits):
        
        from_step = step0
        to_step = step1 - 1
        data[(from_step, to_step)] = []
        for bkey in mesh.blocks():
            step = mesh.block_step(bkey)
            if from_step <= step and step <= to_step:
                data[(from_step, to_step)].append(bkey)

        if step1 == 2 * k or step1 == 3 * k:
            data[(step1, step1)] = {}
            for bkey in mesh.blocks():
                step, substep = mesh.block_step(bkey), mesh.block_substep(bkey)
                if step == step1:
                    if substep in data[(step1, step1)]:
                        data[(step1, step1)][substep].append(bkey)
                    else:
                        data[(step1, step1)][substep] = [bkey]


    for (from_step, to_step), subdata in data.items():

        if isinstance(subdata, list):
            bkeys = subdata
            parameters = [int(radius_out), int(radius_in), int(phi_deg), int(n), int(right), int(left), int(from_step), int(to_step)]
            name = ''
            for param in parameters:
                name += '_' + str(param)
            FILE_3DEC = os.path.join(HERE, 'threedec/dome_herringbone{}.wrl'.format(name))

            masonry_pattern_to_vrml(mesh, bkeys, FILE_3DEC)
            print('Exported {} blocks from step {} to step {} as 3DEC model at {}'.format(len(bkeys), from_step, to_step, FILE_3DEC))

        elif isinstance(subdata, dict):
            for substep, bkeys in subdata.items():
                parameters = [int(radius_out), int(radius_in), int(phi_deg), int(n), int(right), int(left), int(from_step), int(to_step), int(substep)]
                name = ''
                for param in parameters:
                    name += '_' + str(param)
                FILE_3DEC = os.path.join(HERE, 'threedec/dome_herringbone{}.wrl'.format(name))

                masonry_pattern_to_vrml(mesh, bkeys, FILE_3DEC)
                print('Exported {} blocks at step {} and substep {} as 3DEC model at {}'.format(len(bkeys), from_step, substep, FILE_3DEC))


# view

print(mesh.number_of_faces(), len(set([fkey for bkey in mesh.blocks() for fkey in mesh.block_faces(bkey)])))

if view:
    pbm_bkeys = set([bkey for bkey in mesh.blocks() if mesh.block_substep(bkey) > 200])
    viewer = App()
    viewer.add(mesh)
    for bkey in mesh.blocks():
        step = mesh.block_step(bkey)
        substep = mesh.block_substep(bkey)
        substeps = sorted(step2substeps[step])
        max_substep = substeps[-2]
        if step < view_percentage * max_step:
            color = (1.0, 0.0, 1.0) if bkey in lox_bkeys else (min(substep / max(1.0, max_substep), 1.0), step / max_step, 0.0)
            # color = (1.0, 0.0, 1.0) if bkey in pbm_bkeys else (0.0, 0.0, 0.0)
            try:
                viewer.add(mesh.block_volume(bkey), facecolor=color)
            except:
                pass
    viewer.show()
