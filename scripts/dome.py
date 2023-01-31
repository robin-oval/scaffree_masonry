import os

from time import time

from math import floor, pi, cos, sin, asin

from random import random

from compas.geometry import dot_vectors

from compas.datastructures import mesh_weld, meshes_join

from compas_quad.datastructures import QuadMesh

from scaffree_masonry.datastructure import MasonryPattern

from scaffree_masonry.interoperability import masonry_pattern_to_vrml

from compas.utilities import pairwise

from compas_view2.app import App

for n in [8]:#, 12, 20, 32]:
    for phi_deg in [90]:#, 90]:
        for radius_out in [5.0]:#, 6.0, 7.0]:

            # phi_deg = 65 # embrace angle
            phi = phi_deg / 180 * pi
            # radius_out = 5.0 # XY radius of outer boundary of dome (support)
            radius_in = 1.25 # XY radius of inner boundary of dome (oculus)
            radius = radius_out / sin(phi) # radius of sphere

            phi_oculus = asin(radius_in / radius)
            
            hb_thickness, pb_thickness = 0.25, 0.25
            brick_length, brick_width = 0.125, 0.0625
            # brick_length, brick_width = 0.25, 0.15
            
            # n = 12
            d1, d2 = int(2 * pi * radius_out / brick_length), int(radius * (phi - phi_oculus) / brick_width)
            # d1 = int(floor(d1 / n) * n) + 1 # needs to be one above a multiple of n for loxdrome to intersec with one block only
            m = floor((d1 / n - 4) / 2)
            d1 = n * ( 2 * m + 4) + 1
            print('Mesh subdivision with {} meridians and {} parallels'.format(d1, d2))

            right, left = True, True
            toggle_orientation, toggle_offset = 1, 0

            HERE = os.path.dirname(__file__)
            export_3dec = False
            export_json = False

            view = True

            # geometry and mesh
            vertices = []
            for j in range(d2):
                for i in range(d1):
                    phi_polar = 2 * pi * i / (d1 - 1)
                    theta_polar = phi_oculus + j / (d2 - 1) * (phi - phi_oculus)
                    vertices.append([radius * sin(theta_polar) * cos(phi_polar), radius * sin(theta_polar) * sin(phi_polar), radius * cos(theta_polar)])
            faces = [[i + j * d1, i + 1 + j * d1, i + 1 + (j + 1) * d1 , i + (j + 1) * d1] for i in range(d1 - 1) for j in range(d2 - 1)]
            mesh = MasonryPattern.from_vertices_and_faces(vertices, faces)
            mesh = mesh_weld(mesh)
            if dot_vectors(mesh.face_normal(mesh.get_any_face()), [0.0, 0.0, 1.0]) < 0:
                mesh.flip_cycles()
            
            t0 = time()
            mesh.collect_polyedges()
            t1 = time()

            if export_json:
                parameters = [int(radius_out), int(radius_in), int(phi_deg), int(n), int(right), int(left)]
                name = ''
                for param in parameters:
                    name += '_' + str(param)
                FILE_MESH = os.path.join(HERE, 'dome_herringbone_mesh{}.obj'.format(name))
                mesh.to_obj(FILE_MESH)
                print('Exported parameterisation mesh {} to {}'.format(mesh, FILE_MESH))

            # identify supports
            bdry2length = {tuple(bdry): mesh.polyedge_length(bdry) for bdry in mesh.vertices_on_boundaries()}
            supports = list(reversed(max(bdry2length, key=bdry2length.get)))
            edge0 = supports[:2]

            t2 = time()

            # pattern generation
            lox_edges = mesh.loxodrome_patch(supports, n, directions=0)
            faces_to_skip = set([fkey for edge in lox_edges for fkey in mesh.edge_faces(*edge)])
            mesh.staggered_blocks(edge0, flip=True, offset=False, faces_to_skip=faces_to_skip)
            mesh.dummy_blocks()

            t3 = time()

            # block generation
            block2thickness = {bkey: hb_thickness if mesh.block_type(bkey) == 'lox' else pb_thickness for bkey in mesh.blocks()}
            block2ecc = {bkey: block2thickness[bkey] / 2 for bkey in mesh.blocks()}
            block2volume = mesh.block_volumes(block2thickness, triangulate=True)

            t4 = time()

            # assembly generation
            mesh.identify_supports(supports)
            mesh.assembly_type_dome(include_substeps=True)
            max_step = int(max([mesh.block_step(bkey) for bkey in mesh.blocks()]))
            step2maxsubstep = {step: [] for step in range(max_step + 1)}
            for bkey in mesh.blocks():
                step, substep = mesh.block_step(bkey), mesh.block_substep(bkey)
                step2maxsubstep[step].append(substep)
            for step, substeps in step2maxsubstep.items():
                step2maxsubstep[step] = max(substeps)
            print(step2maxsubstep)    

            t5 = time()
            # export
            if export_json:
                parameters = [int(radius_out), int(radius_in), int(phi_deg), int(n), int(right), int(left)]
                name = ''
                for param in parameters:
                    name += '_' + str(param)
                FILE_PATTERN = os.path.join(HERE, 'dome_herringbone_pattern{}.obj'.format(name))
                mesh.to_obj(FILE_PATTERN)
                print('Exported design {} to {}'.format(mesh, FILE_PATTERN))

                # all_block_meshes = meshes_join([mesh.block_volume(bkey) for bkey in mesh.blocks()])
                # FILE_VOLUMES = os.path.join(HERE, 'dome_herringbone_volumes{}.obj'.format(name))
                # all_block_meshes.to_obj(FILE_VOLUMES)
                # print('Exported volumes {} to {}'.format(all_block_meshes, FILE_VOLUMES))

            if export_3dec:
                
                angular_subdivision = 15
                nb_steps = mesh.number_of_assembly_steps()
                k = int(angular_subdivision / phi_deg * nb_steps)
                step_splits = [step for step in range(nb_steps + 2) if step == 0 or step == 1 or step == 2 or step == nb_steps + 1 or step % k == 0]    

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

            t6 = time()

            print(t6 - t5, t5 - t4, t4 - t3, t3 - t2, t2 - t1, t1 - t0)

            try:
                max_step = max([mesh.block_step(bkey) for bkey in mesh.blocks()])
            except:
                max_step = 1

            if view:
                viewer = App()
                viewer.add(mesh)
                for bkey in mesh.blocks():
                    try:
                        step = mesh.block_step(bkey)
                    except:
                        step = 0
                    try:
                        substep = mesh.block_substep(bkey)
                    except:
                        substep = 0
                    color = (1.0, 0.0, 1.0) if mesh.block_type(bkey) == 'lox' else (substep / (step2maxsubstep[step] + 1), step / max_step, 0.0)
                    try:
                        viewer.add(mesh.block_volume(bkey), facecolor=color)
                    except:
                        pass
                viewer.show()

