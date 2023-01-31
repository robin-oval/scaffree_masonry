from __future__ import absolute_import
from __future__ import division
from __future__ import print_function


__all__ = []

def masonry_pattern_to_vrml(mesh, bkeys, path):
    """Export block data from masonry pattern to VRML file.

    Inputs
    ------
    mesh : MasonryPattern
        The MasonryPattern object.
    bkeys : list
        The keys of the blocks to export.
    path : str
        The VRML file path
    """

    with open(path, 'w') as file:

        file.write(
            '#VRML V2.0 utf8 \n Background { \n \t skyColor [ 0.627451 0.627451 0.627451 ] \n }')

        for bkey in bkeys:
            block = mesh.block_volume(bkey)

            # start
            file.write('\n Shape { # \n \t appearance Appearance { \n \t \t material    Material { \n \t \t ambientIntensity 0 \n \t \t diffuseColor     0 1 0 \n \t \t specularColor    0 0 0 \n \t \t emissiveColor    0 0 0 \n \t \t shininess        1 \n \t \t transparency     0 \n \t } \n } # appearance')
            file.write(
                '\n geometry IndexedFaceSet { # triangle mesh \n \t ccw TRUE \n \t convex TRUE \n \t solid  FALSE')

            # coordIndex
            file.write('\n \t coordIndex [ ')
            for i, fkey in enumerate(block.faces()):
                file.write('\n \t \t')
                for vkey in block.face_vertices(fkey):
                    file.write(' ' + str(vkey) + ',')
                file.write(' -1,')
                if i == 0:
                    file.write(' # triangle 0')
            file.write('\n \t ] # {} triangles'.format(
                block.number_of_faces()))

            # texCoordIndex
            file.write('\n \t texCoordIndex [ ')
            for i, fkey in enumerate(block.faces()):
                file.write('\n \t \t')
                for vkey in block.face_vertices(fkey):
                    file.write(' ' + str(vkey) + ',')
                file.write(' -1,')
                if i == 0:
                    file.write(' # triangle 0')
            file.write('\n \t ] # {} triangles'.format(
                block.number_of_faces()))

            # coord
            file.write('\n \t coord Coordinate { point [ ')
            for i, vkey in enumerate(block.vertices()):
                file.write('\n \t \t')
                for xyz in block.vertex_coordinates(vkey):
                    file.write(' ' + str(xyz))
                file.write(',')
                if i == 0:
                    file.write(' # coord point 0')
            file.write('\n \t ] } # ' +
                       str(block.number_of_vertices()) + ' coord points')

            # normal
            file.write('\n \t normal Normal { vector [ ')
            for i, fkey in enumerate(block.faces()):
                file.write('\n \t \t')
                for xyz in block.face_normal(fkey):
                    file.write(' ' + str(xyz))
                file.write(',')
                if i == 0:
                    file.write(' # normal vector 0')
            file.write('\n \t ] } # ' +
                       str(block.number_of_faces()) + ' normal vectors')

            # end
            file.write('\n \t } # geometry \n } #')


def herringbone_to_vrml(mesh, path, from_step=None, to_step=None, from_substep=None, to_substep=None):
    """Export data from herringbone pattern to VRML file.

    Inputs
    ------
    path : str
        VRML file path
    from_step : int
        Initial assembly step from which to indluce the blocks.
    to_step : int
        Final assembly step to which to indluce the blocks.
    """

    from_step = from_step if from_step else float('-inf')
    to_step = to_step if to_step else float('inf')

    from_substep = from_substep if from_substep else float('-inf')
    to_substep = to_substep if to_substep else float('inf')

    count = 0

    with open(path, 'w') as file:

        file.write(
            '#VRML V2.0 utf8 \n Background { \n \t skyColor [ 0.627451 0.627451 0.627451 ] \n }')

        for bkey in mesh.blocks():
            step, substep = mesh.block_step(bkey), mesh.block_substep(bkey)
            if step >= from_step and step <= to_step:
                if substep >= from_substep and substep <= to_substep:
                    count += 1
                    block = mesh.block_volume(bkey)

                    # start
                    file.write('\n Shape { # \n \t appearance Appearance { \n \t \t material 	Material { \n \t \t ambientIntensity 0 \n \t \t diffuseColor     0 1 0 \n \t \t specularColor    0 0 0 \n \t \t emissiveColor    0 0 0 \n \t \t shininess        1 \n \t \t transparency     0 \n \t } \n } # appearance')
                    file.write(
                        '\n geometry IndexedFaceSet { # triangle mesh \n \t ccw TRUE \n \t convex TRUE \n \t solid  FALSE')

                    # coordIndex
                    file.write('\n \t coordIndex [ ')
                    for i, fkey in enumerate(block.faces()):
                        file.write('\n \t \t')
                        for vkey in block.face_vertices(fkey):
                            file.write(' ' + str(vkey) + ',')
                        file.write(' -1,')
                        if i == 0:
                            file.write(' # triangle 0')
                    file.write('\n \t ] # {} triangles'.format(
                        block.number_of_faces()))

                    # texCoordIndex
                    file.write('\n \t texCoordIndex [ ')
                    for i, fkey in enumerate(block.faces()):
                        file.write('\n \t \t')
                        for vkey in block.face_vertices(fkey):
                            file.write(' ' + str(vkey) + ',')
                        file.write(' -1,')
                        if i == 0:
                            file.write(' # triangle 0')
                    file.write('\n \t ] # {} triangles'.format(
                        block.number_of_faces()))

                    # coord
                    file.write('\n \t coord Coordinate { point [ ')
                    for i, vkey in enumerate(block.vertices()):
                        file.write('\n \t \t')
                        for xyz in block.vertex_coordinates(vkey):
                            file.write(' ' + str(xyz))
                        file.write(',')
                        if i == 0:
                            file.write(' # coord point 0')
                    file.write('\n \t ] } # ' +
                               str(block.number_of_vertices()) + ' coord points')

                    # normal
                    file.write('\n \t normal Normal { vector [ ')
                    for i, fkey in enumerate(block.faces()):
                        file.write('\n \t \t')
                        for xyz in block.face_normal(fkey):
                            file.write(' ' + str(xyz))
                        file.write(',')
                        if i == 0:
                            file.write(' # normal vector 0')
                    file.write('\n \t ] } # ' +
                               str(block.number_of_faces()) + ' normal vectors')

                    # end
                    file.write('\n \t } # geometry \n } #')

    print('{} blocks exported'.format(count))

def herringbone_to_3ddat(mesh, path, from_step, to_step):
    """Export data from herringbone pattern to 3ddat file.

    Inputs
    ------
    path : str
        3ddat file path
    from_step : int
        Initial assembly step from which to indluce the blocks.
    to_step : int
        Final assembly step to which to indluce the blocks.
    """

    with open(path, 'w') as file:

        block_count = 0
        face_count = 0

        for bkey in mesh.blocks():
            if mesh.block_step(bkey) >= from_step and mesh.block_step(bkey) <= to_step:
                block_count += 1
                block = mesh.block_volume(bkey)

                file.write('; block {} \n'.format(block_count))
                file.write(
                    'poly reg \t \t {} mat 1 con 1 & \n'.format(block_count))

                for fkey in block.faces():
                    face_count += 1
                    file.write('face ID \t \t {}  '.format(face_count))
                    for i, vkey in enumerate(reversed(block.face_vertices(fkey))):
                        if i != 0:
                            file.write('\t \t \t \t \t')
                        for xyz in block.vertex_coordinates(vkey):
                            file.write('{:.6f} '.format(xyz))
                        file.write('& \n')

        file.write('ret')
