import argparse
from collections import defaultdict
import copy
import re

import numpy as np

import openmc
from openmc.data import get_thermal_name
from openmc.data.neutron import _get_metadata
from .macrobody import Box, RightCircularCylinder as RCC


_CELL1_RE = re.compile(r'\s*(\d+)\s+(\d+)([ \t0-9:#().dDeE\+-]+)((?:\S+\s*=.*)?)')
_CELL_KEYWORDS_RE = re.compile(r'[A-Za-z: \t]+=[-0-9:() \t]+')
_CELL2_RE = re.compile(r'\s*(\d+)\s+like\s+(\d+)\s+but\s+(\S+\s*=.*)')
_SURFACE_RE = re.compile(r'\s*(\*?\d+)(\s*[-0-9]+)?\s+(\S+)((?:\s+\S+)+)')
_MATERIAL_RE = re.compile(r'\s*[Mm](\d+)((?:\s+\S+)+)')
_SAB_RE = re.compile(r'\s*[Mm][Tt](\d+)((?:\s+\S+)+)')
_MODE_RE = re.compile(r'\s*mode(?:\s+\S+)*')
_COMPLEMENT_RE = re.compile(r'(#)(\d+)')
_REPEAT_RE = re.compile(r'(\d+)\s+(\d+)[rR]')
_NUM_RE = re.compile(r'(\d)([+-])(\d)')


def float_(val):
    return float(_NUM_RE.sub(r'\1e\2\3', val))


def parse_cell(line):
    if 'like' in line.lower():
        raise NotImplementedError('like N but form not supported')
        # Handle LIKE n BUT form
        m = _CELL2_RE.match(line.lower())
        if m is not None:
            g = m.groups()
            parameters = (dict(kv.split('=') for kv in g[-1].split())
                          if len(g) == 3 else {})
            return {'id': g[0], 'like': g[1], 'parameters': parameters}

    else:
        # Handle normal form
        m = _CELL1_RE.match(line)
        if m is not None:
            g = m.groups()
            if g[1] == '0':
                density = None
                region = g[2].strip()
            else:
                words = g[2].split()
                density = float_(words[0])
                region = ' '.join(words[1:])
            parameters = {}
            if g[-1]:
                words = (g[-1] + ' after').split('=')
                for before, after in zip(words[:-1], words[1:]):
                    key = before.split()[-1]
                    value = ' '.join(after.split()[:-1])
                    parameters[key] = value
            return {'id': int(g[0]), 'material': int(g[1]), 'density': density,
                    'region': region, 'parameters': parameters}


def parse_surface(line):
    m = _SURFACE_RE.match(line)
    if m is not None:
        g = m.groups()
        surface = {}
        if '*' in g[0]:
            surface['reflective'] = True
            uid = int(g[0][1:])
        else:
            surface['reflective'] = False
            uid = int(g[0])
        surface.update({'id': uid, 'mnemonic': g[2].lower(),
                        'coefficients': g[3].strip()})
        if g[1] is not None:
            if int(g[1]) < 0:
                surface['periodic'] = int(g[1])
                raise NotImplementedError('Periodic boundary conditions not supported')
            else:
                surface['tr'] = int(g[1])
                raise NotImplementedError('Surface transformation not supported')
        return surface
    else:
        raise NotImplementedError("Unable to convert surface card: {}".format(line))


def parse_data(section):
    data = {'materials': defaultdict(dict)}

    lines = section.split('\n')
    for line in lines:
        if _MATERIAL_RE.match(line):
            g = _MATERIAL_RE.match(line).groups()
            spec = g[1].split()
            try:
                nuclides = zip(spec[::2], map(float_, spec[1::2]))
            except Exception:
                raise ValueError('Invalid material specification?')
            uid = int(g[0])
            data['materials'][uid].update({'id': uid, 'nuclides': nuclides})
        elif _SAB_RE.match(line):
            g = _SAB_RE.match(line).groups()
            uid = int(g[0])
            spec = g[1].split()
            data['materials'][uid]['sab'] = spec
        elif line.lower().startswith('mode'):
            words = line.split()
            data['mode'] = words[1:]
        elif line.lower().startswith('kcode'):
            words = line.split()
            data['kcode'] = {'n_particles': int(words[1]),
                             'initial_k': float_(words[2]),
                             'inactive': int(words[3]),
                             'batches': int(words[4])}

    return data


def parse(filename):
    # Find beginning of cell section
    text = open(filename, 'r').read()
    m = re.search(r'^[ \t]*(\d+)[ \t]+', text, flags=re.MULTILINE)
    text = text[m.start():]

    sections = re.split('\n[ \t]*\n', text)
    for i in range(len(sections)):
        # Remove end-of-line comments
        sections[i] = re.sub('\$.*$', '', sections[i], flags=re.MULTILINE)

        # Remove comment cards
        sections[i] = re.sub('^[ \t]*?[cC].*?$\n?', '', sections[i], flags=re.MULTILINE)

        # Turn continuation lines into single line
        sections[i] = re.sub('&.*\n', ' ', sections[i])
        sections[i] = re.sub('\n {5}', ' ', sections[i])

        # Expand repeated numbers
        m = _REPEAT_RE.search(sections[i])
        while m is not None:
            sections[i] = _REPEAT_RE.sub(' '.join((int(m.group(2)) + 1)*[m.group(1)]),
                                         sections[i], 1)
            m = _REPEAT_RE.search(sections[i])

    cells = list(map(parse_cell, sections[0].strip().split('\n')))
    surfaces = list(map(parse_surface, sections[1].strip().split('\n')))
    data = parse_data(sections[2])

    return cells, surfaces, data


def get_openmc_materials(materials):
    openmc_materials = {}
    for m in materials.values():
        if 'id' not in m:
            continue
        material = openmc.Material(m['id'])
        for nuclide, percent in m['nuclides']:
            zaid, xs = nuclide.split('.')
            name, element, Z, A, metastable = _get_metadata(int(zaid), 'mcnp')
            if percent < 0:
                if A > 0:
                    material.add_nuclide(name, abs(percent), 'wo')
                else:
                    material.add_element(element, abs(percent), 'wo')
            else:
                if A > 0:
                    material.add_nuclide(name, percent, 'ao')
                else:
                    material.add_element(element, percent, 'ao')

        if 'sab' in m:
            for sab in m['sab']:
                name, xs = sab.split('.')
                material.add_s_alpha_beta(get_thermal_name(name))
        openmc_materials[m['id']] = material

    return openmc_materials


def get_openmc_surfaces(surfaces):
    # Ensure that autogenerated IDs for surfaces don't conflict
    openmc.Surface.next_id = max(s['id'] for s in surfaces) + 1

    openmc_surfaces = {}
    for s in surfaces:
        if s['mnemonic'] == 'p':
            coeffs = s['coefficients'].split()
            if len(coeffs) == 9:
                raise NotImplementedError('General plane defined by three points not supported')
            else:
                A, B, C, D = map(float_, coeffs)
            surf = openmc.Plane(surface_id=s['id'], A=A, B=B, C=C, D=D)
        elif s['mnemonic'] == 'px':
            x0 = float_(s['coefficients'])
            surf = openmc.XPlane(surface_id=s['id'], x0=x0)
        elif s['mnemonic'] == 'py':
            y0 = float_(s['coefficients'])
            surf = openmc.YPlane(surface_id=s['id'], y0=y0)
        elif s['mnemonic'] == 'pz':
            z0 = float_(s['coefficients'])
            surf = openmc.ZPlane(surface_id=s['id'], z0=z0)
        elif s['mnemonic'] == 'so':
            R = float_(s['coefficients'])
            surf = openmc.Sphere(surface_id=s['id'], R=R)
        elif s['mnemonic'] == 's':
            x0, y0, z0, R = map(float_, s['coefficients'].split())
            surf = openmc.Sphere(surface_id=s['id'], x0=x0, y0=y0, z0=z0, R=R)
        elif s['mnemonic'] == 'sx':
            x0, R = map(float_, s['coefficients'].split())
            surf = openmc.Sphere(surface_id=s['id'], x0=x0, R=R)
        elif s['mnemonic'] == 'sy':
            y0, R = map(float_, s['coefficients'].split())
            surf = openmc.Sphere(surface_id=s['id'], y0=y0, R=R)
        elif s['mnemonic'] == 'sz':
            z0, R = map(float_, s['coefficients'].split())
            surf = openmc.Sphere(surface_id=s['id'], z0=z0, R=R)
        elif s['mnemonic'] == 'c/x':
            y0, z0, R = map(float_, s['coefficients'].split())
            surf = openmc.XCylinder(surface_id=s['id'], y0=y0, z0=z0, R=R)
        elif s['mnemonic'] == 'c/y':
            x0, z0, R = map(float_, s['coefficients'].split())
            surf = openmc.YCylinder(surface_id=s['id'], x0=x0, z0=z0, R=R)
        elif s['mnemonic'] == 'c/z':
            x0, y0, R = map(float_, s['coefficients'].split())
            surf = openmc.ZCylinder(surface_id=s['id'], x0=x0, y0=y0, R=R)
        elif s['mnemonic'] == 'cx':
            R = float_(s['coefficients'])
            surf = openmc.XCylinder(surface_id=s['id'], R=R)
        elif s['mnemonic'] == 'cy':
            R = float_(s['coefficients'])
            surf = openmc.YCylinder(surface_id=s['id'], R=R)
        elif s['mnemonic'] == 'cz':
            R = float_(s['coefficients'])
            surf = openmc.ZCylinder(surface_id=s['id'], R=R)
        elif s['mnemonic'] in ('k/x', 'k/y', 'k/z'):
            coeffs = [float_(x) for x in s['coefficients'].split()]
            x0, y0, z0, R2 = coeffs[:4]
            if len(coeffs) > 4:
                raise NotImplementedError('One-sheet cone not supported')
            if s['mnemonic'] == 'k/x':
                surf = openmc.XCone(surface_id=s['id'], x0=x0, y0=y0, z0=z0, R2=R2)
            elif s['mnemonic'] == 'k/y':
                surf = openmc.YCone(surface_id=s['id'], x0=x0, y0=y0, z0=z0, R2=R2)
            elif s['mnemonic'] == 'k/z':
                surf = openmc.ZCone(surface_id=s['id'], x0=x0, y0=y0, z0=z0, R2=R2)
        elif s['mnemonic'] in ('kx', 'ky', 'kz'):
            coeffs = [float_(x) for x in s['coefficients'].split()]
            x, R2 = coeffs[:2]
            if len(coeffs) > 2:
                raise NotImplementedError('One-sheet cone not supported')
            if s['mnemonic'] == 'kx':
                surf = openmc.XCone(surface_id=s['id'], x0=x, R2=R2)
            elif s['mnemonic'] == 'ky':
                surf = openmc.YCone(surface_id=s['id'], y0=x, R2=R2)
            elif s['mnemonic'] == 'kz':
                surf = openmc.ZCone(surface_id=s['id'], z0=x, R2=R2)
        elif s['mnemonic'] == 'gq':
            a, b, c, d, e, f, g, h, j, k = map(float_, s['coefficients'].split())
            surf = openmc.Quadric(surface_id=s['id'], a=a, b=b, c=c, d=d, e=e,
                                  f=f, g=g, h=h, j=j, k=k)
        elif s['mnemonic'] == 'rcc':
            vx, vy, vz, hx, hy, hz, r = map(float_, s['coefficients'].split())
            if hx == 0.0 and hy == 0.0:
                surf = RCC((vx, vy, vz), hz, r, axis='z')
            elif hy == 0.0 and hz == 0.0:
                surf = RCC((vx, vy, vz), hx, r, axis='x')
            elif hx == 0.0 and hz == 0.0:
                surf = RCC((vx, vy, vz), hy, r, axis='y')
            else:
                raise notImplementedError('RCC macrobody with non-axis-aligned'
                                          'height vector not supported.')
        elif s['mnemonic'] == 'box':
            coeffs = list(map(float_, s['coefficients'].split()))
            surf = Box(coeffs[:3], coeffs[3:6], coeffs[6:9], coeffs[9:])
        else:
            raise NotImplementedError('Surface type "{}" not supported'
                                      .format(s['mnemonic']))

        if s['reflective']:
            surf.boundary_type = 'reflective'

        openmc_surfaces[s['id']] = surf

    return openmc_surfaces


def get_openmc_universes(cells, surfaces, materials):
    openmc_cells = {}
    universes = {}
    root_universe = openmc.Universe(0)
    universes[0] = root_universe

    # Determine maximum IDs so that autogenerated IDs don't conflict
    openmc.Cell.next_id = max(c['id'] for c in cells) + 1
    all_univ_ids = set()
    for c in cells:
        if 'u' in c['parameters']:
            all_univ_ids.add(abs(int(c['parameters']['u'])))
    if all_univ_ids:
        openmc.Universe.next_id = max(all_univ_ids) + 1

    # Since OpenMC doesn't support MCNP's cell-complement notation, we need
    # to replace it with a normal surface-based complement
    cell_by_id = {c['id']: c for c in cells}
    have_complements = True
    used_in_complement = set()
    while have_complements:
        for c in cells:
            match = _COMPLEMENT_RE.search(c['region'])
            if match:
                other_id = int(match.group(2))
                region = cell_by_id[other_id]['region']
                c['region'] = _COMPLEMENT_RE.sub(
                    r'~({})'.format(region),
                    c['region']
                )
                used_in_complement.add(other_id)
                break
        else:
            have_complements = False

    translate_memo = {}
    for c in cells:
        cell = openmc.Cell(cell_id=c['id'])
        region = c['region'].replace(':', '|')

        # Assign region to cell based on expression
        try:
            cell.region = openmc.Region.from_expression(region, surfaces)
        except Exception:
            raise ValueError('Could not parse region: {}'.format(region))

        # Check for unsupported keywords
        if 'trcl' in c['parameters']:
            if cell.id in used_in_complement:
                raise NotImplementedError(
                    'Cannot transform cell (ID={}) that was used '
                    'in a cell complement.'.format(cell.id)
                )

            trcl = c['parameters']['trcl'].strip()
            if not trcl.startswith('('):
                raise NotImplementedError('TRn card not supported.')

            # Drop parentheses
            trcl = trcl[1:-1].split()
            if len(trcl) > 3:
                raise NotImplementedError('Cell rotation not supported.')
            assert len(trcl) == 3
            vector = tuple(float(x) for x in trcl)
            cell.region = cell.region.translate(vector, translate_memo)

        # Add cell to universes if necessary
        if 'u' in c['parameters']:
            if 'lat' not in c['parameters']:
                # Note: a negative universe indicates that the cell is not
                # truncated by the boundary of a higher level cell.
                uid = abs(int(c['parameters']['u']))
                if uid not in universes:
                    universes[uid] = openmc.Universe(uid)
                universes[uid].add_cell(cell)
        else:
            root_universe.add_cell(cell)

        # Look for vacuum boundary condition
        if isinstance(cell.region, openmc.Union):
            if all([isinstance(n, openmc.Halfspace) for n in cell.region]):
                if 'imp:n' in c['parameters'] and c['parameters']['imp:n'] == '0':
                    for n in cell.region:
                        if n.surface.boundary_type == 'transmission':
                            n.surface.boundary_type = 'vacuum'
                    root_universe.remove_cell(cell)
        elif isinstance(cell.region, openmc.Halfspace):
            if 'imp:n' in c['parameters'] and c['parameters']['imp:n'] == '0':
                if cell.region.surface.boundary_type == 'transmission':
                    cell.region.surface.boundary_type = 'vacuum'
                root_universe.remove_cell(cell)

        # Determine material fill if present -- this is not assigned until later
        # in case it's used in a lattice (need to create an extra universe then)
        if c['material'] > 0:
            mat = materials[c['material']]
            if mat.density is None:
                if c['density'] > 0:
                    mat.set_density('atom/b-cm', c['density'])
                else:
                    mat.set_density('g/cm3', abs(c['density']))
            else:
                if mat.density != abs(c['density']):
                    print("WARNING: Cell {} has same material but with a "
                          "different density than that of another.".format(c['id']))
                    mat = mat.clone()
                    if c['density'] > 0:
                        mat.set_density('atom/b-cm', c['density'])
                    else:
                        mat.set_density('g/cm3', abs(c['density']))

        # Create lattices
        if 'fill' in c['parameters']:
            if 'lat' in c['parameters']:
                # Check what kind of lattice this is
                if int(c['parameters']['lat']) == 2:
                    raise NotImplementedError("Hexagonal lattices not supported")

                # Cell filled with Lattice
                uid = int(c['parameters']['u'])
                if uid not in universes:
                    universes[uid] = openmc.RectLattice(uid)
                lattice = universes[uid]

                # Determine dimensions of single lattice element
                if len(cell.region) < 4:
                    raise NotImplementedError('One-dimensional lattices not supported')
                sides = {'x': [], 'y': [], 'z': []}
                for n in cell.region:
                    if isinstance(n.surface, openmc.XPlane):
                        sides['x'].append(n.surface.x0)
                    elif isinstance(n.surface, openmc.YPlane):
                        sides['y'].append(n.surface.y0)
                    elif isinstance(n.surface, openmc.ZPlane):
                        sides['z'].append(n.surface.z0)
                if not sides['x'] or not sides['y']:
                    raise NotImplementedError('2D lattice with basis other than x-y not supported')

                # MCNP's convention is that across the first surface listed is
                # the (1,0,0) element and across the second surface is the
                # (-1,0,0) element
                if sides['z']:
                    v1, v0 = np.array([sides['x'], sides['y'], sides['z']]).T
                else:
                    v1, v0 = np.array([sides['x'], sides['y']]).T

                pitch = abs(v1 - v0)

                def get_universe(uid):
                    if uid not in universes:
                        universes[uid] = openmc.Universe(uid)
                    return universes[uid]

                # Get extent of lattice
                words = c['parameters']['fill'].split()
                if len(words) == 1:
                    # Infinite lattice
                    lattice.pitch = pitch
                    lattice.lower_left = np.min(np.vstack((v0, v1)), axis=0)
                    lattice.dimension = np.ones(pitch.size, dtype=int)
                    universe = get_universe(int(words[0]))
                    if pitch.size == 2:
                        lattice.universes = [[universe]]
                    else:
                        lattice.universes = [[[universe]]]
                    lattice.outer = universe

                    # TODO: Shift center
                else:
                    pairs = re.findall(r'-?\d+\s*:\s*-?\d+', c['parameters']['fill'])
                    i_colon = c['parameters']['fill'].rfind(':')
                    univ_ids = c['parameters']['fill'][i_colon + 1:].split()[1:]

                    if not pairs:
                        raise ValueError('Cant find lattice specification')

                    xmin, xmax = map(int, pairs[0].split(':'))
                    ymin, ymax = map(int, pairs[1].split(':'))
                    zmin, zmax = map(int, pairs[2].split(':'))
                    assert xmax >= xmin
                    assert ymax >= ymin
                    assert zmax >= zmin

                    if pitch.size == 3:
                        index0 = np.array([xmin, ymin, zmin])
                        index1 = np.array([xmax, ymax, zmax])
                    else:
                        index0 = np.array([xmin, ymin])
                        index1 = np.array([xmax, ymax])
                    shape = index1 - index0 + 1

                    # Determine lower-left corner of lattice
                    corner0 = v0 + index0*(v1 - v0)
                    corner1 = v1 + index1*(v1 - v0)
                    lower_left = np.min(np.vstack((corner0, corner1)), axis=0)

                    lattice.pitch = pitch
                    lattice.lower_left = lower_left
                    lattice.dimension = shape

                    # Universe IDs array as ([z], y, x)
                    univ_ids = np.asarray(univ_ids, dtype=int)
                    univ_ids.shape = shape[::-1]

                    # Depending on the order of the surfaces listed, it may be
                    # necessary to flip some axes
                    if (v1 - v0)[0] < 0.:
                        # lattice positions on x-axis are backwards
                        univ_ids = np.flip(univ_ids, axis=-1)
                    if (v1 - v0)[1] < 0.:
                        # lattice positions on y-axis are backwards
                        univ_ids = np.flip(univ_ids, axis=-2)
                    if sides['z'] and (v1 - v0)[2] < 0.:
                        # lattice positions on z-axis are backwards
                        univ_ids = np.flip(univ_ids, axis=-3)

                    # Check for universe ID same as the ID assigned to the cell
                    # itself -- since OpenMC can't handle this directly, we need
                    # to create an extra cell/universe to fill in the lattice
                    if np.any(univ_ids == uid):
                        c = openmc.Cell(fill=mat)
                        u = openmc.Universe(cells=[c])
                        univ_ids[univ_ids == uid] = u.id

                        # Put it in universes dictionary so that get_universe
                        # works correctly
                        universes[u.id] = u

                    # If center of MCNP lattice element is not (0,0,0), we need
                    # to translate the universe
                    center = np.zeros(3)
                    center[:v0.size] = (v0 + v1)/2
                    if not np.all(center == 0.0):
                        for uid in np.unique(univ_ids):
                            # Create translated universe
                            c = openmc.Cell(fill=get_universe(uid))
                            c.translation = -center
                            u = openmc.Universe(cells=[c])
                            universes[u.id] = u

                            # Replace original universes with translated ones
                            univ_ids[univ_ids == uid] = u.id

                    # Get an array of universes instead of IDs
                    lat_univ = np.vectorize(get_universe)(univ_ids)

                    # Fill universes in OpenMC lattice, reversing y direction
                    lattice.universes = lat_univ[..., ::-1, :]

                cell._lattice = True
            else:
                # Cell filled with universes
                if '(' in c['parameters']['fill']:
                    raise NotImplementedError('fill transformation not supported')
                uid = int(c['parameters']['fill'])
                if uid not in universes:
                    for ci in cells:
                        if 'u' in ci['parameters']:
                            if int(ci['parameters']['u']) == uid:
                                if 'lat' in ci['parameters']:
                                    universes[uid] = openmc.RectLattice(uid)
                                else:
                                    universes[uid] = openmc.Universe(uid)
                                break
                cell.fill = universes[uid]

        elif c['material'] > 0:
            cell.fill = mat

        if not hasattr(cell, '_lattice'):
            openmc_cells[c['id']] = cell

    # Expand shorthand notation
    def replace_complement(region, cells):
        if isinstance(region, (openmc.Intersection, openmc.Union)):
            for n in region:
                replace_complement(n, cells)
        elif isinstance(region, openmc.Complement):
            if isinstance(region.node, openmc.Halfspace):
                region.node = cells[region.node.surface.id].region

    for cell in openmc_cells.values():
        replace_complement(cell.region, openmc_cells)
    return universes
