import sys
import matplotlib.pyplot as plt
import numpy as np
import numpy.linalg as LA
import json
from os.path import basename
import argparse
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import itertools

cfg={
        'other_ms': 6,
        'bond_lw': 0.5,
        'defect_ms': 12,     # bright red defect spheres
        'bond_alpha': 0.5,
        'marker_alpha': 0.2,
        'defect_alpha': 0.9
    }


def wrap(x, A):
    '''
    Wraps point 'x' into the cell whose columns are the vectors of 'A'.
    '''
    Ad = np.array(A, dtype=np.float64)
    b = np.mod(np.linalg.solve(Ad, np.asarray(x, dtype=np.float64)), 1)
    return Ad @ b


def plot_unitcell(data):
    A = np.array(data['cell_vectors'])
    aLinv = np.array(data['primitive_cell_vectors'])

    for f in get_faces(A):
        ax.plot(*np.array(f, dtype=np.float64).T, color='k')

    for f in get_faces(aLinv):
        ax.plot(*np.array(f, dtype=np.float64).T, color='green')


def get_faces(M):
    faces = []

    for i in range(3):
        vec_a = M[:, i]
        vec_b = M[:, (i+1) % 3]
        vec_c = M[:, (i+2) % 3]
        f = [[0, 0, 0], vec_a, vec_a + vec_b, vec_b, [0, 0, 0]]
        f = [np.array(x, dtype=np.int64) for x in f]
        faces.append(f)
        faces.append([x + vec_c for x in f])
    return faces


def unwrap(dx, A):
    '''
    Searches for the unitcell-equivalent point that makes 'dx' the shortest vector.
    '''
    candidate_DX = dx
    for idx in itertools.product((-1, 0, 1), (-1, 0, 1), (-1, 0, 1)):
        tmp = dx + A@idx
        if LA.norm(tmp) < LA.norm(candidate_DX):
            candidate_DX = tmp
    return candidate_DX

def link2startstop(x):
    # expects 'x' to be in link format, i.e.
    # x = {'position', 'boundary': [ [r0, m0], [r1, m1])]}

    if len(x["boundary"]) > 2:
        print("Malformed boundary: link should have two or zero ends")
    r0, m0 = x["boundary"][0]
    r1, m1 = x["boundary"][1]
    assert (m0 + m1 == 0), "Malformed link: multipliers do not sum to 0"
    if m0 == -1:
        # swap
        m0, m1 = m1, m0
        r0, r1 = r1, r0

    # it should now be possible to deduce the correct wrapping of p0, p1,
    # such that all are in same cell
    return np.array(r0), np.array(r1)


def link_endpoints(x, A):
    '''
    Wraps the link centroid and both boundary points into the specified cell
    'A', then returns the centroid together with the shortest centroid-relative
    vectors to each endpoint.
    '''
    r0, r1 = link2startstop(x)
    pos = wrap(x['pos'], A)
    dx0 = unwrap(wrap(r0, A) - pos, A)
    dx1 = unwrap(wrap(r1, A) - pos, A)
    return pos, dx0, dx1


def plot_directed_link(x, A):
    pos, dx0, dx1 = link_endpoints(x, A)
    ax.quiver(*(pos + dx0), *(-dx0), color='k')
    ax.quiver(*(pos + dx1), *(-dx1), color='k', arrow_length_ratio=0)


def plot_undirected_link(x, A):
    pos, dx0, dx1 = link_endpoints(x, A)
    ax.quiver(*(pos + dx0), *(-dx0), color='k', arrow_length_ratio=0)
    ax.quiver(*(pos + dx1), *(-dx1), color='k', arrow_length_ratio=0)


def plot_idx(x, i, pos=None):
    if pos is None:
        pos = x['pos']
    ax.text(*pos, "%d" % i)

def plot_pos(x, pos=None):
    if pos is None:
        pos = x['pos']
    ax.text(*pos, "%d %d %d" % tuple(l for l in x['pos']))

def plot_points(data, args):
    point_data = data['points']

    if point_data is None:
        print("No points in latfile.")
        return

    A = np.array(data['cell_vectors'])

    xyz = []
    for i, x in enumerate(point_data):
        if x['coboundary'] is None:
            continue

        pos = wrap(x['pos'], A)
        xyz.append(pos)

        # The point's coboundary lists the links (spin sites) at the corners of
        # this tetrahedron. Draw the pyrochlore bonds as thin black lines
        # between every pair of those sites.
        sites = [pos + unwrap( wrap(link_pos, A) - pos, A)
                 for link_pos, _ in x['coboundary']]
        for a, b in itertools.combinations(sites, 2):
            ax.plot(*np.array([a, b]).T, color=('k', cfg['bond_alpha']),
                                   linewidth=cfg["bond_lw"]
                    )

        if args.show_idx:
            plot_idx(x, i, pos)
        if args.show_pos:
            plot_pos(x, pos)

    # ax.scatter(*np.array(xyz).T, color='r', marker='o')


def plot_links(data, args):
    A = np.array(data['cell_vectors'])
    link_data = data['links']

    if link_data is None:
        print("No links in latfile.")
        return

    # In the literal pyrochlore, the lattice sites sit at the link centres, so
    # render each link as a black sphere rather than an edge between endpoints.
    xyz = []
    for i, x in enumerate(link_data):
        pos = wrap(x['pos'], A)
        xyz.append(pos)

        if args.show_idx:
            plot_idx(x, i, pos)
        if args.show_pos:
            plot_pos(x, pos)

    ax.scatter(*np.array(xyz).T, color='k', marker='o', s=cfg["other_ms"],
               alpha=cfg["marker_alpha"], depthshade=False)


def find_link(linkpos, link_data):
    for link in link_data:
        if link['pos'] == linkpos:
            return link
    raise LookupError(f"No link at {linkpos}")


def plot_plaqs(data, args):
    link_data = data["links"]
    plaq_data = data["plaqs"]
    
    if plaq_data is None:
        print("No plaquettes in latfile.")
        return
    A = np.array(data['cell_vectors'])

    for i, x in enumerate(plaq_data):
        # re-wrap the plaquette centroid into the specified cell
        pos_w = wrap(x['pos'], A)

        # The pyrochlore plaquette is the polygon spanning the centres of its
        # boundary links (which are the pyrochlore sites), not the diamond
        # link endpoints.
        verts = []
        for link_pos, m in x['boundary']:
            link = find_link(link_pos, link_data)
            # wrap the link centre into the cell, then take the shortest
            # centroid-relative vector
            verts.append(unwrap(wrap(link['pos'], A) - pos_w, A))
        verts = np.array(verts)

        # The boundary links arrive in an arbitrary order, which produces a
        # self-intersecting polygon and render artifacts. Sort the vertices by
        # angle within the plaquette plane, then fan-triangulate from the
        # (interior) centroid so the fill is always a clean set of triangles.
        normal = LA.svd(verts)[2][2]
        e1 = verts[0] - np.dot(verts[0], normal) * normal
        e1 /= LA.norm(e1)
        e2 = np.cross(normal, e1)
        verts = verts[np.argsort(np.arctan2(verts @ e2, verts @ e1))]

        n = len(verts)
        triangles = [[pos_w, pos_w + verts[j], pos_w + verts[(j + 1) % n]]
                     for j in range(n)]
        fill = Poly3DCollection(triangles, alpha=0.35)
        fill.set_facecolor('purple')  # Set surface color
        ax.add_collection3d(fill)

        # Outline just the hexagon boundary (no interior spokes).
        # outline = Poly3DCollection([pos_w + verts], facecolor='none',
        #                            edgecolor=('k', cfg['bond_alpha']),
        #                            linewidth=cfg["bond_lw"])
        # ax.add_collection3d(outline)

        if args.show_idx:
            plot_idx(x, i, pos_w)
        if args.show_pos:
            plot_pos(x, pos_w)


def plot_vols(data, args):
    vol_data = data['vols']
    if vol_data is None:
        print("No volumes in latfile.")
        return
    A = np.array(data['cell_vectors'])
    xyz = []
    for i, x in enumerate(vol_data):
        pos = wrap(x['pos'], A)
        xyz.append(pos)

        if args.show_idx:
            plot_idx(x, i, pos)
        if args.show_pos:
            plot_pos(x, pos)

    ax.scatter(*np.array(xyz).T, color='b', marker='o')


def plot_deleted(data, args):
    # Render the spins removed in the very first dilution stage as bright red
    # spheres. These are stored separately from the surviving links (which
    # plot_links draws as black pyrochlore sites).
    A = np.array(data['cell_vectors'])
    locs = data.get('deleted_spin_locs')

    if not locs:
        print("No deleted_spin_locs in latfile (needs --save_lattice with C++ "
              "recording of the initial dilution).")
        return

    xyz = [wrap(p, A) for p in locs]
    ax.scatter(*np.array(xyz).T, color='red', marker='o',
               s=cfg["defect_ms"] ,
               alpha=cfg["defect_alpha"], depthshade=False)


func_to_run = {
        'points': plot_points,
        'links': plot_links,
        'plaqs': plot_plaqs,
        'vols': plot_vols,
        'deleted': plot_deleted
        }

ap = argparse.ArgumentParser()
ap.add_argument("file", help="The .lat.json file specifying the lattice",
                type=str)
ap.add_argument("objects", nargs='+', choices=func_to_run.keys())
ap.add_argument("--show_idx", action='store_true')
ap.add_argument("--show_pos", action='store_true')
ap.add_argument("--undirected", action='store_true')
ap.add_argument("--save")

args = ap.parse_args()

data = None
with open(args.file, 'r') as f:
    data = json.load(f)

fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.axis('off')

for arg in args.objects:
    func_to_run[arg](data, args)

def parse_filename(fname):
    retval = {}
    for tok in basename(fname).split(';'):
        t = tok.split('=')
        if len(t) != 2:
            continue
        retval[t[0]] = t[1]
    return retval

opts = parse_filename(args.file)
try:
    ax.set_title(f'{float(opts['p'])*100}% dilution, removed {opts['nn']}-neighbours')
except KeyError:
    ax.set_title(basename(args.file))

if args.save:
    ax.set_title('')
    print("Saving to "+args.save)
    fig.savefig(args.save,dpi=300,transparent=True)
else:
    plt.show()
