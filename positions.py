#!/usr/bin/env python3

# Transform haskell output into less sucky json and compute positions.
# (implicit units of C-C bond length)

# NOTE: This thing is a *right mess*.  It started out general and then
#       a bunch of hacks were piled on to make it better for hexagonal.
#       Some of it's translated from *Haskell*.
#       Half of it could probably be deleted and would never be missed

from moire.exact import MoirePattern, find_nicer_cell
from sympy import S, Matrix, sqrt, oo, rot_axis3, atan2, simplify

import functools
import itertools
import sys
import json

class Layout:
    ABC_SCALE = 'a-b-c-n-d-k'
    VOLUME_DISAMBIG = 'volume-disambig'
    all = [ABC_SCALE, VOLUME_DISAMBIG]

class Input:
    HASKELL = 'haskell' # legacy format with beta and scale factor:  [BETA, [A, B, C], {"numerator": N, "denominator": D, "rootNum": R}]
    ABC = 'abc'  # [A, B, C]
    RUST_PAIRS = 'rust-pairs'  # a pair of solutions [[A1, B1, C1], [A2, B2, C2]] that are theta and 60-theta
    all = [RUST_PAIRS, HASKELL, ABC]

def main():
    global PARANOID

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--stream-in',
        action='store_true',
        help="interpret the input file as JSON Lines format rather than standard JSON"
    )
    parser.add_argument(
        '--stream-out',
        action='store_true',
        help="outputting json dicts one by one instead of a list."
        " The generated file is technically invalid json, but jq will gobble it up."
    )
    parser.add_argument('--paranoid', '-P', default=0, action='count')
    parser.add_argument('--carefree', '-C', default=0, action='count')
    parser.add_argument(
        '--min-volume', '-n', default=0, type=int, help='inclusive')
    parser.add_argument(
        '--max-volume', '-N', default=None, type=int, help='exclusive')
    parser.add_argument(
        '--key-layout', '-k', default=Layout.ABC_SCALE, choices=Layout.all)
    parser.add_argument(
        '--input-format', default=Input.RUST_PAIRS, choices=Input.all)
    args = parser.parse_args()

    PARANOID = args.paranoid - args.carefree

    if sys.stdin.isatty:
        print("Reading from standard input...", file=sys.stderr)

    def load_input():
        if args.stream_in:
            for line in sys.stdin:
                if line.strip():
                    yield json.loads(line)
        else:
            yield from json.load(sys.stdin)

    run_sol = functools.partial(main_, min_volume=args.min_volume, max_volume=args.max_volume, key_layout=args.key_layout)

    # beta = 3, r = 1
    add_parts_for_hex = lambda abc: (3, abc, {'numerator': 1, 'denominator': 1, 'rootNum': 1})
    if args.input_format == Input.HASKELL:
        def handle_item(x):
            yield run_sol(x)
    elif args.input_format == Input.ABC:
        def handle_item(abc):
            yield run_sol(add_parts_for_hex(abc))
    elif args.input_format == Input.RUST_PAIRS:
        def handle_item(xs):
            a, b = xs
            volume = a[2] if a[2] % 2 == 1 else a[2] // 2
            if args.min_volume is not None and not (args.min_volume <= volume):
                return
            if args.max_volume is not None and not (volume <= args.max_volume):
                return
            yield run_sol(add_parts_for_hex(a), partner=b)
            yield run_sol(add_parts_for_hex(b), partner=a)
    else:
        assert False

    it = (out_item for soln in load_input() for out_item in handle_item(soln))
    dump = lambda x: json.dump(x, sys.stdout)

    UNIQUE_KEYS = set()
    def check_key_uniqueness(x):
        key = x['key']['string']
        if key in UNIQUE_KEYS:
            raise AssertionError('Duplicate key: {!r}'.format(key))
        UNIQUE_KEYS.add(key)
        return x

    def passthru(testFunc, x):
        testFunc(x)
        return x

    if PARANOID >= 1:
        it = (check_key_uniqueness(passthru(do_special_case_sanity_checks, x)) for x in it)

    if args.stream_out:
        for x in it:
            if x: dump(x)
    else: dump([x for x in it if x])


def main_(soln, key_layout, min_volume, max_volume, partner=None):
    beta, (a, b, c), rparts = soln
    rn = int(rparts['numerator'])
    rd = int(rparts['denominator'])
    rk = int(rparts['rootNum'])
    r = S(rn) / rd * sqrt(rk) # scale factor (for layers of different cell sizes, e.g. graphene and h-BN)

    # family matrices and sites ought to be added to input
    #  before we start trying anything but hex
    assert beta == 3

    # HACK: show the ABC_SCALE key for progress feedback. (VOLUME_LETTER key cannot be computed yet)
    print("handling {}-{}-{}-{}-{}-{}".format(a,b,c,rn,rd,rk), file=sys.stderr)

    sites = [[0, 0], [S(2) / 3, S(1) / 3]]
    A = Matrix([[1, 0], [-S(1) / 2, sqrt(3) / 2]])

    M = S(r) / c * Matrix([[a, -b * sqrt(beta)], [b * sqrt(beta), a]])
    AM = mdot(A, M.T)
    AMM = mdot(A, M.T, M.T)

    kw = {
        'max_volume': max_volume,
        'min_volume': min_volume,
    }
    positions = do_multi_layer(sites, [A, AM], **kw)

    def compute_key():
        if key_layout == Layout.ABC_SCALE:
            key_parts = [a, b, c, rn, rd, rk]
            key_string = '-'.join(map(str, map(int, key_parts)))

        elif key_layout == Layout.VOLUME_DISAMBIG:
            assert rn == rd == rk == 1, "nontrivial scale not supported"
            assert partner is not None, "need partner; use rust-pairs format"

            # This encoding fixes the non-uniqueness of VOLUME_LETTER by incorporating the
            #  value of 'c - a' for the solution with angle < 30.

            # (why 'c - a'? because it's usually a tiny number)

            import math
            if (a,b,c) == (1,1,2): # 60 degrees exact
                letter = 'b'
            else: # we can trust floating point precision for the rest
                letter = chr(ord('a') + int(math.acos(a/c) // (math.pi / 6)))

            if letter == 'a':
                effective_a = a
                effective_c = c
            else:
                assert letter == 'b'
                effective_a, _, effective_c = partner

            v = int(positions['meta']['volume'][0])

            # Now to demonstrate that, unlike VOLUME_LETTER, this encoding is unique.
            #
            # Well... not quite. The following is taken as conjecture:
            #
            #    Solution triples (a, b, c) for r == 1 whose rotation angles are in [0.0, 30.0 deg]
            #    are uniquely labeled by (c-a, v), where v is the solution volume.
            #
            # I am as yet unable to prove this... BUT! I verified it by brute force for
            # all solutions with c < 10_000_000, which gives angles down to 0.018121 degrees
            assert v < 10_000_000, "conjecture not proven for such large volume"

            disambig = effective_c - effective_a
            key_parts = [v, disambig, letter]
            key_string = f'{v:03d}-{disambig}-{letter}'
        else:
            raise RuntimeError("incomplete switch for Layout")

        return dict(layout=key_layout, parts=key_parts, string=key_string)

    key = compute_key()

    return {
        'key': key,
        'solution': {
            'abc': [a, b, c],
            'β': beta,
            'families': [
                [[1, 2], [1, 2]],
                [[1, 2], [1, 2]],
            ],
            'r': {
                'square': [rn * rn * rk, rd * rd],
                'exact': [rn, rd, rk],
                'approx': float(r),
            },
        },
        'positions': positions,
    }

def checked_div(a, b):
    assert a % b == 0
    return a // b

def is_square(x):
    f = float(x)
    assert x == f, "too big for float mantissa"
    return round(f**0.5)**2 == x

assert list(map(is_square, range(10))) == list(map(bool, [1, 1, 0, 0, 1, 0, 0, 0, 0, 1]))

# NOTE: This used to be used to match e.g. 3 or more layers but
#       currently it's only ever used on two
def do_multi_layer(basicSites, units, *, max_volume, min_volume):
    from functools import reduce
    # make a cell commensurate of all layers
    SC = reduce(
        (lambda A,B: MoirePattern.from_cells(A,B).commensurate_cell()),
        units,
    )
    SC = find_nicer_cell(SC)
    # rotate/reflect basis for lammps
    (SC, trans) = lammps_friendly_cell(SC)
    SC = no_60s_allowed(SC)

    units = [mdot(A, trans) for A in units]
    coeffs = [mdot(SC, A.inv()) for A in units]
    volumes = [abs(C.det()) for C in coeffs]
    if max_volume is not None and S(max_volume) <= S(volumes[0]):
        return None
    if S(volumes[0]) < S(min_volume):
        return None

    # FIXME:  HACK:
    # check that the cell is "standard" for hexagonal;
    # both vectors should be of equal length
    assert sqnorm(SC[:2]) == sqnorm(SC[2:])

    maxIndices = [supercellMaxIndices(C) for C in coeffs]

    if PARANOID >= 1:
        validate_standard_hex_cell(SC)
        for A in units:
        # for (A, sites, latts) in zip(units, allSites, allLatts):
            # validate_hexagonal_shape(SC, latts)
            # validate_honeycomb_shape(SC, sites)
            validate_standard_hex_cell(A)

    uniter = lambda f, it: [f(x) for x in it]
    uniter2 = lambda f, it: [[f(x) for x in xs] for xs in it]
    unmat = lambda f, M: uniter2(f, M.tolist())

    return {
        'lattice': unmat(float, SC),
        'layer': [
            {
                'frac-sites': uniter2(float, basicSites),
                'frac-lattice': unmat(float, C.inv()),
                'repeat': uniter(int, idx),
            } for (C, idx) in zip(coeffs, maxIndices)
        ],
        'meta': {
            'layer': [
                {
                    'cart': {'approx': unmat(float, A)},
                    'frac': {'approx': unmat(float, C.inv())},
                } for (A, C) in zip(units, coeffs)
            ],
            'coeff': [ unmat(int, C) for C in coeffs ],
            'volume': [ abs(int(C.det())) for C in coeffs ],
        },
    }


def cut_out_third_layer(d):
    # We fully specify how the dictionary should be transformed
    #  to be sure we don't miss anything.
    ACTION_KEEP = object()
    ACTION_TAKE_2_OF_3 = object()
    SPEC = {
        'lattice': ACTION_KEEP,
        'layer': ACTION_TAKE_2_OF_3,
        'meta': {
            'layer': ACTION_TAKE_2_OF_3,
            'coeff': ACTION_TAKE_2_OF_3,
            'volume': ACTION_TAKE_2_OF_3,
        },
    }

    def transform_by_spec(spec, d):
        if spec is ACTION_KEEP:
            return d
        elif spec is ACTION_TAKE_2_OF_3:
            assert isinstance(d, list)
            assert len(d) == 3
            return d[:2]
        elif isinstance(spec, dict):
            return zip_dict_with(transform_by_spec, spec, d)
        else:
            assert False, "complete switch"

    return zip_dict_with(transform_by_spec, SPEC, d)

def zip_dict_with(func, d1, d2):
    assert isinstance(d1, dict)
    assert isinstance(d2, dict)
    if set(d1) != set(d2):
        raise ValueError("dict keysets not parallel")
    return { k:func(d1[k], d2[k]) for k in d1 }


def lammps_friendly_cell(m):
    m = Matrix(m)
    # This cannot be achieved solely through a unitary transform;
    # An additional cartesian transformation must be performed, which is
    #  returned in 'trans' (and must be applied to A and B to preserve
    #  the supercell matrices C and D)
    # rotate a to +x
    rot = rot_axis3(-atan2(m[1], m[0]))

    # sympy rotation matrix is apparently [[+cos,+sin],[-sin,cos]].
    # I can't say I've *ever* seen this convention before... :/
    rot = rot.T

    rot.col_del(2)
    rot.row_del(2)

    trans = simplify(rot.T)
    m = mdot(m, trans).as_mutable()

    # maybe negate second basis vector
    # (this is unitary and can be left out of 'trans')
    if m[3] < 0:
        m[2] = -m[2]
        m[3] = -m[3]
    # check
    assert m[1] == 0, str(m)
    assert m[0] > 0, str(m)
    assert m[3] > 0, str(m)
    return (m, trans)

HNF_SEARCH_BUILT = False
def run_hnf_search():
    import subprocess
    from subprocess import Popen, PIPE

    global HNF_SEARCH_BUILT
    if not HNF_SEARCH_BUILT:
        subprocess.run(['cargo', 'build', '--release'], cwd='rust/hnf-search', check=True)
        HNF_SEARCH_BUILT = True

    return Popen('rust/hnf-search/target/release/hnf-search', stdin=PIPE, stdout=PIPE, stderr=PIPE)

# Input:  Integer supercell matrix
# Output: (iMax, jMax), the max number of unique images along each axis
#         (equivalently, the diagonal of the HNF of C)
def supercellMaxIndices(m):
    # HACK
    # I never thought this little rust thing would see the light of day (especially since it's
    # "yet another implementation" of an algorithm that shows up 50 bajillion times in this code
    # base), but when it comes to unit cells that are tens of thousands of atoms,
    # python is *just too slow*.
    import subprocess

    p = run_hnf_search()
    mInv = m.inv().tolist()
    mInv = [[x.as_numer_denom() for x in row] for row in mInv]
    mInv = [["{}/{}".format(n,d) for (n,d) in row] for row in mInv]
    mInv = '[{}]'.format(','.join('[{}]'.format(','.join(row)) for row in mInv))
    (out, err) = p.communicate(mInv.encode('utf-8'))
    import json
    [[c00,_],[_,c11]] = json.loads(out.decode('utf-8'))
    return c00,c11

    # SLOW PYTHON VERSION

    # getfrac = fracModMatrix(m)
    # seen = set()

    # # how many points to a row? (all rows will share the same number, because the lengths
    # #  of parallel lines between two parallel lines are equal)
    # for j in itertools.count(1):
    #     if getfrac((0, j)) == (0, 0):
    #         nj = j
    #         break

    # for i in itertools.count(0):
    #     # find first new column in this row
    #     for j in range(nj):
    #         if getfrac((i, j)) not in seen:
    #             break
    #     # went a whole period; there's nothing new left
    #     else:
    #         ni = i
    #         break

    #     ijs = [(i, j) for j in range(j, j + nj)]
    #     fracs = set(getfrac(ij) for ij in ijs)
    #     assert len(ijs) == len(fracs)
    #     seen.update(fracs)
    # return (ni, nj)

# =================================
# NOTE: These all take a single vector as their second (curried) argument.

# IMPORTANT: These are for row-based vector formalisms
def mulMatrix(m):
    return lambda ij: tuple((Matrix([list(ij)]) * m))


def truedivMatrix(m):
    mInv = m.inv()
    return lambda ij: tuple((Matrix([list(ij)]) * mInv))


def fracModMatrix(m):
    mInv = m.inv()
    return lambda ij: tuple((Matrix([list(ij)]) * mInv).applyfunc(lambda x: x % 1))


def modMatrix(m):
    mInv = m.inv()
    return lambda ij: tuple((Matrix([list(ij)]) * mInv).applyfunc(lambda x: x % 1) * m)


def fracModMatrixMany(m, ij):
    # One might think that, no matter how inefficient the implementation of matrix
    #  multiplication is in sympy, it could not possibly be so inefficient that there
    #  would be any tangible speedup by switching to numpy dtype=object arrays.
    #
    # ...One would, of course, be dead wrong.
    import numpy as np
    mInv = m.inv()
    mInv = np.array(mInv.tolist())
    ij = np.array(list(map(list, ij)))
    ij = mdot(ij, mInv) % 1
    return list(map(tuple, ij.tolist()))
# =================================


def diff(xy1, xy2):
    assert len(xy1) == len(xy2)
    return tuple(a - b for (a, b) in zip(xy1, xy2))


def norm(x):
    return sqrt(sqnorm(x))


def sqnorm(x):
    return dot(x, x)


def dot(x, y):
    (x, y) = (list(x), list(y))
    assert len(x) == len(y)
    return sum(a * b for (a, b) in zip(x, y))


def cross2(xy1, xy2):
    assert len(xy1) == 2 and len(xy2) == 2
    return xy1[0] * xy2[1] - xy1[1] * xy2[0]

def mdot(*args):
    """ A single consistent syntax for matrix multiplication.
    Of course, such was the point of the matrix multiplication operator,
    but for some reason numpy decided to stop supporting that for arrays
    of objects. """
    import numpy as np
    import sympy
    import operator
    from functools import reduce
    is_numpy = lambda x: isinstance(x, np.ndarray)
    is_sympy = lambda x: isinstance(x, (sympy.Matrix, sympy.ImmutableMatrix))

    if is_numpy(args[0]):
        if not all(map(is_numpy, args)):
            raise TypeError(set(map(type, args)))
        return reduce(np.dot, args)

    elif is_sympy(args[0]):
        if not all(map(is_sympy, args)):
            raise TypeError(set(map(type, args)))
        return reduce(operator.mul, args)

    else:
        raise TypeError(type(args[0]))


# Ensure that the cell is the one used on the bilbao server,
# which apparently has an interior angle of 120 degrees, not 60.
# (this distinction is important, because it affects the location of K)
def validate_bilbao_hex_cell(m):
    top, bot = m[:2], m[2:]
    assert sqnorm(top) == sqnorm(bot)
    assert dot(top, bot) / sqnorm(
        bot) == -S(1) / 2  # Note: +1/2 would imply the angle is 60


# precondition:   M has hexagonal symmetry with cell angle == 60 or 120
# postcondition:  M has hexagonal symmetry with cell angle == 120
def no_60s_allowed(M):
    M = Matrix(M)
    if dot(M[2:], M[:2]) / dot(M[2:], M[2:]) == +S(1) / 2:
        M[2] -= M[0]
        M[3] -= M[1]
    assert dot(M[2:], M[:2]) / dot(M[2:], M[2:]) == -S(1) / 2
    return M

def iterrows(M):
    for i in M.rows:
        yield list(M[i,:])

def validate_standard_hex_cell(M):
    # vectors of equal length
    assert sqnorm(M[0,:]) == sqnorm(M[1,:]), M
    # obtuse angle
    assert dot(M[0,:], M[1,:]) < 0, M
    # 120 degrees; results in the property that (a + b) is similar to a and b
    assert sqnorm(M[0,:] + M[1,:]) == sqnorm(M[0,:]), M

def do_special_case_sanity_checks(d):
    if d['key']['string'] == '1-0-1-1-1-1':
        p = d['positions']
        for C in p['meta']['coeff']:
            assert C == [[1,0],[0,1]]
        A,*Bs = p['meta']['layer']
        for B in Bs:
            assert A == B

if __name__ == '__main__':
    main()
