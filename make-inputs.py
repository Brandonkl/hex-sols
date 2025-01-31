#!/usr/bin/env python3

# The role of this script is to take the single file produced by positions.py
# and generate MANY files, in a hierarchial structure to be used by Shakefile-inner.
#
# Primary responsibilities include:
#
#  * what parameters are our computations parameterized over?
#    (i.e. how many levels of directory nesting, and what do they mean?)
#  * which structures do we care about?
#
# ...though I'd be lying if I claimed you could "just edit this file."
# Have fun fixing the Shakefiles!!
#
# As a general rule of thumb, if it needs numpy, it shouldn't be done in this file.

import os
import sys
import json
import yaml # easier to work with in haskell...
import itertools

def main():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('-S', '--params', default='general-spatial-params.yaml')
    parser.add_argument('-P', '--positions', default='positions.json')
    parser.add_argument('-I', '--ignore-from', default=[], action='append', help='file with lines listing keys to ignore')
    parser.add_argument('-W', '--whitelist-from', default=[], action='append', help='whitelist file (do only these keys)')
    parser.add_argument('-o', '--output-dir', default='data', help='directory to populate with output')
    # HACK
    parser.add_argument('--layer-index', default=None,
        help='python3 list expression indicating layer indices to draw from.'
        ' Default is to use all layers as supplied, e.g. "[0,1,2,3]" for an input with 4 layers.')
    parser.add_argument('--layer-shift', default=None,
        help='python3 expression for a list of length-2 lists;'
        ' one for each layer in --layer-index.')
    args = parser.parse_args()
    
    spatial_params_base = yaml.load(open(args.params), Loader=yaml.FullLoader)
    # target-supercell: Lower bounds on supercell sidelength, in units of the single-layer graphene cell
    DEFAULT_TARGET_SUPERCELL = {
        'phonopy': 0,  # used to set phonopy.supercell_dim in sp2's config.json
        'wobble': 0,   # used to set display size for animations
    } # (a lower bound of 0 will produce 1x1x1)
    target_supercell = spatial_params_base.pop('target-supercell')
    if set(target_supercell) != set(DEFAULT_TARGET_SUPERCELL):
        parser.error('incorrect set of keys in target-supercell ({} != {})'.format(set(target_supercell), set(DEFAULT_TARGET_SUPERCELL)))

    arg_layer_index = None if args.layer_index is None else eval(args.layer_index)
    arg_layer_shift = None if args.layer_shift is None else eval(args.layer_shift)

    ignore = set()
    whitelist = set()
    for path in args.ignore_from:
        with open(path) as f:
            ignore.update([x.strip() for x in f if x.strip()])
    for path in args.whitelist_from:
        with open(path) as f:
            whitelist.update([x.strip() for x in f if x.strip()])

    ensure_dir(args.output_dir)

    for d in json.load(open(args.positions)):
        soln_id = d['key']['string']
        if soln_id in ignore: continue
        if whitelist and soln_id not in whitelist: continue

        soln_dir = ensure_subdir(args.output_dir, soln_id)

        # HACK: this used to be something like {'ab': { ... }, 'abc': { ... }}
        #       but now it is just { ... }
        for positions_d in [d['positions']]:
            # layer_style_dir = ensure_subdir(soln_dir, layer_style)

            # reduce number of places where we need to rename a variable...
            final_dir = soln_dir

            # Input to the assemble script, with reasonably configurable params
            with open(os.path.join(final_dir, 'spatial-params.yaml'), 'w') as f:
                yaml.dump(spatial_params_base, f)

            # Another input to the assemble script, describing the structure in essence
            with open(os.path.join(final_dir, 'layers.yaml'), 'w') as f:
                layers = positions_d['layer']

                if arg_layer_index is not None:
                    layers = [layers[i] for i in arg_layer_index]
                if arg_layer_shift is not None:
                    assert len(arg_layer_shift) == len(layers)
                    assert all(len(x) == 2 for x in arg_layer_shift)
                    layers = [dict(shift=x, **d) for (d, x) in zip(layers, arg_layer_shift)]
                yaml.dump({
                    'lattice': positions_d['lattice'],
                    'layer': layers,
                }, f)

            # Make a flatter positions.json that collapses ['positions']['abc'] into the main tree.
            # As of me writing this comment, positions.json is significantly less important
            #  than it used to be, and mostly exists for direct user reference now.
            # (though the shakefile might peek inside at the coefficient matrices)
            with open(os.path.join(final_dir, 'positions.json'), 'w') as f:
                new_d = dict(d)
                new_d.pop('positions')
                new_d = dict(new_d, **positions_d)
                new_d.pop('layer') # this already has its own file
                json.dump(new_d, f)

            with open(os.path.join(final_dir, 'supercells.json'), 'w') as f:
                def best_supercell(target):
                    # (work with volume for sake of exactness)
                    volume = positions_d['meta']['volume'][0]
                    for f in itertools.count(1):
                        if f*f*volume >= target*target:
                            return [f,f,1]
                supercells = map_dict(best_supercell, target_supercell)

                json.dump(supercells, f)

def map_dict(func, d):
    assert isinstance(d, dict)
    return {k:func(v) for (k,v) in d.items()}

def ensure_dir(path):
    try: os.mkdir(path)
    except FileExistsError: pass

def ensure_subdir(dirname, basename):
    path = os.path.join(dirname, basename)
    ensure_dir(path)
    return path

main()
