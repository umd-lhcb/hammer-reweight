#!/usr/bin/env python
#
# Author: Yipeng Sun
# Last Change: Fri May 13, 2022 at 01:43 PM -0400

import yaml

from argparse import ArgumentParser


#######################
# Command line parser #
#######################

def parse_input():
    parser = ArgumentParser(
        description='generate HAMMER parameters.')

    parser.add_argument('input', help='specify input YAML.')

    parser.add_argument('-m', '--modes', nargs='+',
                        default=['BtoDBGLVar', 'BtoDstarBGLVar'],
                        help='specify FF params to print')

    return parser.parse_args()

#################
# Model helpers #
#################

def gen_cov_mat_BtoDBGL(m_cov, v_err):
    pass


###########
# Helpers #
###########

def print_param_general(ff_alias, param, val):
    print(f'  ham.setOptions("{ff_alias}", {{{param}, {val}}})')


def print_param_ff_var(process, model, shifts, params, comments):
    print('  FF variation params:')
    for single, c in zip(shifts, comments):
        eigen_vector_spec = [f'{{"{p}", {s}}}' for p, s in zip(params, single)]
        eigen_vector_str = '{' + ', '.join(eigen_vector_spec) + '}'

        print(f'    ham.setFFEigenvectors{{"{process}", "{model}", {eigen_vector_str}}}; // {c};')


########
# Main #
########

if __name__ == '__main__':
    args = parse_input()

    with open(args.input) as f:
        cfg = yaml.safe_load(f)

    for ff_model in [m for m in args.modes if m in cfg]:
        print(f'Handling {ff_model}...')
        for param, val in cfg[ff_model].items():
            if param.startswith('_'):
                continue  # temp variables, skip them

            if not isinstance(val, str):
                print_param_general(ff_model, param, val)
            else:
                pass
