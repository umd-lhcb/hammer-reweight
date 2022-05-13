#!/usr/bin/env python
#
# Author: Yipeng Sun
# Last Change: Fri May 13, 2022 at 04:00 AM -0400

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

#################
# Model helpers #
#################

def gen_cov_mat(m_cov, v_err):
    pass


###########
# Helpers #
###########

def print_param_shift(process, model, shifts, params, comments):
    print('  Param shifts:')
    for single, c in zip(shifts, comments):
        eigen_vector_spec = [f'{{"{p}", {s}}}' for p, s in zip(params, single)]
        eigen_vector_str = '{' + ', '.join(eigen_vector_spec) + '}'

        print(f'    ham.setFFEigenvectors{{"{process}", "{model}", {eigen_vector_str}}} // {c};')


