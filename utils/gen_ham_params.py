#!/usr/bin/env python
#
# Author: Yipeng Sun
# Last Change: Mon May 30, 2022 at 05:00 PM -0400

import yaml
import numpy as np

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

###########
# Helpers #
###########

def print_param_general(ff_alias, param, val):
    print(f'  ham.setOptions("{ff_alias}, {{{param}, {val}}}")')


def print_param_ff_var(process, model, shifts, params, comments):
    print('  FF variation params:')
    for single, c in zip(shifts, comments):
        eigen_vector_spec = [f'{{"{p}", {s}}}' for p, s in zip(params, single)]
        eigen_vector_str = '{' + ', '.join(eigen_vector_spec) + '}'

        print(f'    ham.setFFEigenvectors{{"{process}", "{model}", {eigen_vector_str}}}; // {c};')


def fmt_dict_as_cpp_map(dct):
    output = '{'
    for key, value in dct.items():
        elem = f'{{"{key}", {value}}}, '
        output += elem
    return output[:-2] + '},'


def eval_fake_sandbox(code, add_vars):
    loc = dict()
    sandbox = {k: v for k, v in globals().items()}
    sandbox.update(add_vars)
    exec(f'__result = {code}', sandbox, loc)
    return loc['__result']


#################
# Model helpers #
#################

def gen_param_var(process, model, m_corr, v_err, param_names, add_params):
    m_corr = np.matrix(m_corr)
    v_err = np.array(v_err)

    # make the covariance matrix based on the correlation and errors
    m_cov = np.einsum('ij,i,j->ij', m_corr, v_err, v_err)
    v_eigen, m_eigen = np.linalg.eig(m_cov)
    # m_eigen is the matrix formed by eigen vectors; v_eigen is the vector
    #   formed by eigenvalues
    # m_eigen is essentially M in the ANA
    # now we need to find C in the ANA to construct A
    m_c = np.einsum('ij,i->ij', np.eye(v_eigen.size), (1 / np.sqrt(v_eigen)))
    m_a = np.einsum('ik,kj', m_c, m_eigen)
    # NOTE: uncomment the line below to reproduce RD+'s numbers
    #m_a = m_eigen
    # now find the inverse of A, A^-1 tells us how to transform from
    #  an ERROR eigenbasis to our NORMAL FF paramterization basis (HAMMER basis)
    m_a_inv = np.linalg.inv(m_a)

    # define variations in the ERROR eigenbasis:
    # m_variation_ham: each column represents a variation of +1sigma for a
    # ERROR eigen vector
    # the output will be: u1p, u1m, u2p, u2m, ...
    print()
    print(f'{process}{model} FF variations:')
    for i in range(len(param_names)):
        var_values = m_a_inv[:,i]
        var_pos = {'delta_'+k: v for k, v in zip(param_names, var_values)}
        for name, coeff in add_params.items():
            var_pos['delta_'+name] = np.dot(np.array(coeff), var_values)

        var_neg = {k: -v for k, v in var_pos.items()}
        print(f'  {fmt_dict_as_cpp_map(var_pos)}')
        print(f'  {fmt_dict_as_cpp_map(var_neg)}')


def gen_param_shifted_BtoDBGL(
        process, model, m_corr, v_err, add_params, ap, a0):
    m_corr = np.matrix(m_corr)
    v_err = np.array(v_err)

    m_cov = np.einsum('ij,i,j->ij', m_corr, v_err, v_err)
    v_eigen, m_eigen = np.linalg.eig(m_cov)
    m_c = np.einsum('ij,i->ij', np.eye(v_eigen.size), (1 / np.sqrt(v_eigen)))
    m_a = np.einsum('ik,kj', m_c, m_eigen)
    m_a_inv = np.linalg.inv(m_a)

    print()
    print(f'{process}{model} shifted central values:')
    for i in range(5):
        var_values = m_a_inv[:,i]
        # compute shift for ap's
        var_ap = var_values[:3]
        var_ap = np.append(var_ap, 0.0)
        # compute shift for a0'
        var_a0 = var_values[3:]
        var_a00 = np.dot(np.array(add_params['a00']), var_values)
        var_a0 = np.insert(var_a0, 0, var_a00)
        var_a0 = np.append(var_a0, 0.0)

        print(f' shifting {i+1}-th param in + direction....')
        print_param_general(process+model, 'ap', list(np.array(ap) + var_ap))
        print_param_general(process+model, 'a0', list(np.array(a0) + var_a0))

        print(f' shifting {i+1}-th param in - direction....')
        print_param_general(process+model, 'ap', list(np.array(ap) - var_ap))
        print_param_general(process+model, 'a0', list(np.array(a0) - var_a0))



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
                eval_fake_sandbox(val, cfg[ff_model])
