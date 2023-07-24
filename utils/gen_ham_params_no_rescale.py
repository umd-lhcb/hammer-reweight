#!/usr/bin/env python
#
# Author: Yipeng Sun, Alex Fernez

# Compared to gen_ham_params.py, this script generates variations with rescaling removed
# Motivation: rescaling the FFs has the effect of simply multiplying the diff cross section
# by an overall factor, which won't affect template shapes, and in fact allowing this 
# effective variation in the fit can lead to stability issues
# Math: assuming FF parameters f have (potentially non-diagonal) covariance U, the first step is to
# diagonalize to get the variations in an uncorrelated basis: U=VLV^T (in the uncorrelated basis
# 1sigma variations of f' look like f' -> f' +- sigma'_i*e_i [e_i this ith unit vector, sigma'_i is
# the sqrt of the ith eigenvalue found when diagonalizing], corresponding to 
# f=Af' -> f +- sigma'_i*v_i [v_i is the ith column of V; the ith eigenvector found from 
# diagonalizing]). Next, the overlap of the variations with the rescaling direction (the direction
# of the nominal FF params, n) can be removed, creating w_i=v_i-(v_i.n)n, which are now 
# unnormalized, non-orthogonal, and (potentially) don't correspond to variations in an
# uncorrelated basis. To account for this, define normalized vectors x_i s.t. sigma'_i*w_i = s_i*x_i
# (s_i is the length of the ith variation with rescaling removed) and an effective covariance matrix
# U'=XS^2X^T (where the columns of X are the x_i and S is a diagonal matrix with the s_i on the 
# diag; note that X^T != X^-1 here!). Finally, again diagonalize: U'=AZA^T (now A^T=A^-1), where Z
# will have one 0 on the diagonal corresponding to an eigenvector (column of A) equal to n. The 
# columns of A, excluding the one equal to n, multiplied by the sqrt of the corresponding eigenvalue
# in Z, are the final variations with rescaling removed

import yaml
import numpy as np

from argparse import ArgumentParser


#######################
# Command line parser #
#######################

def parse_input():
    parser = ArgumentParser(description="generate HAMMER parameters and variations, \
                                         excluding rescaling.")

    parser.add_argument("input", help="specify input YAML.")

    parser.add_argument(
        "-m",
        "--modes",
        nargs="+",
        default=[
            "BtoDBGLVar",
            "BtoDstarBGLVar",
            "BtoD0starBLRVar",
            "BtoD1starBLRVar",
            "BtoD1BLRVar",
            "BtoD2starBLRVar",
        ],
        help="specify FF params to print",
    )

    return parser.parse_args()

###########
# Helpers #
###########

# need to determine if some numpy scalars/vectors are equal, modulo some floating point arithmetic
# and overall minus signs
def almost_equal(a, b, z=1e-14):
    return np.linalg.norm(a-b) < z or np.linalg.norm(a+b) < z

def print_param_general(ff_alias, param, val, add_hammer=False):
    if add_hammer:
        print(f'  ham.setOptions(scheme + ": {{{param}: {val}}}");')
    else:
        print(f'    "{{{param}: {val}}}",')

def print_param_ff_var(process, model, shifts, params, comments):
    print("  FF variation params:")
    for single, c in zip(shifts, comments):
        eigen_vector_spec = [f'{{"{p}", {s}}}' for p, s in zip(params, single)]
        eigen_vector_str = "{" + ", ".join(eigen_vector_spec) + "}"

        print(
            f'    ham.setFFEigenvectors{{"{process}", "{model}", {eigen_vector_str}}}; // {c};'
        )

def fmt_dict_as_cpp_map(dct):
    output = "{"
    for key, value in dct.items():
        elem = f'{{"{key}", {value}}}, '
        output += elem
    return output[:-2] + "},"

def eval_fake_sandbox(code, add_vars):
    loc = dict()
    sandbox = {k: v for k, v in globals().items()}
    sandbox.update(add_vars)
    exec(f"__result = {code}", sandbox, loc)
    return loc["__result"]

def remove_rescale(Vars, nom):
    # take in n variations that include rescaling component, return n-1 variations with rescaling
    # removed (for notation, see Math comments at beginning of script)
    nom = np.array(nom)/np.linalg.norm(nom) # just need the direction
    Vars_minus_rescale = Vars - (np.array([np.dot(col,nom)*nom for col in Vars.T])).T # sigma'_i*w_i
    X = (np.array([col/np.linalg.norm(col) for col in Vars_minus_rescale.T])).T
    S = np.diag([np.linalg.norm(col) for col in Vars_minus_rescale.T])
    Up = np.matmul(X, np.matmul(S**2, X.T))
    z, A = np.linalg.eig(Up) # z is a vector (of variances) here
    # do a quick check that there's exactly 1 zero eval with corresponding evec = nom
    nom_entry = -1
    for i in range(len(z)):
        val = z[i]
        if almost_equal(val, 0):
            if not nom_entry == -1:
                print(f'\n\nMore than one 0 eigenval found when removing rescale!!! ({val})\n\n')
            nom_entry = i
            if not almost_equal(nom, A.T[nom_entry], 1e-8): # a little more forgiving for vectors
                print(f'\n\nFound 0 eigenval with eigenvec not equal to nom direc!!!\n{nom}\
                      \n{A.T[nom_entry]}\n\n')
    # assuming nothing gets printed out, can safely return variations with rescaling removed!
    z = np.concatenate((z[:nom_entry],z[nom_entry+1:]))
    A = (np.concatenate((A.T[:nom_entry],A.T[nom_entry+1:]))).T
    return z**(1/2)*A

#################
# Model helpers #
#################

def gen_param_shifted_BtoDBGL(
    process, model, m_corr, v_err, add_params, ap, a0, verbose=True
):
    # build the covariance matrix
    C = np.array(m_corr)
    Std = np.diag(np.array(v_err))
    U = np.matmul(Std, np.matmul(C, Std))
    # extract uncorrelated variations
    L, V = np.linalg.eig(U) # evals (variances), evecs
    Vars = L**(1/2)*V
    # remove rescaling
    Vars = remove_rescale(Vars, ap[:-1]+a0[1:-1]) # ap/a0 have 0 appended as last element, a00 fixed

    print()
    print(f"{process}{model} 1sigma variated (with rescaling removed) FF values:")
    for i in range(5-1):
        var_values = Vars.T[i]
        # compute shift for ap's
        var_ap = var_values[:3]
        var_ap = np.append(var_ap, 0.0)
        # compute shift for a0'
        var_a0 = var_values[3:]
        var_a00 = np.dot(np.array(add_params["a00"]), var_values)
        var_a0 = np.insert(var_a0, 0, var_a00)
        var_a0 = np.append(var_a0, 0.0)

        if verbose:
            print(f"  // shifting in {i+1}-th direction (+)...")
            print("  {")
        print_param_general(process + model, "ap", list(np.array(ap) + var_ap))
        print_param_general(process + model, "a0", list(np.array(a0) + var_a0))
        if verbose:
            print("  },")

        if verbose:
            print(f"  // shifting in {i+1}-th direction (-)...")
            print("  {")
        print_param_general(process + model, "ap", list(np.array(ap) - var_ap))
        print_param_general(process + model, "a0", list(np.array(a0) - var_a0))
        if verbose:
            print("  },")

def gen_param_shifted_BtoDstarBGL(
    process,
    model,
    m_corr,
    avec,
    bvec,
    cvec,
    dvec,
    aerr,
    berr,
    cerr,
    derr,
    scale,
    verbose=True,
):
    # build the covariance matrix
    C = np.array(m_corr)
    Std = np.diag(np.array(aerr + berr + cerr + derr) * scale) # scale error w/ Vcb * eta_ew as well!
    U = np.matmul(Std, np.matmul(C, Std))
    # extract uncorrelated variations
    L, V = np.linalg.eig(U) # evals (variances), evecs
    Vars = L**(1/2)*V
    # remove rescaling
    Vars = remove_rescale(Vars, avec+bvec+cvec+dvec)

    print()
    print(f"{process}{model} 1sigma variated (with rescaling removed) FF values:")
    for i in range(12-1):
        var_values = Vars.T[i]
        var_a = var_values[:3]
        var_b = var_values[3:6]
        var_c = var_values[6:9]
        var_d = var_values[9:12]

        if verbose:
            print(f"  // shifting in {i+1}-th direction (+)...")
            print("  {")
        for name, nom, delta in zip(
            ["avec", "bvec", "cvec", "dvec"],
            [avec, bvec, cvec, dvec],
            [var_a, var_b, var_c, var_d],
        ):
            print_param_general(process + model, name, list(np.array(nom) + delta))
        if verbose:
            print("  },")

        if verbose:
            print(f"  // shifting in {i+1}-th direction (-)...")
            print("  {")
        for name, nom, delta in zip(
            ["avec", "bvec", "cvec", "dvec"],
            [avec, bvec, cvec, dvec],
            [var_a, var_b, var_c, var_d],
        ):
            print_param_general(process + model, name, list(np.array(nom) - delta))
        if verbose:
            print("  },")

def gen_param_shifted_BtoDstarstarBLR(
    process, model, m_corr, v_nom, v_err, param_names, verbose=True
):
    # if don't fully trust paper constraints, optionally increase 1sigma variations by a factor
    delta_factor = 2

    # build the covariance matrix
    C = np.array(m_corr)
    Std = np.diag(np.array(v_err))
    U = np.matmul(Std, np.matmul(C, Std))
    # extract uncorrelated variations
    L, V = np.linalg.eig(U) # evals (variances), evecs
    Vars = L**(1/2)*V
    # remove rescaling
    Vars = remove_rescale(Vars, v_nom)

    print()
    print(f"{process}{model} 1sigma variated (with rescaling removed) FF values:")
    for i in range(len(v_nom)-1):
        var_values = Vars.T[i]

        if verbose:
            print(f"  // shifting in {i+1}-th direction (+)...")
            print("  {")
        for name, nom, delta in zip(param_names, v_nom, var_values):
            print_param_general(process + model, name, nom + delta_factor*delta)
        if verbose:
            print("  },")

        if verbose:
            print(f"  // shifting in {i+1}-th direction (-)...")
            print("  {")
        for name, nom, delta in zip(param_names, v_nom, var_values):
            print_param_general(process + model, name, nom - delta_factor*delta)
        if verbose:
            print("  },")

########
# Main #
########

if __name__ == "__main__":
    args = parse_input()

    with open(args.input) as f:
        cfg = yaml.safe_load(f)

    for ff_model in [m for m in args.modes if m in cfg]:
        print(f"Handling {ff_model}...")
        for param, val in cfg[ff_model].items():
            if param.startswith("_"):
                continue  # temp variables, skip them

            if not isinstance(val, str):
                print_param_general(ff_model, param, val, True)
            else:
                result = eval_fake_sandbox(val, cfg[ff_model])
                if result:
                    print_param_general(ff_model, param, result, True)
                    cfg[ff_model][param] = result
        print()
