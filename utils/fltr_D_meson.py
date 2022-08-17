#!/usr/bin/env python
#
# Author: Yipeng Sun
# Last Change: Thu Sep 30, 2021 at 03:14 PM +0200

import uproot
import numpy as np

from argparse import ArgumentParser
from pyTuplingUtils.io import read_branches, read_branches_dict


def parse_input():
    parser = ArgumentParser(description="filter particles and print its truth info.")

    parser.add_argument("ntp", help="input ntuple.")

    parser.add_argument(
        "-t", "--tree", default="TupleB0/DecayTree", help="specify tree name."
    )

    parser.add_argument(
        "-i",
        "--id",
        nargs="+",
        type=int,
        default=[10421, 10411],
        help="specify particle true ID.",
    )

    parser.add_argument(
        "-b",
        "--base-br",
        default="d_meson1_true_id",
        help="specify the true ID branch.",
    )

    parser.add_argument(
        "-w", "--weight-br", default="wff", help="specify the weight branch."
    )

    parser.add_argument(
        "-o", "--ok-br", default="ham_ok", help="specify the cut branch."
    )

    parser.add_argument(
        "-B",
        "--add-br",
        nargs="+",
        default=[
            "d_meson1_true_m",
        ],
        help="specify additional branches.",
    )

    parser.add_argument(
        "-N",
        "--name",
        nargs="+",
        default=["inv.mass."],
        help="specify printout names of the additional branches.",
    )

    return parser.parse_args()


def fltr_func(true_ids):
    def inner(x):
        for i in true_ids:
            if abs(x) == i:
                return True
        return False

    return np.vectorize(inner)


if __name__ == "__main__":
    args = parse_input()
    ntp = uproot.open(args.ntp)

    br_id, br_wt, br_ok = read_branches(
        ntp, args.tree, [args.base_br, args.weight_br, args.ok_br]
    )
    br_addon = read_branches_dict(ntp, args.tree, args.add_br)

    print(
        "Keep particles with these true IDs: {}".format(
            ", ".join([str(i) for i in args.id])
        )
    )
    fltr = fltr_func(args.id)
    sel = fltr(br_id)

    br_id = br_id[sel]
    br_wt = br_wt[sel]
    br_ok = br_ok[sel]

    for key in br_addon:
        br_addon[key] = br_addon[key][sel]

    # Let's sort the array by its weight, in descending order
    order = np.flipud(np.argsort(br_wt))

    name_dict = dict(zip(br_addon, args.name))

    for idx in order:
        if not br_ok[idx]:
            continue

        print("ID: {}, Weight: {:.3f}".format(br_id[idx], br_wt[idx]))
        for key in br_addon:
            print("  {}: {:.3f}".format(name_dict[key], br_addon[key][idx]))
