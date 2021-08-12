#!/usr/bin/env python
#
# Author: Zishuo Yang, Yipeng Sun
# License: GPLv2
# Based on:
#   https://github.com/ZishuoYang/my-hammer-reweighting/blob/master/plot_ratio.py
# Last Change: Fri Aug 13, 2021 at 01:42 AM +0200

import ROOT as rt
rt.PyConfig.IgnoreCommandLineOptions = True  # Don't hijack argparse!
rt.PyConfig.DisableRootLogon = True  # Don't read .rootlogon.py

from argparse import ArgumentParser
from os.path import basename, splitext


#################################
# Command line arguments parser #
#################################

def parse_input():
    parser = ArgumentParser(description='''
Plot fit variables of R(D(*)) with and without FF reweighting.''')

    parser.add_argument('-n', '--ntuple',
                        required=True,
                        help='''
specify path to ntuple that contains HAMMER weights.''')

    parser.add_argument('-t', '--tree',
                        required=True,
                        help='''
specify tree name in HAMMER-weight ntuple.''')

    parser.add_argument('-o', '--output-path',
                        default='./gen',
                        help='''
specify output path for plots.''')

    parser.add_argument('--ff-weight',
                        default='w_ff',
                        help='''
specify branch name of the FF weight.''')

    parser.add_argument('--cuts',
                        default='flag_ham_ok',
                        help='''
specify cuts to the plots.''')

    parser.add_argument('--vars',
                        nargs='+',
                        default=['q2_true'],
                        help='''
specify variables to plot.''')

    parser.add_argument('--bin-ranges',
                        nargs='+',
                        default=['(15,0,12)', '(15,-0.2,10)', '(15,0.,2.1)'],
                        help='''
specify number of bins and x ranges.''')

    return parser.parse_args()


###########
# Helpers #
###########

def get_filename(path):
    return basename(splitext(path)[0])


def plot_ratio(tree, output_path, ntp_name,
               var, bin_range, cuts, weight):
    canvas = rt.TCanvas('canvas', 'canvas')
    plot_suffix = get_filename(ntp_name)

    tree.Draw('{}>>h1{}'.format(var, bin_range), cuts, "goff")
    h1 = rt.gDirectory.Get('h1')
    h1.SetMarkerColor(rt.kBlue)
    h1.SetLineColor(rt.kBlue)

    tree.Draw('{}>>h2{}'.format(var, bin_range), '{}'.format(weight), "goff")
    h2 = rt.gDirectory.Get('h2')
    h2.SetMarkerColor(rt.kRed)
    h2.SetLineColor(rt.kRed)

    rp = rt.TRatioPlot(h1, h2, 'divsym')
    rp.Draw()

    canvas.Update()
    canvas.Print('{}/{}_{}.png'.format(output_path, plot_suffix, var))


if __name__ == '__main__':
    args = parse_input()

    ntuple = rt.TFile(args.ntuple)
    tree = ntuple.Get(args.tree)

    rt.gROOT.SetBatch(rt.kTRUE)  # Don't output anything on screen

    for var, bin_range in zip(args.vars, args.bin_ranges):
        plot_ratio(tree, args.output_path, args.ntuple,
                   var, bin_range, args.cuts, args.ff_weight)
