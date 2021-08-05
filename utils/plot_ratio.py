#!/usr/bin/env python
#
# Author: Zishuo Yang, Yipeng Sun
# License: GPLv2
# Based on:
#   https://github.com/ZishuoYang/my-hammer-reweighting/blob/master/plot_ratio.py
# Last Change: Thu Aug 05, 2021 at 08:51 PM +0200

import ROOT as rt

from argparse import ArgumentParser


#################################
# Command line arguments parser #
#################################

def parse_input():
    parser = ArgumentParser(description='''
Plot fit variables of R(D(*)) with and without FF reweighting.''')

    parser.add_argument('-o', '--output-path',
                        nargs='?',
                        default='./gen',
                        help='''
specify output path for plots.''')

    parser.add_argument('-d', '--data-ntuple',
                        nargs='?',
                        required=True,
                        help='''
specify path to ntuple that contains fit variables.''')

    parser.add_argument('-t', '--data-tree',
                        nargs='?',
                        required=True,
                        help='''
specify tree name in fit-variable ntuple.''')

    parser.add_argument('--ff-weight',
                        nargs='?',
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
                        default=['(80,3,11)', '(80,-0.2,10)', '(80,0.,2.1)'],
                        help='''
specify number of bins and x ranges.''')

    return parser.parse_args()


###########
# Helpers #
###########

def plot_ratio(tree, output_path,
               var, weight, title,
               bin_range, cuts):
    # rt.gStyle.SetOptStat(0)
    canvas = rt.TCanvas('canvas', 'A ratio plot')

    tree.Draw('{}>>h1{}'.format(var, bin_range), cuts, 'goff')
    h1 = rt.gDirectory.Get('h1')
    h1.SetMarkerColor(rt.kBlue)
    h1.SetLineColor(rt.kBlue)

    tree.Draw('{}>>h2{}'.format(var, bin_range), '({})*{}'.format(cuts, weight)
              , 'goff')
    h2 = rt.gDirectory.Get('h2')
    h2.SetMarkerColor(rt.kRed)
    h2.SetLineColor(rt.kRed)

    rp = rt.TRatioPlot(h1, h2, 'divsym')
    rp.Draw()

    canvas.Update()

    canvas.Print('{}/{}.png'.format(output_path, title))


if __name__ == '__main__':
    args = parse_input()

    data_ntuple = rt.TFile(args.data_ntuple)
    data_tree = data_ntuple.Get(args.data_tree)

    rt.gROOT.SetBatch(rt.kTRUE)  # Don't output anything on screen

    for var, bin_range in zip(args.vars, args.bin_ranges):
        plot_ratio(data_tree, args.output_path,
                   var, args.ff_weight, var,
                   bin_range, args.cuts)
