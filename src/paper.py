#!/usr/bin/env python3

"""
Python Align Pair-End Reads (PAPER)
"""

import argparse
from classes import Aligner

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Python Align Pair-End Reads")
    parser.add_argument('-a', dest='align', action='store_true', help='Align pair-end reads')
    parser.add_argument('-i', dest='data_dir', default='', help='Directory of sequence files')
    parser.add_argument('-o', dest='out_dir', default='', help='Directory for results')
    parser.add_argument('-c', dest='cpu', default=None, help='Number of processors')
    parser.add_argument('-m', dest='mem', default=None, help='Amount of memory in G')

    args = parser.parse_args()

    if args.align:
        aligner = Aligner(target_dir=args.data_dir, out_dir=args.out_dir, mem=args.mem, cpu=args.cpu)
        aligner.trinity()
