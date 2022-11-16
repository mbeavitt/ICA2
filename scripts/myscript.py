#!/bin/python3
# myscript.py written on 16/11/2022 by B214618. Currently this is simply a base script to build the main script on.

import configparser
import argparse
import sys

# Initializing a parser object, to facilitate command line arguments to be passed to the program
parser = argparse.ArgumentParser(
    prog = "B214618's Protein conservation tool",
    description = "This program <content to be added>",
    epilog = "<content to be added>")

# Takes an argument of an optional infile/outfile. if none is provided, it takes the standard input. 'nargs='?'' specifies that these are optional.
parser.add_argument('infile', nargs='?', type=argparse.FileType('r'), default=sys.stdin)
parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'), default=sys.stdout)

args = parser.parse_args()

print(args.infile, args.outfile)

for line in args.infile:
    args.outfile.write(line)
