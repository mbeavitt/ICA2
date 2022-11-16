#!/bin/python3
# myscript.py written on 16/11/2022 by B214618. Currently this is simply a base script to build the main script on.

import configparser
import argparse
import sys
import subprocess
from urllib import request

# Initializing a parser object, to facilitate command line arguments to be passed to the program
parser = argparse.ArgumentParser(
    prog = "B214618's Protein conservation tool",
    description = "This program <content to be added>",
    epilog = "<content to be added>")

parser.add_argument('--protein', dest='protein', default="pyruvate dehydrogenase")
parser.add_argument('--group', dest='grouping', default="Ascomycota")
parser.add_argument('--database', dest='database', default="protein")

# Takes an argument of an optional infile/outfile. if none is provided, it takes the standard input. 'nargs='?'' specifies that these are optional.
#parser.add_argument('infile', nargs='?', type=argparse.FileType('r'), default=sys.stdin)
#parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'), default=sys.stdout)

args = parser.parse_args()

print(args.protein, args.grouping)

#esearch_output = subprocess.check_output("esearch -db " + "\"" + args.database + "\"" +  " -query " + "\"" + args.protein + "\"", shell=True)
with request.urlopen("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nucleotide&term=Cosmoscarta&usehistory=y") as f:
     print(f.read().decode('utf-8'))

#print(esearch_output)
