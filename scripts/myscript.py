#!/bin/python3
# myscript.py written on 16/11/2022 by B214618. Currently this is simply a base script to build the main script on.

import configparser
import argparse
import sys
import subprocess
import requests
import re


# Initializing a parser object, to facilitate command line arguments to be passed to the program
args_parser = argparse.ArgumentParser(
    prog = "B214618's Protein conservation tool",
    description = "This program <content to be added>",
    epilog = "<content to be added>")

args_parser.add_argument('--protein', dest='protein', default="pyruvate dehydrogenase")
args_parser.add_argument('--group', dest='grouping', default="Ascomycota")
args_parser.add_argument('--database', dest='database', default="protein")

# Takes an argument of an optional infile/outfile. if none is provided, it takes the standard input. 'nargs='?'' specifies that these are optional.
#parser.add_argument('infile', nargs='?', type=argparse.FileType('r'), default=sys.stdin)
#parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'), default=sys.stdout)

args = args_parser.parse_args()
print(args.protein, args.grouping)

#esearch_output = subprocess.check_output("esearch -db " + "\"" + args.database + "\"" +  " -query " + "\"" + args.protein + "\"", shell=True)
with requests.get("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nucleotide&term=Cosmoscarta&usehistory=y") as f:
    xml_file_byte = f.content
    xml_file = xml_file_byte.decode('utf-8')
    webenv_match = re.search('(?<=<WebEnv>)(.*)(?=<\/WebEnv>)', xml_file)
    querykey_match = re.search('(?<=<QueryKey>)(.*)(?=<\/QueryKey>)', xml_file)
    if webenv_match:
        web_env = webenv_match.group(1)
        print('WebEnv key:' + web_env)
    if querykey_match:
        query_key = querykey_match.group(1)
        print('Query key:' + query_key)
#print(esearch_output)
