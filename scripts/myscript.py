#!/bin/python3
# myscript.py written on 16/11/2022 by B214618. Currently this is simply a base script to build the main script on.

import configparser
import argparse
import sys
import subprocess
#import requests
import re
import pandas as pd
import pprint
import os
import statistics
from collections import Counter

#####################################################################################
# A function to write out the group summary information to a file, then             #
# return the contents of the file. The information is:                              #
#    - Most common genuses                                                          #
#    - most common species                                                          #
#    - median sequence length                                                       #
#    - conserved sequence of group                                                  #
#####################################################################################

def groupDisplay(seq_data, cons_output):
    groupset = list(set(seq_data['Group_ID']))
    groupset.sort()
    group_options = groupset
    cons_list = []
    for consensus_seq in cons_output:
        cons_list.append(consensus_seq)
    with open("group_summary.txt", "w") as group_summary:
        group_summary.write("\nSUMMARY OF GROUPS:\n\n")
        for group in groupset:
            group_summary.write(f"\n\nGROUP {group} SUMMARY:\n")
            group_summary.write("\n MOST COMMON GENUSES: \n")
            group_genus = seq_data.loc[seq_data['Group_ID'] == group]['Genus']
            group_name = seq_data.loc[seq_data['Group_ID'] == group]['Full_Name']
            group_len = seq_data.loc[seq_data['Group_ID'] == group]['Sequence'].tolist()

            for genus, count in Counter(group_genus).most_common(10):
                group_summary.write(genus + ": " + str(count) + '\n')

            group_summary.write("\n MOST COMMON SPECIES: \n")
            for name, count in Counter(group_name).most_common(10):
                group_summary.write(name + ": " + str(count) + '\n')

            group_summary.write("\n MEDIAN SEQUENCE LENGTH: \n")
            temp_lens = []

            for sequence in group_len:
                temp_lens.append(len(sequence))
            group_summary.write(str(statistics.median(temp_lens)))
            group_summary.write('\n')
            group_summary.write("\n GROUP CONSERVED SEQUENCE: \n")
            group_cons = cons_list[int(group)].decode("utf-8")
            group_summary.write('\n'.join(group_cons.split('\n')[1:]))
            group_summary.write('\n\n\n-------------------------------------------------------')

    cons_summary = open("group_summary.txt", "r").read()

    return cons_summary, group_options

#####################################################################################
# A basic function for creating groupwise MSAs and then getting a consensus sequence#
# for them. I guess I don't really need a function for this, but it feels 'tidy'    #
# the cons program was a bit annoying to work with. Its intro messages were being   #
# routed to stderr and needed to be silenced, and you couldn't have a regular       #
# stdout if no out file was passed, necessitating writing of yet another file       #
#####################################################################################


def groupwiseMSA(group_filenames):
    print("Generating secondary alignments")
    cons_output = []
    for file in group_filenames:
        subprocess.check_output(f"clustalo --auto --force --threads={args.threads} --outfmt=msf -i {file}.fa -o {file}.msf", shell=True)
        outfile = subprocess.check_output(f"cons -sequence {file}.msf -outseq /dev/stdout", stderr=subprocess.DEVNULL, shell=True)
        cons_output.append(outfile)
    return cons_output



#####################################################################################
# A function to wrangle gpc file format into a dataframe. The method for getting    #
# the protein name may seem redundant, but it's because I want to leave space       #
# for adding the option to potentially select multiple related proteins in future   #
#####################################################################################

def gpcWrangle(sequence_data):
    print('Wrangling data...')
    seq_data_list = sequence_data.split('\n')
    seq_data_list = list(filter(None, seq_data_list))

    list_of_rows=[]
    # Slightly roundabout way of getting values to populate my dataframe...
    for entry in seq_data_list:
        fields = entry.split('\t')
        accession_code = fields[0]
        protein_name = args.grouping
        full_name = fields[1]
        genus = full_name.split()[0]
        sequence = fields[2]
    # Creating a list of 'dataframe rows'
        list_of_rows.append([accession_code, protein_name, genus, full_name, sequence])

    #print(list_of_rows)
    # Creating the dataframe, with relevant column names...
    seq_data = pd.DataFrame(list_of_rows, columns = ['Accession', 'Protein name', 'Genus', 'Full_Name', 'Sequence'])

    return seq_data


#####################################################################################
# A function to write out fasta formatted sequences of each group in the query      #
# and return a list of the group file names.                                        #
#####################################################################################

def groupFasta(cluster_dict):
    print('Generating groupwise fasta files...')
    group_filenames = []
    for key in cluster_dict.keys():
        index_list = cluster_dict[key]
        for index in index_list:
            seq_data.loc[int(index), 'Group_ID'] = key
        group_fasta = ""
        group_data = seq_data.loc[seq_data['Group_ID'] == key]
        for accession, sequence in zip(group_data.Accession, group_data.Sequence):
            group_fasta = group_fasta + ">" + accession + "\n" + sequence + "\n"
            with open(f"group_{key}.fa", "w") as groupalign:
                groupalign.write(group_fasta)
        group_filenames.append(f"group_{key}")
    return group_filenames


#####################################################################################
# A function for creating a dictionary where each key corresponds to a multiple     #
# alignment cluster, and each "value" is a list of index positions in a             #
# dataframe that can be used to source data about the members of the cluster        #
#####################################################################################

def clusterIndexer(clusterfile):
    print("Building groupings...")
    key_value_tuples = []
    cluster_dict = {}
    for file in clusterfile:
    # Extracting the bits I want using regex
        key = re.search(r'(?<=Cluster)(.*)(?=\:)', file).group(1).strip()
        value = re.search(r'(?<=index)(.*)(?=\()', file).group(1).strip()
        key_value_tuples.append((key, value))

    # Creating a dictionary of index lists
    for key, value in key_value_tuples:
        cluster_dict.setdefault(key, []).append(value)
    # Checking that there are no 'groups' of 1. If there are, they are removed.
    if len(cluster_dict[key]) == 1:
        cluster_dict.pop(key, None)
    if len(cluster_dict.keys()) == 1:
        print("Skip to final stage <to be coded>")
    if len(cluster_dict.keys()) == 0:
        print("Something has gone terribly wrong. Try expanding your search...?")
        sys.exit()
    return cluster_dict


#####################################################################################
# A function for choosing groups formed by clustal omega. The idea is that this     #
# removes the dependence on species choice/sequence length for narrowing down the   #
# comparison, and lets the user compare based on groups of similar sequences        #
# instead. This should allow much finer control over the process, and give insights #
# into similarity *across* species and genus as well as within.                     #
#####################################################################################

def groupChoose(group_options):

    possible_choices = []
    print('Please pick a group number, e.g. for \"Group 1\" enter \"1\". To finish, enter \"f\":')
    for option in group_options:
        possible_choices.append(option)
        print(" -- Group " + option)
    print("\n -- ALL GROUPS: a")
    print(" -- FINISH: f")
    print(" -- CLEAR CHOICES: c")
    user_choice = None
    selected = []
    while user_choice not in possible_choices:
        user_input = input('Selection: ')
        if user_input == "c":
            selected = []
            print("\nGroups in selection:")
            print("No groups selected")
            continue
        if user_input == "f":
            is_empty = (len(set(selected)) == 0)
            if is_empty == False:
                print(set(selected))
            else:
                print("No groups selected")
            break
        if user_input == "a":
            selected = possible_choices
            print("\nGroups in selection:")
            print(set(selected))
            continue
        if str.isdigit(user_input) == True:
            input_number = int(user_input)
        else:
            print('Please choose a valid group number')
            continue
        if input_number > -1 and input_number < len(possible_choices):
            selected.append(possible_choices[input_number])
            print("\nGroups in selection:")
            print(set(selected))
        else:
            print('Please choose a valid group number')

    return selected


#### INITIALIZING AN ARGUMENT PARSER AND POPULATING IT ####

# Initializing a parser object, to facilitate command line arguments to be passed to the program
args_parser = argparse.ArgumentParser(
    prog = "B214618's Protein conservation tool",
    description = "This program <content to be added>",
    epilog = "<content to be added>")

# Adding arguments...
args_parser.add_argument('--protein', dest='protein', required=True, help = "use --protein to define protein query (e.g. \"Pyruvate dehydrogenase\")")
args_parser.add_argument('--group', dest='grouping', required=True, help = "use --group to define group query (e.g. \"Ascomycota\" or \"txid4890\")")
args_parser.add_argument('--database', dest='database', default="protein", help = "NCBI database to query (e.g. \"protein\")")
args_parser.add_argument('--winsize', dest='winsize', default=10, help = "--winsize used to set plotcon winsize, default is 10")
args_parser.add_argument('--cluster-size', dest='cluster', type=int, help = "--cluster-size takes a number, and is used to determine the granularity of the resulting primary alignment. The smaller the number, the more groups will result. If --cluster-size is not defined (recommended), the clustalo program will calculate default values based on the attributes of the search.")
args_parser.add_argument('--threads', dest='threads', default=50, type=int, help = "--threads takes an integer, and is used to determine how many threads to use in all processes that allow thread number choices (currently just clustalo). Default is 50.")
args_parser.add_argument('--save-group-summary', dest='savesum', action='store_true', help = "--save-group-summary will save the group summary page in this program's output to a .txt file called group_summary.txt.")
args_parser.add_argument('--force', dest='force', action='store_true', help = "Use --force to overwrite files. Default is false.")


# assigning parsed args to variable
args = args_parser.parse_args()


#### MAIN CODE ####


# Fetching sequence query info
search_query = subprocess.check_output(f"esearch -db \"{args.database}\" -query \"{args.protein}[Protein Name] AND {args.grouping}[Organism] NOT partial[Properties]\" | xtract -pattern ENTREZ_DIRECT -element Count", shell=True).decode('utf-8').strip()
if search_query == 0:
    print("No sequeces found for these search terms. Please try again, or broaden your search parameters")
    sys.exit()
print("\nThe number of sequences in your search are: " + str(search_query) + "." + "\n Would you like to continue...? Please note sequence entries >1000 may take a while to process.")
continue_ornot = None
while continue_ornot not in {"y", "n"}:
    continue_ornot = input("Please enter y/n:")

if continue_ornot == "y":
    print("Processing...")
    pass
else:
    sys.exit()

# Fetch sequence data using gpc format
print("Fetching sequences...")
sequence_data = subprocess.check_output(f"esearch -db \"{args.database}\" -query \"{args.protein}[Protein Name] AND {args.grouping}[Organism] NOT partial[Properties]\" | efetch -format gpc | xtract -pattern INSDSeq -element INSDSeq_accession-version INSDSeq_organism INSDSeq_sequence", shell = True).decode('utf-8')

seq_data = gpcWrangle(sequence_data)

fasta_string = ""
for accession, sequence in zip(seq_data.Accession, seq_data.Sequence):
    fasta_string = fasta_string + ">" + accession + "\n" + sequence + "\n"

# Writing fasta file for primary MSA
with open(f"fasta_formatted.fa", "w") as fasta_formatted_file:
    fasta_formatted_file.write(fasta_string)


# Running the primary MSA
print("Running primary MSA...")
subprocess.run(f"clustalo --force --auto  --threads={args.threads} --clustering-out=clusterfile.txt --outfmt=msf -i fasta_formatted.fa -o align.msf", shell=True)
clusterfile = open("clusterfile.txt", "r").read().split('\n')
clusterfile =  list(filter(None, clusterfile))


# Input > Process > Output... Repeat!
cluster_dict = clusterIndexer(clusterfile)
group_filenames = groupFasta(cluster_dict)
cons_output = groupwiseMSA(group_filenames)
cons_summary, group_options = groupDisplay(seq_data, cons_output)
print(cons_summary)
user_selection = groupChoose(group_options)

#print(key_value_tuples)
#print(seq_data)
#test_index_list = cluster_dict['0']
#print(seq_data.iloc[test_index_list])

# An attempt to design a grouping system around percent identity - abandoned in favour of distance grouping
#pi_dist_mat = open("distmat.txt", "r").read().split('\n')
#pi_dist_mat.pop(0)
#pi_dist_mat.pop(-1)
#colnames = []
#count = 0
#for i in pi_dist_mat:
#    pi_dist_mat[count] = i.split()
#    colnames.append(pi_dist_mat[count][0])
#    pi_dist_mat[count].pop(0)
#    count += 1
#pi_matrix_df = pd.DataFrame(pi_dist_mat, columns = [colnames], index = [colnames])

#print("running plotcon...")
#subprocess.run(f"plotcon align.msf -winsize {args.winsize} -graph png", shell=True)

# Takes an argument of an optional infile/outfile. if none is provided, it takes the standard input. 'nargs='?'' specifies that these are optional.
#parser.add_argument('infile', nargs='?', type=argparse.FileType('r'), default=sys.stdin)
#parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'), default=sys.stdout)

# It felt like cheating to just be passing commands to the shell, so I decided to work out how to send and receive requests directly in python instead!
#print(f"http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db={database}&term={grouping}[Organism]+{protein}[Protein+Name]&usehistory=y")
#with requests.get(f"http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db={database}&term={grouping}[Organism]+AND{protein}[Protein+Name]&usehistory=y") as f:
# INTERESTING URL FORMAT: https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?retmax=0&usehistory=y&db=protein&term=pyruvate%20dehydrogenase%5BPROT%5D%20AND%20Ascomycota%20%5BORGANISM%5D%20NOT%20partial%20%5BProperties%5D
#    xml_file = f.content.decode('utf-8')
#    webenv_match = re.search('(?<=<WebEnv>)(.*)(?=<\/WebEnv>)', xml_file)
#    querykey_match = re.search('(?<=<QueryKey>)(.*)(?=<\/QueryKey>)', xml_file)
#    if webenv_match:
#        web_env = webenv_match.group(1)
#        print('WebEnv key:' + web_env)
#    if querykey_match:
#        query_key = querykey_match.group(1)
#        print('Query key:' + query_key)

#print(f"http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db={database}&query_key={query_key}&WebEnv={web_env}&rettype=uilist&retmode=text")
#with requests.get(f"http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db={database}&query_key={query_key}&WebEnv={web_env}&rettype=uilist&retmode=text") as f:
#    my_file = f.content.decode('utf-8')
#    print(my_file)
