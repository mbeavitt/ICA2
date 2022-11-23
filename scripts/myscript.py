#!/bin/python3
# myscript.py written on 16/11/2022 by B214618. Currently this is simply a base script to build the main script on.

import configparser
import argparse
import sys
import subprocess
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
    groupset = list(set(seq_data["Group_ID"]))
    groupset.sort(key=int)
    print(groupset)
    group_options = groupset
    cons_list = []
    for consensus_seq in cons_output:
        cons_list.append(consensus_seq)
    with open("group_summary.txt", "w") as group_summary:
        group_summary.write("\nSUMMARY OF GROUPS:\n\n")
        for group in groupset:
            group_summary.write(f"\n\n\n\nGROUP {group} SUMMARY:\n")
            group_summary.write("\n MOST COMMON GENUSES: \n")
            group_genus = seq_data.loc[seq_data["Group_ID"] == group]["Genus"]
            group_name = seq_data.loc[seq_data["Group_ID"] == group]["Full_Name"]
            group_len = seq_data.loc[seq_data["Group_ID"] == group]["Sequence"].tolist()

            for genus, count in Counter(group_genus).most_common(10):
                group_summary.write(genus + ": " + str(count) + "\n")

            group_summary.write("\n MOST COMMON SPECIES: \n")
            for name, count in Counter(group_name).most_common(10):
                group_summary.write(name + ": " + str(count) + "\n")

            group_summary.write("\n MEDIAN SEQUENCE LENGTH: \n")
            temp_lens = []

            for sequence in group_len:
                temp_lens.append(len(sequence))
            group_summary.write(str(statistics.median(temp_lens)))
            group_summary.write("\n")
            group_summary.write("\n GROUP CONSERVED SEQUENCE: \n")
            group_cons = cons_list[int(group)].decode("utf-8")
            group_summary.write("\n".join(group_cons.split("\n")[1:]))
            group_summary.write(
                "\n\n\n-------------------------------------------------------"
            )

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
    print("Generating secondary alignments...")
    cons_output = []
    for file in group_filenames:
        subprocess.check_output(
            (
                f"clustalo --auto --force --threads={args.threads} --outfmt=msf -i"
                f" {file}.fa -o {file}.msf"
            ),
            shell=True,
        )
        outfile = subprocess.check_output(
            f"cons -sequence {file}.msf -outseq /dev/stdout",
            stderr=subprocess.DEVNULL,
            shell=True,
        )
        cons_output.append(outfile)
    return cons_output

#####################################################################################
# A function to wrangle gpc file format into a dataframe. The method for getting    #
# the protein name may seem redundant, but it's because I want to leave space       #
# for adding the option to potentially select multiple related proteins in future   #
#####################################################################################

def gpcWrangle(sequence_data):
    print("Wrangling data...")
    seq_data_list = sequence_data.split("\n")
    seq_data_list = list(filter(None, seq_data_list))

    list_of_rows = []
    # Slightly roundabout way of getting values to populate my dataframe...
    for entry in seq_data_list:
        fields = entry.split("\t")
        accession_code = fields[0]
        protein_name = args.grouping
        full_name = fields[1]
        genus = full_name.split()[0]
        sequence = fields[2]
        # Creating a list of 'dataframe rows'
        list_of_rows.append([accession_code, protein_name, genus, full_name, sequence])

    # print(list_of_rows)
    # Creating the dataframe, with relevant column names...
    seq_data = pd.DataFrame(
        list_of_rows,
        columns=["Accession", "Protein name", "Genus", "Full_Name", "Sequence"],
    )

    return seq_data

#####################################################################################
# A function to write out fasta formatted sequences of each group in the query      #
# and return a list of the group file names.                                        #
#####################################################################################

def groupFasta(cluster_dict):
    print("Generating groupwise fasta files...")
    group_filenames = []
    for key in cluster_dict.keys():
        index_list = cluster_dict[key]
        for index in index_list:
            seq_data.loc[int(index), "Group_ID"] = key
        group_fasta = ""
        group_data = seq_data.loc[seq_data["Group_ID"] == key]
        group_data.dropna()
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
        key = re.search(r"(?<=Cluster)(.*)(?=\:)", file).group(1).strip()
        value = re.search(r"(?<=index)(.*)(?=\()", file).group(1).strip()
        key_value_tuples.append((key, value))

    # Creating a dictionary of index lists
    for key, value in key_value_tuples:
        cluster_dict.setdefault(key, []).append(value)
    # Checking that there are no 'groups' of 1. If there are, they are removed.
    for key in list(cluster_dict.keys()):
        if len(cluster_dict[key]) == 1:
            cluster_dict.pop(str(key), None)
    # Renaming keys in order, using count from 0 upwards... These four lines of code
    # took me TWO DAYS to come up with. The bugs I encountered that were a result of this problem drove me MAD!!
    count = 0
    for key in list(cluster_dict.keys()):
        cluster_dict[str(count)] = cluster_dict.pop(key, None)
        count += 1
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
    print(
        'Please pick a group number or option, e.g. for "Group 1" enter "1", To finish, enter'
        ' "f":'
    )
    for option in group_options:
        possible_choices.append(option)
        print(" -- Group " + option)
    print("\n -- ALL GROUPS: a")
    print(" -- FINISH: f")
    print(" -- CLEAR CHOICES: c")
    user_choice = None
    selected = []
    while user_choice not in possible_choices:
        user_input = input("Selection: ")
        if user_input == "c":
            selected = []
            print("\nGroups in selection:")
            print("No groups selected")
            continue
        if user_input == "f":
            is_empty = len(set(selected)) == 0
            if is_empty == False:
                print("You have chosen groups...")
                for choice in set(selected):
                    print(choice + ", ", end="")
                print("\n")
            else:
                print("No groups selected")
            break
        if user_input == "a":
            selected = possible_choices
            print("\nGroups in selection:")
            print(list(set(selected)).sort())
            continue
        if str.isdigit(user_input) == True:
            input_number = int(user_input)
        else:
            print("Please choose a valid group number")
            continue
        if input_number > -1 and input_number < len(possible_choices):
            selected.append(possible_choices[input_number])
            print("\nGroups in selection:")
            print(list(set(selected).sort())
        else:
            print("Please choose a valid group number or an option")

    return selected

#### INITIALIZING AN ARGUMENT PARSER AND POPULATING IT ####

# Initializing a parser object, to facilitate command line arguments to be passed to the program
args_parser = argparse.ArgumentParser(
    prog="B214618's Protein conservation tool",
    description="This program <content to be added>",
    epilog="<content to be added>",
)

# Adding arguments...
args_parser.add_argument(
    "--protein",
    dest="protein",
    required=True,
    help='use --protein to define protein query (e.g. "Pyruvate dehydrogenase")',
)
args_parser.add_argument(
    "--group",
    dest="grouping",
    required=True,
    help='use --group to define group query (e.g. "Ascomycota" or "txid4890")',
)
args_parser.add_argument(
    "--database",
    dest="database",
    default="protein",
    help='NCBI database to query (e.g. "protein")',
)
args_parser.add_argument(
    "--winsize",
    dest="winsize",
    default=10,
    help="--winsize used to set plotcon winsize, default is 10",
)
args_parser.add_argument(
    "--cluster-size",
    dest="cluster",
    type=int,
    help=(
        "--cluster-size takes a number, and is used to determine the granularity of the"
        " resulting primary alignment. The smaller the number, the more groups will"
        " result. If --cluster-size is not defined (recommended), the clustalo program"
        " will calculate default values based on the attributes of the search."
    ),
)
args_parser.add_argument(
    "--threads",
    dest="threads",
    default=50,
    type=int,
    help=(
        "--threads takes an integer, and is used to determine how many threads to use"
        " in all processes that allow thread number choices (currently just clustalo)."
        " Default is 50."
    ),
)
args_parser.add_argument(
    "--force",
    dest="force",
    action="store_true",
    help="Use --force to overwrite files. Default is false.",
)

# assigning parsed args to variable
args = args_parser.parse_args()

#### MAIN CODE ####

# Fetching sequence query info
search_query = (
    subprocess.check_output(
        (
            f'esearch -db "{args.database}" -query "{args.protein}[Protein Name] AND'
            f' {args.grouping}[Organism] NOT partial[Properties]" | xtract -pattern'
            " ENTREZ_DIRECT -element Count"
        ),
        shell=True,
    )
    .decode("utf-8")
    .strip()
)
if search_query == 0:
    print(
        "No sequeces found for these search terms. Please try again, or broaden your"
        " search parameters"
    )
    sys.exit()
print(
    "\nThe number of sequences in your search are: "
    + str(search_query)
    + "."
    + "\n Would you like to continue...? Please note sequence entries >1000 may take a"
    " while to process."
)
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
sequence_data = subprocess.check_output(
    (
        f'esearch -db "{args.database}" -query "{args.protein}[Protein Name] AND'
        f' {args.grouping}[Organism] NOT partial[Properties]" | efetch -format gpc |'
        " xtract -pattern INSDSeq -element INSDSeq_accession-version INSDSeq_organism"
        " INSDSeq_sequence"
    ),
    shell=True,
).decode("utf-8")

#### ---- DEBUGGING SECTION - to save time on NCBI queries, saves to file instead, and subsequently reads from that. ---- ####
#         run once with first two lines uncommented, then comment out subprocess above and uncomment third line.
#
# with open("sequence_data_save.txt", "w") as seq_data_save:
#    seq_data_save.write(sequence_data)
# sequence_data = open("sequence_data_save.txt", "r").read()
#
#
#### ----/DEBUGGING SECTION ---- ####

seq_data = gpcWrangle(sequence_data)

fasta_string = ""
for accession, sequence in zip(seq_data.Accession, seq_data.Sequence):
    fasta_string = fasta_string + ">" + accession + "\n" + sequence + "\n"

# Writing fasta file for primary MSA
with open(f"fasta_formatted.fa", "w") as fasta_formatted_file:
    fasta_formatted_file.write(fasta_string)

# Running the primary MSA
print("Running primary MSA...")
subprocess.run(
    (
        "clustalo --force --auto "
        f" --threads={args.threads} --clustering-out=clusterfile.txt --outfmt=msf -i"
        " fasta_formatted.fa -o align.msf"
    ),
    shell=True,
)
clusterfile = open("clusterfile.txt", "r").read().split("\n")
clusterfile = list(filter(None, clusterfile))

# Input > Process > Output... Repeat!
cluster_dict = clusterIndexer(clusterfile)
group_filenames = groupFasta(cluster_dict)
cons_output = groupwiseMSA(group_filenames)
seq_data = seq_data.dropna()
cons_summary, group_options = groupDisplay(seq_data, cons_output)
print(cons_summary)
user_selection = groupChoose(group_options)
