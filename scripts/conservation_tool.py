#!/bin/python3
# myscript.py written on 16/11/2022 by B214618.

#### ---- IMPORTING MODULES ---- ####

import configparser
import argparse
import sys
import subprocess
import re
import pandas as pd
import os
import statistics
from collections import Counter
from shutil import rmtree
import numpy as np
#import matplotlib.pyplot as plt
#from matplotlib.backends.backend_pdf import PdfPages


# A temporary solution to display the dataframe at the end - in the "final version" I hope to output a PDF report.
pd.set_option("display.max_rows", 1000)

#### ---- FUNCTIONS ---- ####

#####################################################################################
# An interactive prompt that guides the user through removing files from previous   #
# analysis sessions. I did think about implementing a method for the user to create #
# new folders with unique names, but that could quickly become messy, and it's not  #
# something I've seen other command line programs do before. --force tends to be    #
# more common in my limited experience.                                             #
#####################################################################################


def checkDirs(dirs_list):
    # Command line argument will bypass the interactive elements
    if args.force == True:
        for dir in dirs_list:
            # Checking nothing exists before making a directory...
            if not os.path.exists(dir):
                os.mkdir(dir)
            # A previous directory is found! Removing it, and making a new one.
            else:
                rmtree(dir)
                os.mkdir(dir)
        return

    continue_pipeline = True
    for dir in dirs_list:
        if os.path.exists(dir):
            print(
                "\nYou have already run an analysis. Please save your data, remove the"
                " folders, and try again.\nAlternatively, would you like to overwrite"
                " your files...?"
            )
            continue_pipeline = False
            break

    if continue_pipeline == True:
        for dir in dirs_list:
            os.mkdir(dir)
    else:
        continue_ornot = None
        while continue_ornot not in {"y", "n"}:
            continue_ornot = input("Overwrite?  y/n:")

        if continue_ornot == "y":
            print("Files will be overwritten.\nSearching NCBI...")
            for dir in dirs_list:
                if not os.path.exists(dir):
                    os.mkdir(dir)
                else:
                    rmtree(dir)
                    os.mkdir(dir)
        else:
            sys.exit()


#####################################################################################
# A function to wrangle gpc file format into a dataframe. The method for getting    #
# the protein name may seem redundant, but it's because I want to leave space       #
# for adding the option to potentially select multiple related proteins in future   #
# NOTE: This functionality has been added now (per se) so this needs fixed          #
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
        # This next line needs to be changed ASAP. It should extract the protein name from the record, but I don't know how to do that without mistakes yet. This is fine for one protein, but when the general protein search argument is specified, the output will not be correct.
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
    # Placeholder for potential to do something different if there's only one group. Didn't get around to coding this unfortunately
    if len(cluster_dict.keys()) == 1:
        print("")
    if len(cluster_dict.keys()) == 0:
        print("Something has gone terribly wrong. Try expanding your search...?")
        sys.exit()
    return cluster_dict


#####################################################################################
# A function to write out fasta formatted sequences of each group in the query      #
# and return a list of the group file names.                                        #
#####################################################################################


def groupFasta(cluster_dict):
    print("Generating groupwise fasta files...")
    group_filenames = []
    for key in cluster_dict.keys():
        # Indexing the dataframe and writing groupwise multiple fasta files for the group MSA stage
        index_list = cluster_dict[key]
        for index in index_list:
            seq_data.loc[int(index), "Group_ID"] = key
        group_fasta = ""
        group_data = seq_data.loc[seq_data["Group_ID"] == key]
        group_data.dropna()
        for accession, sequence in zip(group_data.Accession, group_data.Sequence):
            group_fasta = group_fasta + ">" + accession + "\n" + sequence + "\n"
            with open(f"{fasta_path}group_{key}.fa", "w") as groupalign:
                groupalign.write(group_fasta)
        group_filenames.append(f"group_{key}")
    return group_filenames


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
    # A slightly messy way of accomodating different arguments passed. I'm sure there's a better way of doing this.
    # This handles custom cluster sizes
    if args.cluster != None:
        for file in group_filenames:
            subprocess.check_output(
                (
                    "clustalo "
                    "--auto "
                    "--force "
                    f"--threads={args.threads} "
                    f"--cluster-size={args.cluster}"
                    " --outfmt=msf "
                    f"-i {fasta_path}{file}.fa "
                    f"-o {msa_path}{file}.msf"
                ),
                shell=True,
            )
            outfile = subprocess.check_output(
                # /dev/stdout used to route output to variable.
                f"cons -sequence {msa_path}{file}.msf -outseq /dev/stdout",
                stderr=subprocess.DEVNULL,
                shell=True,
            )
            cons_output.append(outfile)
        return cons_output
    else:
        # If no custom cluster size is specified...
        for file in group_filenames:
            subprocess.check_output(
                (
                    "clustalo "
                    "--auto "
                    "--force "
                    f"--threads={args.threads} "
                    "--outfmt=msf "
                    f"-i {fasta_path}{file}.fa "
                    f"-o {msa_path}{file}.msf"
                ),
                shell=True,
            )
            outfile = subprocess.check_output(
                # /dev/stdout used to route output to variable.
                f"cons -sequence {msa_path}{file}.msf -outseq /dev/stdout",
                stderr=subprocess.DEVNULL,
                shell=True,
            )
            cons_output.append(outfile)
        return cons_output


#####################################################################################
# A function to write out the group summary information to a file, then             #
# return the contents of the file. The information is:                              #
#    - Most common genuses                                                          #
#    - most common species                                                          #
#    - median sequence length                                                       #
#    - conserved sequence of group                                                  #
#####################################################################################


def groupDisplay(seq_data, cons_output):
    if args.nogrouping:
        pass
    else:
        # Getting a set of all the possible groups and sorting them
        groupset = list(set(seq_data["Group_ID"]))
        groupset.sort(key=int)
        group_options = groupset
        cons_list = []

        # Do I really need to do this??? Isn't it a list already...???? Leaving this in because I don't have time to test.
        for consensus_seq in cons_output:
            cons_list.append(consensus_seq)

    # Writing summaries! All pretty basic, just using a context manager to write a bunch of stuff to a text file
    # really wanted to do pdfs, but ran out of time to learn how to do that unfortunately. The Counter module came in clutch here!
    if args.nogrouping:
        with open(f"{summary_path}nogrouping_summary.txt", "w") as group_summary:
            group_summary.write("\nSUMMARY OF SEARCH ANALYSIS (NO GROUPING):\n\n")
            group_summary.write("\nMOST COMMON GENUSES: \n")
            for genus, count in Counter(seq_data["Genus"]).most_common(10):
                group_summary.write(genus + ": " + str(count) + "\n")

            group_summary.write("\nMOST COMMON SPECIES: \n")
            for name, count in Counter(seq_data["Full_Name"]).most_common(10):
                group_summary.write(name + ": " + str(count) + "\n")

            group_summary.write("\nMEDIAN SEQUENCE LENGTH: \n")
            temp_lens = []
            for sequence in seq_data["Sequence"].tolist():
                temp_lens.append(len(sequence))
            group_summary.write(str(statistics.median(temp_lens)))
            group_summary.write("\n")
            group_summary.write("\nCONSERVED SEQUENCE: \n")
            group_cons = cons_output.decode("utf-8")
            group_summary.write("\n".join(group_cons.split("\n")[1:]))

        cons_summary = open(f"{summary_path}nogrouping_summary.txt", "r").read()

        return cons_summary

    else:
        with open(f"{summary_path}group_summary.txt", "w") as group_summary:
            group_summary.write("\nSUMMARY OF GROUPS:\n\n")
            for group in groupset:
                group_summary.write(f"\n\n\n\nGROUP {group} SUMMARY:\n")
                group_summary.write("\nMOST COMMON GENUSES: \n")
                group_genus = seq_data.loc[seq_data["Group_ID"] == group]["Genus"]
                group_name = seq_data.loc[seq_data["Group_ID"] == group]["Full_Name"]
                group_len = seq_data.loc[seq_data["Group_ID"] == group][
                    "Sequence"
                ].tolist()

                for genus, count in Counter(group_genus).most_common(10):
                    group_summary.write(genus + ": " + str(count) + "\n")

                group_summary.write("\nMOST COMMON SPECIES: \n")
                for name, count in Counter(group_name).most_common(10):
                    group_summary.write(name + ": " + str(count) + "\n")

                group_summary.write("\nMEDIAN SEQUENCE LENGTH: \n")
                temp_lens = []

                for sequence in group_len:
                    temp_lens.append(len(sequence))
                group_summary.write(str(statistics.median(temp_lens)))
                group_summary.write("\n")
                group_summary.write("\nGROUP CONSERVED SEQUENCE: \n")
                group_cons = cons_list[int(group)].decode("utf-8")
                group_summary.write("\n".join(group_cons.split("\n")[1:]))
                group_summary.write(
                    "\n\n\n-------------------------------------------------------"
                )

        cons_summary = open(f"{summary_path}group_summary.txt", "r").read()

        return cons_summary, group_options


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
        'Please pick a group number or option, e.g. for "Group 1" enter "1", To finish,'
        ' enter "f":'
    )
    # Getting options and showing the user
    for option in group_options:
        possible_choices.append(option)
        print(" -- Group " + option)
    print("\n -- ALL GROUPS: a")
    print(" -- FINISH: f")
    print(" -- CLEAR CHOICES: c")
    print(" -- VIEW CONSERVATION PLOTS: p\n")
    user_choice = "None"
    selected = []
    # Getting user input, handling menu choices.
    # This bit allows the user to only choose from the list of possible picks.
    while user_choice not in possible_choices:
        user_input = input("Selection: ")
        if user_input == "c":
            selected = []
            print("\nGroups in selection:")
            print("No groups selected")
            continue
        # Sending the user to another menu for plotcon
        if user_input == "p":
            chooseConsPlot(group_options)
            print(
                'Please pick a group number or option, e.g. for "Group 1" enter "1", To'
                ' finish, enter "f":'
            )
            for option in group_options:
                print(" -- Group " + option)
            print("\n -- ALL GROUPS: a")
            print(" -- FINISH: f")
            print(" -- CLEAR CHOICES: c")
            print(" -- VIEW CONSERVATION PLOTS: p\n")
            continue
        # Exit statement
        if user_input == "f":
            is_empty = len(set(selected)) == 0
            if is_empty == False:
                print("You have chosen group(s)...")
                for choice in set(selected):
                    print(choice + ", ", end="")
                print("\n")
            else:
                print("No groups selected")
                os.exit()
            break
        if user_input == "a":
            selected = possible_choices
            print("\nGroups in selection:")
            print(set(selected))
            continue
        if str.isdigit(user_input) == True:
            input_number = int(user_input)
        else:
            print("Please choose a valid group number")
            continue
        if input_number > -1 and input_number < len(possible_choices):
            selected.append(possible_choices[input_number])
            print("\nGroups in selection:")
            print(set(selected))
        # Invalid option check
        else:
            print("Please choose a valid group number or an option")

    return selected


#####################################################################################
# A function for allowing the user to view conservation plots of the groups of      #
# interest, one at a time. Ideally I'd like to find some other more elegant way     #
# of doing this, but for now this will suffice.                                     #
#####################################################################################


def chooseConsPlot(group_options):
    choices = []
    print(
        "\nPlease pick a group number to produce and save a conservation plot for that"
        ' group. To go back to the previous screen, enter "q":\n Please note, in order'
        " to continue after producing a plot, you must close the plot window."
    )
    # Setting possible options
    for option in group_options:
        choices.append(option)
        print(" -- Group " + option)
    print("\n -- GO BACK: q")
    user_choice = None
    while user_choice not in choices:
        # Taking user input
        user_input = None
        user_input = input("Selection: ")
        if user_input == "q":
            print("\n\n")
            break
        # Checking user input and setting it as a variable
        if str.isdigit(user_input) == True:
            input_number = int(user_input)
        else:
            print("Please choose a valid group number")
            continue
        # Calling plotcon commands from the comand line
        if input_number > -1 and input_number < len(choices):
            subprocess.call(
                (
                    f"plotcon -sequences {msa_path}group_{input_number}.msf "
                    f"-winsize {args.winsize} "
                    "-graph x11"
                ),
                shell=True,
            )
            subprocess.call(
                (
                    f"plotcon -sequences {msa_path}group_{input_number}.msf "
                    f"-winsize {args.winsize} "
                    "-graph png"
                ),
                shell=True,
            )
            print("Plot saved")


#####################################################################################
# A function for searching the local prosite database using the groups acquired     #
# from the user's input.                                                            #
#####################################################################################


def prositeGroupSearch(user_selection=""):
    if args.nogrouping:
        for accession, sequence in zip(seq_data["Accession"], seq_data["Sequence"]):
            # Making individual fasta files as variables, then writing them to a file.
            # The patmatmotif analysis is carried out once every loop on the single file.
            fasta_formatted = ">" + accession + "\n" + sequence
            with open(f"{prosite_path}{accession}.fa", "w") as fasta_format_file:
                fasta_format_file.write(fasta_formatted)
            subprocess.run(
                (
                    f"patmatmotifs -sequence {prosite_path}{accession}.fa "
                    f"-outfile {prosite_path}{accession}.patmatmotif"
                ),
                stderr=subprocess.DEVNULL,
                shell=True,
            )
            # Opening the resultant file and using regex to add it to the dataframe.
            with open(f"{prosite_path}{accession}.patmatmotif", "r") as patmat:
                contents = patmat.read()
                list_of_matches = re.findall(r"(?<=Motif = )(.+)", contents)
                list_of_matches = ", ".join(list_of_matches)
                try:
                    seq_data.loc[
                        seq_data["Accession"] == accession, "Prosite_matches"
                    ] = list_of_matches
                except:
                    pass

    else:
        for selection in user_selection:
            for accession, sequence in zip(
                seq_data.loc[seq_data["Group_ID"] == selection]["Accession"],
                seq_data.loc[seq_data["Group_ID"] == selection]["Sequence"],
            ):
                fasta_formatted = ">" + accession + "\n" + sequence
                with open(f"{prosite_path}{accession}.fa", "w") as fasta_format_file:
                    fasta_format_file.write(fasta_formatted)
                subprocess.run(
                    (
                        f"patmatmotifs -sequence {prosite_path}{accession}.fa "
                        f"-outfile {prosite_path}{accession}.patmatmotif"
                    ),
                    stderr=subprocess.DEVNULL,
                    shell=True,
                )
                with open(f"{prosite_path}{accession}.patmatmotif", "r") as patmat:
                    contents = patmat.read()
                    list_of_matches = re.findall(r"(?<=Motif = )(.+)", contents)
                    list_of_matches = " ".join(list_of_matches)
                    try:
                        seq_data.loc[
                            seq_data["Accession"] == accession, "Prosite_matches"
                        ] = list_of_matches
                    except:
                        pass


#### INITIALIZING AN ARGUMENT PARSER AND POPULATING IT ####

# Initializing a parser object, to facilitate command line arguments to be passed to the program
args_parser = argparse.ArgumentParser(
    prog="B214618's Protein conservation tool",
    description=(
        "This program takes an input of a protein and a taxonomic group, and creates a"
        " number of groups using clustal omega's alignments. The user will then choose"
        " the group(s) they wish to query against the prosite database for conserved"
        " domain hits."
    ),
    epilog="",
)

# Adding arguments...
args_parser.add_argument(
    "--protein",
    dest="protein",
    required=True,
    help='use --protein to define protein query (e.g. "Pyruvate dehydrogenase")',
)
args_parser.add_argument(
    "--general-protein-search",
    dest="general",
    action="store_true",
    help=(
        "This will search the whole NCBI record for the protein name of interest,"
        " instead of just the protein_name section"
    ),
)
args_parser.add_argument(
    "--group",
    dest="grouping",
    required=True,
    help='use --group to define group query (e.g. "Ascomycota" or "txid4890")',
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
    default=20,
    type=int,
    help=(
        "--threads takes an integer, and is used to determine how many threads to use"
        " in all processes that allow thread number choices (currently just clustalo)."
        " Default is 20."
    ),
)
args_parser.add_argument(
    "--force",
    dest="force",
    action="store_true",
    help="Use --force to overwrite files. Default is false.",
)
args_parser.add_argument(
    "--no-grouping",
    dest="nogrouping",
    action="store_true",
    help=(
        "if --no-grouping is passed as an argument, only the first sequence alignment "
        "will be performed."
    ),
)

# assigning parsed args to variable
args = args_parser.parse_args()

#### CHECKING THAT TOOLS ARE INSTALLED ####

# maybe a bit of a roundabout way of doing it...
# I thought about running the user through an installation, but I could not reliably
# get the damn program to install. Sometimes it would just say "Unable to download EDirect archive",
# and you just had to try again later, and other times it would "complete" the process, but nothing
# would have been installed...! Additionally, I think it's good practice not to have a program
# install things by itself, in case perhaps the EDirect URL is changed, or the user's system
# is not the same as your own. Perhaps they would like to install it in an application folder,
# not their homespace or the current directory.

try:
    subprocess.Popen("efetch", stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
except FileNotFoundError:
    print(
        "You must install EDirect before continuing - see"
        " https://www.ncbi.nlm.nih.gov/books/NBK179288/"
    )

try:
    subprocess.Popen("clustalo", stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
except FileNotFoundError:
    print(
        "You must install Clustal Omega before continuing - see"
        " http://www.clustal.org/omega/"
    )

try:
    subprocess.Popen("needle", stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
except FileNotFoundError:
    print(
        "You must install EMBOSS Tools before continuing - see"
        " https://emboss.sourceforge.net/download/"
    )

#### MAIN CODE ####
# Checking directories...

checkDirs(["Fasta_files", "MSA_files", "Summary_files", "Fasta_files/prosite_files"])
fasta_path = "Fasta_files/"
msa_path = "MSA_files/"
summary_path = "Summary_files/"
prosite_path = "Fasta_files/prosite_files/"

# Fetching sequence query info
if args.general:
    search_query = (
        subprocess.check_output(
            (
                f'esearch -db "protein" -query "{args.protein} AND'
                f' {args.grouping}[Organism] NOT partial[Properties]" | xtract -pattern'
                " ENTREZ_DIRECT -element Count"
            ),
            shell=True,
        )
        .decode("utf-8")
        .strip()
    )
else:
    search_query = (
        subprocess.check_output(
            (
                f'esearch -db "protein" -query "{args.protein}[Protein Name] AND'
                f' {args.grouping}[Organism] NOT partial[Properties]" | xtract -pattern'
                " ENTREZ_DIRECT -element Count"
            ),
            shell=True,
        )
        .decode("utf-8")
        .strip()
    )

if search_query == None:
    print(
        "Something has gone wrong retrieving your search from NCBI. You might have to"
        " try again later, or check your network."
    )

if search_query == 0:
    print(
        "No sequeces found for these search terms. Please try again, or broaden your"
        " search parameters"
    )
    sys.exit()

# A quick few lines to avoid accidental processing of huge volumes of data...
print(
    "\nThe number of sequences in your search are: "
    + str(search_query)
    + "."
    + "\nWould you like to continue...? Please note sequence entries >1000 may take a"
    " while to process."
)

continue_ornot = "None"
while continue_ornot not in {"y", "n"}:
    continue_ornot = input("Please enter y/n:")
    if continue_ornot == "y":
        pass
    elif continue_ornot == "n":
        sys.exit()
    else:
        print("choose y/n")


# Fetching sequence data using gpc format. Fasta was used previously, but unfortunately
# I had a couple nasty accidents with regex gone wrong, so defaulted to xml style files
# using xtract for robustness. The annoying consequence of this is much higher waiting
# times for sequences to download...

print("Fetching sequences...")
if args.general:
    sequence_data = subprocess.check_output(
        (
            'esearch -db "protein"'
            f' -query "{args.protein} AND'
            f' {args.grouping}[Organism] NOT partial[Properties]"'
            "  | efetch -format gpc "
            "  | xtract -pattern INSDSeq "
            "  -element "
            "  INSDSeq_accession-version "
            "  INSDSeq_organism "
            "  INSDSeq_sequence "
        ),
        shell=True,
    ).decode("utf-8")
else:
    sequence_data = subprocess.check_output(
        (
            'esearch -db "protein"'
            f' -query "{args.protein}[Protein Name] AND'
            f' {args.grouping}[Organism] NOT partial[Properties]" '
            "  | efetch -format gpc "
            "  | xtract -pattern INSDSeq "
            "  -element "
            "  INSDSeq_accession-version "
            "  INSDSeq_organism "
            "  INSDSeq_sequence"
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
with open(f"{fasta_path}fasta_formatted.fa", "w") as fasta_formatted_file:
    fasta_formatted_file.write(fasta_string)

# Running the primary MSA
print("Running primary MSA...")

# There is probably a much much more elegant way of handling this, but for now...
if args.nogrouping:
    if args.cluster != None:
        subprocess.run(
            (
                "clustalo "
                "--force "
                "--auto  "
                f"--cluster-size={args.cluster} "
                f"--threads={args.threads} "
                "--outfmt=msf "
                f"-i {fasta_path}fasta_formatted.fa "
                f"-o {msa_path}primary_align.msf"
            ),
            shell=True,
        )
    else:
        subprocess.run(
            (
                "clustalo "
                "--force "
                "--auto  "
                f"--threads={args.threads} "
                "--outfmt=msf "
                f"-i {fasta_path}fasta_formatted.fa "
                f"-o {msa_path}primary_align.msf"
            ),
            shell=True,
        )

else:
    if args.cluster != None:
        subprocess.run(
            (
                "clustalo "
                "--force "
                "--auto  "
                f" --cluster-size={args.cluster} "
                f" --threads={args.threads} "
                f" --clustering-out={summary_path}clusterfile.txt "
                " --outfmt=msf "
                f"-i {fasta_path}fasta_formatted.fa "
                f"-o {msa_path}primary_align.msf"
            ),
            shell=True,
        )
        clusterfile = open(f"{summary_path}clusterfile.txt", "r").read().split("\n")
        clusterfile = list(filter(None, clusterfile))
    else:
        subprocess.run(
            (
                "clustalo "
                "--force "
                "--auto  "
                f" --threads={args.threads} "
                f"--clustering-out={summary_path}clusterfile.txt "
                "--outfmt=msf "
                f"-i {fasta_path}fasta_formatted.fa "
                f"-o {msa_path}primary_align.msf"
            ),
            shell=True,
        )
        clusterfile = open(f"{summary_path}clusterfile.txt", "r").read().split("\n")
        clusterfile = list(filter(None, clusterfile))

# Input > Process > Output... Repeat!
# No grouping option chosen then...:
if args.nogrouping:
    # Basically doing the same as in groupwiseMSA, but for only one group
    cons_output = subprocess.check_output(
        # /dev/stdout used to route output to variable.
        f"cons -sequence {msa_path}primary_align.msf -outseq /dev/stdout",
        stderr=subprocess.DEVNULL,
        shell=True,
    )
    cons_summary = groupDisplay(seq_data, cons_output)
    print(cons_summary)
    prositeGroupSearch()
    df_tosave = seq_data[["Accession", "Full_Name", "Prosite_matches"]]
    df_tosave.dropna()
    print(df_tosave)
    df_tosave.to_csv(f"{summary_path}prosite_summary.csv", encoding="utf-8")
# Default option...:
else:
    cluster_dict = clusterIndexer(clusterfile)
    group_filenames = groupFasta(cluster_dict)
    cons_output = groupwiseMSA(group_filenames)

    # Not entirely sure I actually need to run this command, will test at a later date. Program works fine with it though.
    # The idea was to remove single groups which would then be represented as na in the dataframe, but I might have
    # Made that all work within the groupFasta function already.
    seq_data = seq_data.dropna()

    # Continuing with the process...
    cons_summary, group_options = groupDisplay(seq_data, cons_output)
    print(cons_summary)
    print("THIS SUMMARY HAS BEEN AUTOMATICALLY SAVED IN THE Summary_files FOLDER!\n")
    print("Please choose one or more groups to query for conserved domains...:\n")
    user_selection = groupChoose(group_options)
    prositeGroupSearch(user_selection)
    # Final output! This dataframe output version will hopefully be improved upon, and be a PDF report instead. Need to work out how to do that.
    df_temp = seq_data[seq_data["Group_ID"].isin(user_selection)]
    df_tosave = seq_data[["Accession", "Full_Name", "Prosite_matches"]]
    df_tosave = df_tosave.dropna()
    # Printing to screen
    print(df_tosave)
    df_tosave.to_csv(f"{summary_path}prosite_summary.csv", encoding="utf-8")

# This approach was abandoned due to amount of time required to learn the plotly PdfPages documentation! :(
#    # Generating a plot of the table and saving it to PDF
#    fig, ax = plt.subplots()
#    fig.patch.set_visible(False)
#    ax.axis('off')
#    ax.axis('tight')
#    ax.table(cellText=df_tosave.values, colLabels=df_tosave.columns, loc='center')
#    fig.tight_layout()
#
#    pdf_output = PdfPages(f"{summary_path}prosite_summary.pdf")
#    pdf_output.savefig(fig, bbox_inches='tight')
#    pdf_output.close()
