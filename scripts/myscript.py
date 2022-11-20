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

# Argument to forcibly overwrite files
args_parser.add_argument('--force', dest='force', action='store_true', help = "Use --force to overwrite files. Default is false.")


# assigning parsed args to variable
args = args_parser.parse_args()

# Creating empty dataframe to hold sequence info
seq_data = pd.DataFrame(columns = ['Accession', 'Protein name', 'Genus', 'Species', 'Sequence'])

# Fetching sequences the easy (boring) way
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
print("Fetching sequences...")
fasta_seqs = subprocess.check_output(f"esearch -db \"{args.database}\" -query \"{args.protein}[Protein Name] AND {args.grouping}[Organism] NOT partial[Properties]\" | efetch -format fasta", shell = True).decode('utf-8')
print("Building dataframe...")
fasta_seqs_list = fasta_seqs.split('>')
#print("got sequences!")
# Filtering pesky empty list strings
fasta_seqs_list = list(filter(None, fasta_seqs_list))

list_of_rows=[]
# Slightly roundabout way of getting values to populate my dataframe...
for entry in fasta_seqs_list:
    accession_code = entry.split()[0]
    protein_name = re.search(r'(?<=' + accession_code + r')(.*)(?=\[)', entry).group(1).strip()
    species_block = re.search(r'(?<=\[)(.*)(?=\])', entry).group(1).split()
    genus = species_block[0]
    species = species_block[1]
    sequence = ''.join(entry.split('\n')[1:])
# Creating a list of 'dataframe rows'
    list_of_rows.append([accession_code, protein_name, genus, species, sequence])

#print(list_of_rows)
# Creating the dataframe, with relevant column names...
seq_data = pd.DataFrame(list_of_rows, columns = ['Accession', 'Protein name', 'Genus', 'Species', 'Sequence'])
seq_data['Binomial'] = seq_data['Genus'].map(str) + " " + seq_data['Species'].map(str)
#print(seq_data)
#print(seq_data['Genus'].unique())
#print(seq_data['Binomial'].unique())

# I would really rather use variables for the whole process instead of writing to a .fa file, but I don't know how. Maybe this really is the best way to go about it!
file_for_msa = open("seqs.fa", "w")
file_for_msa.write(fasta_seqs)
file_for_msa.close()

print("Running primary MSA...")
# modified from python docs examples. Might be a more robust way of doing it?
subprocess.run(f"clustalo --force --threads={args.threads} --clustering-out=clusterfile.txt --outfmt=msf -i seqs.fa -o align.msf", shell=True)
clusterfile = open("clusterfile.txt", "r").read().split('\n')
clusterfile =  list(filter(None, clusterfile))

#print(clusterfile)

#A couple of fun little for loops that iterate through the cluster file, picking out relevant values using regex and making them into a dictionary where the key is the group number and the value is a list of indexes for my seq_data dataframe
# Initialising some variables...
print("Building groupings...")
key_value_tuples = []
cluster_dict = {}
for file in clusterfile:
#    Extracting the bits I want using regex
    key = re.search(r'(?<=Cluster)(.*)(?=\:)', file).group(1).strip()
    value = re.search(r'(?<=index)(.*)(?=\()', file).group(1).strip()
    key_value_tuples.append((key, value))

# Creating a dictionary of index lists
for key, value in key_value_tuples:
    cluster_dict.setdefault(key, []).append(value)

#print(cluster_dict.keys())
#Processing the resultant groups to present to the user
for key in cluster_dict.keys():
    index_list = cluster_dict[key]
#    print("GROUP " + key + ":")
    fasta_string_group = ""
    group_df = seq_data.iloc[index_list]
#    Concatenating accession numbers and sequences into a fasta formatted string variable
    for accession, sequence in zip(group_df.Accession, group_df.Sequence):
#        fasta_string_group = fasta_string_group + ">" + accession + " GROUP " + key + "\n" + sequence + "\n"
        fasta_string_group = fasta_string_group + ">" + accession + "\n" + sequence + "\n"
    with open(f"group_{key}_msa.fa", "w") as groupalign:
        groupalign.write(fasta_string_group)
    print(f"Re-aligning group {key}...")
    subprocess.run(f"clustalo --auto --force --threads={args.threads} --outfmt=msf -i group_{key}_msa.fa -o group_{key}_msa.msf", shell=True)
# Calling the conserved sequence generator emboss program, shutting up its annoying "error" messages by sending them to /dev/null... Why does it send the message "Create a consensus sequence from a multiple alignment" to stderr???? Surely it should be to stdout?
    subprocess.call(f"cons -sequence group_{key}_msa.msf -outseq group_{key}_cons.txt", stderr=subprocess.DEVNULL, shell=True)
#    subprocess.run(f"cons {fasta_string_group}"

#def groupDisplay(consensus_list)

def groupChoose(group_options):

    possible_choices = []
    print('Please pick a group number, e.g. for \"Group 1\" enter \"1\". To finish, enter \"f\":')
    for option in group_options:
        possible_choices.append(option)
        print(" -- Group " + option)
    print("\n -- ALL GROUPS: a")
    print(" -- EXIT: f")
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

consensus_list = []
group_options = {}
for key in cluster_dict.keys():
    group_options[key] = f"Group {key}"
    consensus_list.append(open(f"group_{key}_cons.txt", "r").read())

count = 0
for cons in consensus_list:
    consensus_list[count] = "\n".join(cons.split("\n")[1:])
    print(consensus_list[count] + '\n')
    count += 1


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
