#!/bin/python3

from subprocess import check_output

test_output = check_output('esearch -db protein -query "pyruvate dehydrogenase" | efilter -organism "Durotheca rogersii" | efetch -format fasta', shell = True)
print(test_output.decode('utf-8').strip('\n'))
