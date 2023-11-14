#!~/.conda/envs/NADcapturedSeq/bin/python

import pandas as pd
import sys

print("Extracting the alignment rates from {}".format(sys.argv[1]))

names = []
reads = []
alignment_rate = []
type_ = []
spike_RNA = []

with open(sys.argv[1], "r") as f:
    # split the file into chunks, each chunk corresponds to one alignment
    chunk = ""
    for line in f.readlines():      
        chunk = chunk + line

        # once the end of one chunk is reached, extract the data
        if line.endswith("rate\n"):

            # split in lines
            splitted_chunk = chunk.split("\n")

            # extract name 
            name = splitted_chunk[0].split(" ")[3]

            # extract spike RNA name
            spike_RNA.append(splitted_chunk[0].split(" ")[-2])

            # extract reads
            reads.append(splitted_chunk[1].split(" ")[0])

            # swith behaviour depending on name and add type info
            splitted_name = name.split(".")

            # add name info 
            names.append(splitted_name[0].strip("merged_trimmed_"))

            if splitted_name[2] == "extendedFrags":
                type_.append("merged")
                alignment_rate.append(splitted_chunk[6].split(" ")[0])


            elif splitted_name[2]== "notCombined_1":
                type_.append("paired")
                alignment_rate.append(splitted_chunk[15].split(" ")[0])

            # reset the chunk to save storage
            chunk = ""

# create a dataframe for convenient conversion to csv
df = pd.DataFrame({
    "sample":names, 
    "reads_amount": reads, 
    "alignment_rate": alignment_rate, 
    "type": type_, 
    "spike_RNA": spike_RNA
    })
df.to_csv("hisat2_spike_alignments.csv")

print("Dataframe successfully written to hisat2_spike_alignments.csv")



