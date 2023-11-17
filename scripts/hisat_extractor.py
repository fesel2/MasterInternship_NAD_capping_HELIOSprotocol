#!~/.conda/envs/NADcapturedSeq/bin/python

import pandas as pd
import sys

print("Extracting the alignment rates from {}".format(sys.argv[1]))

names = []
reads = []
aligned_reads = []
overall_alignment = []
type_ = []
spike_RNA = []

with open(sys.argv[1], "r") as f:
    # split the file into chunks, each chunk corresponds to one alignment
    chunk = ""
    for line in f:

        # throw error message when an alignment failed
        if line.startswith("Error"):
            raise ValueError("""\nHisat output includes errors, Tty to fix them first.
                             \nContent of the corrupted line: {}""".format(line))
        else:
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
            names.append(splitted_name[0].replace("merged_trimmed_", ""))

            if splitted_name[2] == "extendedFrags":
                type_.append("merged")
                aligned_reads.append(splitted_chunk[4].split("(")[0].strip(" "))
                overall_alignment.append(splitted_chunk[6].split(" ")[0])

            elif splitted_name[2]== "notCombined_1":
                # extracts the counts for the three conditions: concordantly,
                # discordantly and neither nor (but divided by two, ~hisat2 discussion)
                type_.append("paired")
                aligned_reads_1 = int(splitted_chunk[4].split("(")[0].strip(" "))
                aligned_reads_2 = int(splitted_chunk[8].split("(")[0].strip(" "))
                aligned_reads_3 = int(splitted_chunk[13].split("(")[0].strip(" "))
                aligned_reads.append(aligned_reads_1 + aligned_reads_2 + aligned_reads_3/2)
                overall_alignment.append(splitted_chunk[15].split(" ")[0])

            # reset the chunk to save storage
            chunk = ""
        

# create a dataframe for convenient conversion to csv
df = pd.DataFrame({
    "sample":names, 
    "reads_amount": reads, 
    "aligned_reads": aligned_reads,
    "alignment_rate": overall_alignment,
    "type": type_, 
    "spike_RNA": spike_RNA
    })

df_merged = df[df["type"] == "merged"]
df_merged.to_csv("hisat2_spike_alignments_merged.csv", index=False)

df_merged = df[df["type"] == "paired"]
df_merged.to_csv("hisat2_spike_alignments_paired.csv", index=False)

print("Dataframes successfully written to hisat2_spike_alignments_merged\
       and hisat2_spike_alignments_paired.csv")



