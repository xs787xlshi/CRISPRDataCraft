# Import fasta run minced
#!/usr/bin/env python

#########################################
# IGDB_analyze_codon_shift.py
# Author: Xiaoli Shi
# Date: 2024.11.15
# (c) 2024 IGDB. All Rights Reserved.
#########################################



import os
import argparse
import gzip
import glob
import subprocess as sb
from pathlib import Path
import zipfile
from collections import defaultdict
import time
import sys


# both cis and trans strand
# stop_codons = ['TAA', 'TAG', 'TGA','TTA','CTA','TCA']
# only cis strand
stop_codons = ['TAA', 'TAG', 'TGA']

def check_stop_codon(dna_segment):
    for i in range(0, len(dna_segment), 3):
        codon = dna_segment[i:i+3]
        if codon in stop_codons:
            return True
    return False


def reverse_complement(sequence):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N','-':'-'}
    return ''.join(complement[nucleotide] for nucleotide in sequence[::-1])

# alleles_frequency_table
def codon_shift_parse(path_align, outfile1, outfile2):
    path_str = path_align+"/IGDB_on_*/Alleles_frequency_table.zip"
    output_file = outfile1
    outfile = open(output_file, 'w')

    # read_Status: WT, SNP, SNP_STOP,INDEL_3N+1, INDEL_3N+2, INDEL_3N, INDEL_3N_STOP 
    outfile.write("Aligned_Sequence,Reference_Sequence,lib_name,Reference_Name,n_deleted,n_inserted,n_mutated,#Reads,%Reads,Read_Status\n")
    mut_type_list = ["WT","SNP","SNP_STOP","3N+1", "3N+2", "3N", "3N_STOP"]
    mut_info_dic = defaultdict(lambda: defaultdict())
    file_list = glob.glob(path_str)
    if not file_list:
        print("Error: No files found in the specified path.")
        sys.exit(1)
    for filename in file_list:
        # print(filename)
        word = filename.split('/')[-2].split("_")
        lib_name = word[2]+'_'+word[3]+'_'+word[4]
        for item in mut_type_list:
            mut_info_dic[lib_name][item] = 0
        print(lib_name)
        with zipfile.ZipFile(filename, 'r') as zip_ref:
            # loop through each file in the zip archive
            for file_info in zip_ref.infolist():
                # read the file content as a single string
                file_content = zip_ref.read(file_info.filename)
                # if the file content is in bytes, decode it to a string
                file_content = file_content.decode('utf-8')
                # split the content into lines and read line by line
                for line in file_content.splitlines():
                    (Aligned_Sequence, Reference_Sequence, Reference_Name, Read_Status, n_deleted, n_inserted, n_mutated, Reads_numb, Reads_perct) = line.split(",")
                    if Aligned_Sequence == "Aligned_Sequence":
                        continue
                    if Read_Status == "UNMODIFIED":
                        read_status = "WT"
                    elif int(n_deleted) == 0 and int(n_inserted) == 0:##Without Indel
                        # SNP
                        mut_seq = Aligned_Sequence.upper()
                        # rev_mut_seq = reverse_complement(mut_seq)
                        if check_stop_codon(mut_seq):
                            read_status = "SNP_STOP"
                        else:
                            read_status = "SNP"
                    else: # with indel
                        # INDEL
                        dif = int(n_inserted)-int(n_deleted)
                        if dif%3 == 1:
                            read_status = "3N+1"
                        elif dif%3 == 2:
                            read_status = "3N+2"
                        else:
                            mut_seq = Aligned_Sequence.replace("-", "").upper()
                            if check_stop_codon(mut_seq):
                                read_status = "3N_STOP"
                            else:
                                read_status = "3N"
                    outfile.write(Aligned_Sequence+','+ Reference_Sequence+','+lib_name+','+Reference_Name+','+n_deleted+','+n_inserted+','+n_mutated+','+Reads_numb+','+Reads_perct+","+read_status+'\n')
                    mut_info_dic[lib_name][read_status] += int(Reads_numb)
    outfile.close()
    output_file = outfile2
    outfile = open(output_file, 'w')
    outfile.write("lib_name,WT,SNP,SNP_STOP,3N+1,3N+2,3N,3N_STOP,Total\n")
    sum = 0
    for lib_name in sorted(mut_info_dic.keys()):
        outfile.write(lib_name)
        for item in mut_type_list:
            outfile.write(','+str(mut_info_dic[lib_name][item]))
            sum += mut_info_dic[lib_name][item]
        outfile.write(','+str(sum)+'\n')
        sum = 0
    outfile.close()

if __name__ == "__main__":
    # start time
    start_time = time.perf_counter()
    # initialize parser
    parser = argparse.ArgumentParser(description="Tools to summarize codon shift based on alignments of target and mutants",formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-pa', '--path_alignment', type=str,  help='Path of the alignment files', default='', required=True)
    parser.add_argument('-o1', '--output1', type=str,  help='Output file of codon shift status of all reads\nThere are seven shift types: WT, SNP, SNP_stop, IND_3N+1, IND_3N+2, IND_3N, 3N_stop', required=True)
    parser.add_argument('-o2', '--output2', type=str,  help='Output file of integrated codon shift status of reads with the same barcode', required=True)

    # read arguments from command line
    args = parser.parse_args()
    config = vars(args)
    print(config)

    path_align = config['path_alignment']
    outfile1 = config['output1']
    outfile2 = config['output2']

    print(path_align, outfile1, outfile2)

    codon_shift_parse(path_align, outfile1, outfile2)

    # get end time
    end_time = time.perf_counter()
    # calculating running time
    runTime = end_time - start_time
    # output running time
    print("Elapsed Time is: ", runTime, "Second\n")

def deletion_parse(mut_sequence):
    dna_len = len(mut_sequence)
    for idx in range(0, dna_len,3):
        spos = idx
        epos = min(spos+3,dna_len)
        mut_codon = mut_sequence[spos:epos]
        del_count = sum(c== '-' for c in mut_codon)
        if mut_codon in stop_codons:
            return "INDEL_STOP"
        elif del_count == 1:
            return "INDEL_3N+2"
        elif del_count == 2:
            return "INDEL_3N+1"
    return "INDEL_3N"

def insertion_parse(mut_sequence, ref_sequence):
    dna_len = len(ref_sequence)
    for idx in range(0, dna_len,3):
        spos = idx
        epos = min(spos+3,dna_len)
        mut_codon = mut_sequence[spos:epos]
        if mut_codon in stop_codons:
            return "INDEL_STOP"
        ref_codon = ref_sequence[spos:epos]
        ins_count = sum(c=='-' for c in ref_codon)
        if ins_count == 1:
            return "INDEL_3N+1"
        elif ins_count == 2:
            return "INDEL_3N+2"
    return "INDEL_3N"
