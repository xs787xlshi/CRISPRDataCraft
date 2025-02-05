#!/usr/bin/env python


#########################################
# IGDB_take_desired_mutants.py
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
import time
import sys

# Alleles_frequency_table

def reverse_complement(sequence):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return ''.join(complement[nucleotide] for nucleotide in sequence[::-1])

def read_desired_seq(args):
    infile_name = args.desired_seq_file
    desired_seq_dic = {}
    try:
        with open(infile_name, 'r') as file:
            for line in file:
                columns = line.strip().split()
                if len(columns) == 2:
                    if columns[0].upper() not in desired_seq_dic:
                        desired_seq_dic[columns[0].upper()] = []
                    desired_seq_dic[columns[0].upper()].append(columns[1].upper())
                else:
                    print(f"Ignoring invalid line: {line.strip()}")
    except FileNotFoundError:
        print(f"file {file_path} not found.")
    except Exception as e:
        print(f"Error occurred while reading the file {e}")
    return desired_seq_dic


def read_alleles_frequency_table(args, desired_seq_dic):
    path_aft =  args.path_allele_frequency_table+"/IGDB_on_*/Alleles_frequency_table.zip"
    # path_aft =  args.path_allele_frequency_table+"/IGDB_on_TATPSPSS1_9_7A/Alleles_frequency_table.zip"
    output_file = args.output
    outfile = open(output_file, 'w')
    outfile.write("Aligned_Sequence,Reference_Sequence,lib_name,Reference_Name,Read_Status,n_deleted,n_inserted,n_mutated,#Reads,%Reads,Desired\n")
    desired_dic = {}
    file_list = glob.glob(path_aft)
    if not file_list:
        print("Error: No files found in the allele frequency table path.")
        sys.exit(1) 
    for filename in file_list:
        word = filename.split('/')[-2].split("_")
        lib_name = word[2]+'_'+word[3]+'_'+word[4]
        barcode_name = word[2]+'_'+word[3]
        if barcode_name.upper() not in desired_seq_dic:
            print(f"There is no desired sequence for this barcode {barcode_name}")
            continue
        print(lib_name)
        # desired_seq = 'AACAACGTAAGCGCTTACGCAC'
        # rc_desired_seq = reverse_complement(desired_seq)
        desired_list = desired_seq_dic[barcode_name]
        rc_desired_list = []
        for item in desired_list:
            rc_desired_list.append(reverse_complement(item))
        with zipfile.ZipFile(filename, 'r') as zip_ref:
            # loop through each file in the zip archive
            for file_info in zip_ref.infolist():
                # read the file content as a single string
                file_content = zip_ref.read(file_info.filename)
                # if the file content is in bytes, decode it to a string
                file_content = file_content.decode('utf-8')
                # split the content into lines and read line by line
                line_count = 0
                desired_seq_perc_check = 0
                out_list = []
                hit_list = []
                for line in file_content.splitlines():
                    word = line.split(",")
                    if word[0] == "Aligned_Sequence":
                        continue
                    line_count += 1
                    if line_count > args.top_desired_seq: # only take top args.top alignments
                        break
                    read_seq = word[0].replace('-', '')  # get the raw mutant sequence removing '-'
                    marker = 0
                    for idx in range(len(desired_list)): # scan all desired sequence list
                        if desired_list[idx] in read_seq or rc_desired_list[idx] in read_seq:
                            desired_seq_perc_check += float(word[-1])
                            hit_list.append("Y")
                            marker = 1
                    if marker == 0:
                        hit_list.append("N")
                    out_list.append(line)
            if desired_seq_perc_check >= args.desired_seq_perc_sum: # summation of all desired sequence is higher than 20%
                if lib_name not in desired_dic:
                    desired_dic[lib_name] = []
                out_str  = ''
                for idx1 in range(0,len(out_list)):
                    word = out_list[idx1].split(",")
                    out_str += word[0]+','+word[1]+','+lib_name
                    for idx2 in range(2,len(word)):
                        out_str += ','+word[idx2]
                    out_str += ','+hit_list[idx1]+"\n"
                out_str += "\n\n"
                desired_dic[lib_name].append(out_str)
    for lib_name in sorted(desired_dic.keys()):
        for item in desired_dic[lib_name]:
            outfile.write(item)
    outfile.close()

if __name__ == "__main__":
    # start time
    start_time = time.perf_counter()
    # initialize parser
    parser = argparse.ArgumentParser(description="Extract segments containing the desired edited site.",formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-aft', '--path_allele_frequency_table', type=str,  help='Path of the allele frequency table', default='', required=True)
    parser.add_argument('-o', '--output', type=str,  help='Output file directory', required=True)

    parser.add_argument('-dsf','--desired_seq_file', type=str,  help='A file containing two columns: barcode and desired_seq', default='', required=True)
    parser.add_argument('-top', '--top_desired_seq', type=int,  help='Retrieve the top n sequences with the highest prevalence',default = 10)

    parser.add_argument('-pct', '--desired_seq_perc_sum', type=int,  help='The percentage of the desired sequence must not be lower than the given value', default = 20)

    # read arguments from command line
    args = parser.parse_args()
    config = vars(args)
    print(config)

    # check the parameter input information 
    if not os.path.exists(args.path_allele_frequency_table):
        print(f"input the path of the allele frequency table.")
        sys.exit(1)

    desired_seq_dic = read_desired_seq(args)
    '''
    for barcode_name, list in desired_seq_dic.items():
        for item in list:
            print(barcode_name, item)
    exit(0)
    '''
    read_alleles_frequency_table(args, desired_seq_dic)

    # get end time
    end_time = time.perf_counter()
    # calculating running time
    runTime = end_time - start_time
    # output running time
    print("Elapsed Time is: ", runTime, "Second\n")
