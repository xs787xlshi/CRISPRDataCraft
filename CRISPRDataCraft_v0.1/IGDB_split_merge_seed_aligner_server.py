# -*- coding: utf-8 -*-

'''
IGDB_Split_Merge_Seed_Aligner_server.py
Author: Xiaoli Shi
Date: 2025.01.14
Pipeline for the analysis of genome editing outcomes from deep sequencing data
(c) 2024 IGDB. All Rights Reserved.
'''

import difflib
import os
import sys
from copy import deepcopy
from datetime import datetime
import subprocess as sb
import glob
import gzip
import unicodedata
import string
import re
import zipfile
import traceback
import logging
import warnings
import argparse
import time

from collections import OrderedDict
import numpy as np
import pandas as pd
import logging

import multiprocessing as mp
import signal


from functools import partial

_jp=lambda OUTPUT_DIRECTORY, filename: os.path.join(OUTPUT_DIRECTORY, filename)

def filter_fastqs(merged_fastq, output_filename, min_av_read_qual=None, min_bp_qual_or_N=None):
    if not os.path.exists(merged_fastq):
        raise Exception("merged fastq file '"+ merged_fastq +"' does not exist.")

    # Create filehandles
    if merged_fastq.endswith('.gz'):
        f1_in = io.BufferedReader(gzip.open(merged_fastq, 'rb'))
        f1_out_filename=merged_fastq.replace('.fastq', '').replace('.gz', '')+'_filtered.fastq.gz'
        if output_filename:
            f1_out_filename = output_filename
    else:
        f1_in = open(merged_fastq, 'rb')
        f1_out_filename=merged_fastq.replace('.fastq', '')+'_filtered.fastq'
        if output_filename:
            f1_out_filename = output_filename

    if f1_out_filename.endswith('.gz'):
        f1_out = gzip.open(f1_out_filename, 'wt')
    else:
        f1_out = open(f1_out_filename, 'w')

    idLine = f1_in.readline().rstrip().decode('utf-8')
    while idLine:
        seqLine = f1_in.readline().rstrip()
        plusLine = f1_in.readline().rstrip()
        qualLine = f1_in.readline().rstrip()
        npQualLine = np.frombuffer(qualLine, dtype=np.uint8)-33 # assume illumina 1.7
        mean = np.mean(npQualLine)
        if mean >= min_av_read_qual:
            '''
            npSeqLine = np.frombuffer(seqLine, 'c').copy()
            # replace nucleotide with low quality by "N" 
            npSeqLine[npQualLine < min_bp_qual_or_N] = 'N'
            f1_out.write("%s\n%s\n%s\n%s\n"%(idLine, npSeqLine.tobytes().decode('utf-8'), plusLine.decode('utf-8'), qualLine.decode('utf-8')))
            '''
            f1_out.write("%s\n%s\n%s\n%s\n"%(idLine, seqLine.decode('utf-8'), plusLine.decode('utf-8'), qualLine.decode('utf-8')))
        idLine = f1_in.readline().rstrip().decode('utf-8')



class SplitBarcodeCleanReads:
    def __init__(self, barcode, fastq1, fastq2, amplicon, out_type, output_dir, n_processes_for_pooled):
        self.barcode = barcode
        self.fastq1 = fastq1
        self.fastq2 = fastq2
        self.amplicon = amplicon
        self.out_type = out_type
        self.output_dir = output_dir
        self.n_processes_for_pooled = n_processes_for_pooled

        self.out_obj_dic = {}
        self.QUALCUTOFF = 0

        self.snp_site_dic = {}
        self.amplicon_site_coverage = {}

        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)

        self.of_stat = open(self.output_dir + "/IGDB_reads_stat.csv", "w")
        self.of_stat.write("#Sample\tTotal_reads\tAlign_reads\tAlign_ratio\tUnAlign_read\n")
        '''
        self.of_result = open(self.output_dir+"/IGDB_"+out_type+"_results.csv", "w")
        if out_type == "CBE":
            self.of_result.write("#Sample\tAlign_reads\tC>T\tC>T(%)\tonly_C\tonly_C(%)\tInsertion\tInsertion(%)\tDeletion\tDeletion(%)\tInDel(%)\tEditing_reads\tEditing(%)\n")
        elif out_type == "ABE":
            self.of_result.write("#Sample\tAlign_reads\tA>G\tA>G(%)\tonly_A\tonly_A(%)\tInsertion\tInsertion(%)\tDeletion\tDeletion(%)\tInDel(%)\tEditing_reads\tEditing(%)\n")
        if out_type == "STEME":
            self.of_result.write("#Sample\tAlign_reads\tC>T\tC>T(%)\tonly_C\tonly_C(%)\tA>G\tA>G(%)\tonly_A\tonly_A(%)\tBoth_C>T&A>G\tBoth_C>T&A>G(%)\tInsertion\tInsertion(%)\tDeletion\tDeletion(%)\tInDel(%)\tEditing_reads\tEditing(%)\n")
        '''
    def read_barcode(self):
        name_list = []
        bcl_list = []
        bcr_list = []
        with open(self.barcode, 'r') as file:
            first_line = file.readline()
            if re.search(r'\s+',first_line):
                delimiter = r'\s+'
            elif ',' in first_line:
                delimiter = ','
            else:
                raise ValueError("Unsupported delimiter in the barcode file")
            file.seek(0)
            for line in file:
                if len(line.strip()) == 0:
                    continue
                if delimiter == r'\s+':
                    (bc_name, bc_l, bc_r) = re.split(delimiter, line.strip())
                else:
                    (bc_name, bc_l, bc_r) = line.strip().split(delimiter)
                bl_tmp = bc_l.upper().replace(' ','')
                br_tmp = bc_r.upper().replace(' ','')
                name_list.append(bc_name)
                bcl_list.append(bl_tmp)
                bcr_list.append(br_tmp)
                if not os.path.exists(self.output_dir+'/'+bc_name):
                    os.makedirs(self.output_dir+'/'+bc_name)
                # if (not os.path.exists(self.output_dir+'/'+bc_name+'/R1.fastq')) or (os.path.exists(self.output_dir+'/'+bc_name+'/R1.fastq') and os.path.getsize(self.output_dir+'/'+bc_name+'/R1.fastq') == 0):

                self.out_obj_dic[bc_name] = []
                if os.path.exists(self.output_dir+'/R1_unBarcode.fastq') and os.path.getsize(self.output_dir+'/R1_unBarcode.fastq') > 0 and \
                   os.path.exists(self.output_dir+'/R1_unBarcode.fastq') and os.path.getsize(self.output_dir+'/R1_unBarcode.fastq') > 0:
                    pass
                else:
                    R1_fq = open(self.output_dir+'/'+bc_name+'/R1.fastq', "w")
                    R2_fq = open(self.output_dir+'/'+bc_name+'/R2.fastq', "w")
                    self.out_obj_dic[bc_name].append(R1_fq)
                    self.out_obj_dic[bc_name].append(R2_fq)
        '''
        for item in name_list:
            print(item)
        exit(0)
        '''
        return name_list, bcl_list, bcr_list

    def get_total_reads(self):
        with gzip.open(self.fastq1, 'rt') as f1:
            total_reads = sum(1 for _ in f1) // 4
        return total_reads


    def split_barcode_clean_read(self, name_list, bcl_list, bcr_list):
        print(f"Reading {fastq1} and {fastq2} files...")

        # p = sb.Popen(('z' if self.fastq1.endswith('.gz') else '' ) +"cat < %s | wc -l" % self.fastq1, shell=True, stdout=sb.PIPE)
        # n_reads = int(float(p.communicate()[0])/4.0)

        bc_stat = {}
        bc_stat["unbc"] = 0

        # check if unbarcode file is exist
        if os.path.exists(self.output_dir+'/R1_unBarcode.fastq') and os.path.getsize(self.output_dir+'/R1_unBarcode.fastq') > 0 and \
           os.path.exists(self.output_dir+'/R2_unBarcode.fastq') and os.path.getsize(self.output_dir+'/R2_unBarcode.fastq') > 0:
            print("Split barcode files exist and are not empty. If you want to redo split barcode  please delete R1_unBarcode.fastq and R2_unBarcode.fastq.")
            return

        with gzip.open(self.fastq1,'rt') as f1,\
             gzip.open(self.fastq2,'rt') as f2,\
             open(self.output_dir+'/R1_unBarcode.fastq', "w") as unbc_fq1, \
             open(self.output_dir+'/R2_unBarcode.fastq', "w") as unbc_fq2:
            n_reads = self.get_total_reads()
            for i in range(0, n_reads):
                flag = 0
                long1 = f1.readline() + f1.readline() + f1.readline() + f1.readline()
                long2 = f2.readline() + f2.readline() + f2.readline() + f2.readline()
                for k in range(0, len(name_list)):
                    bc_l = bcl_list[k]
                    bc_r = bcr_list[k]
                    bc_name = name_list[k]
                    if (bc_l in long1 and bc_r in long2) or (bc_l in long2 and bc_r in long1):
                        # print(i)
                        flag = 1
                        if bc_name not in bc_stat:
                            bc_stat[bc_name] = 0
                        bc_stat[bc_name] += 1
                        self.out_obj_dic[bc_name][0].write(long1)
                        self.out_obj_dic[bc_name][1].write(long2)
                        break
                if flag == 0:
                    bc_stat["unbc"] += 1
                    unbc_fq1.write(long1)
                    unbc_fq2.write(long2)
        f1.close()
        f2.close()
        unbc_fq1.close()
        unbc_fq2.close()
        for bc_name in self.out_obj_dic:
            self.out_obj_dic[bc_name][0].close()
            self.out_obj_dic[bc_name][1].close()

        bc_stat_of = open(self.output_dir+'/Barcode.stat.out', "w")
        out_sort = {}
        for k in range(0, len(name_list)):
            bc_name = name_list[k]
            bc_l = bcl_list[k]
            bc_r = bcr_list[k]
            if bc_name not in bc_stat:
                bc_stat[bc_name] = 0
                tmp_perc = 0.0
            else:
                tmp_perc = 100*bc_stat[bc_name]/n_reads
            tmp_str = bc_name+'\t'+bc_l+'\t'+bc_r+'\t'+str(bc_stat[bc_name])+'\t'+str("%.2f"%(tmp_perc))
            if tmp_perc not in out_sort:
                out_sort[tmp_perc] = []
            out_sort[tmp_perc].append(tmp_str)
        for k in reversed(sorted(out_sort)):
            for item in out_sort[k]:
                bc_stat_of.write(item+"%\n")
        tmp_perc = 100*bc_stat["unbc"]/n_reads
        bc_stat_of.write('unbarcode\tunbarcode\t'+str(bc_stat["unbc"])+'\t'+str("%.2f"%(tmp_perc))+"%\n")
        bc_stat_of.write('All READS\t'+str(n_reads)+'\n')
        bc_stat_of.close()

    def merge_paired_reads(self, name_list):
        def subprocess_popen(statement):
            p = sb.Popen(statement, shell=True, stdout=sb.PIPE)  # implement shell command and define the output format
            while p.poll() is None:  # check if the process is completed (Popen.poll() Check if the subprocess (command) is completed，return NONE if not finished yet，return state code if finished）
                if p.wait() != 0:  # check if the process is succeful Popen.wait() Waiting for the end of subprocess and return the stat code；if the process is not finished until the timeout designate seconds，then through out a 'TimeoutExpired' Error.)
                    print("The command implement is failed，Check the connection of device.")
                    return False
                else:
                    re = p.stdout.readlines()  # get the original results
                    result = []
                    for i in range(len(re)):  # change code for the original results，convert to utf8 code and remove \n
                        res = re[i].decode('utf-8').strip('\r\n')
                        result.append(res)
                    return result


        word = self.output_dir.split("/")
        for bc_name in name_list:
            prefix = word[-1]
            tmp_R = self.output_dir+"/"+bc_name

            left_reads = tmp_R+"/R1.fastq"
            right_reads = tmp_R+"/R2.fastq"
            if os.path.exists(left_reads) and os.path.getsize(left_reads) > 0 and\
               os.path.exists(right_reads) and os.path.getsize(right_reads) > 0:
                # statment = "flash2 -d "+self.output_dir+"/merge -M 150 -O -o "+self.output_dir+"_"+bc_name+" ./"+self.output_dir+"/"+bc_name+"/R1.fastq ./"+self.output_dir+"/"+bc_name+"/R2.fastq"
                merged_outfile = self.output_dir+'/'+prefix+"_"+bc_name+".extendedFrags.fastq"
                if os.path.exists(merged_outfile) and os.path.getsize(merged_outfile) > 0:
                    pass
                else:
                    statment = "flash2 -d "+self.output_dir+"/merge -M 150 -o "+prefix+"_"+bc_name+" "+tmp_R+"/R1.fastq "+tmp_R+"/R2.fastq"
                    print(statment)
                    subprocess_popen(statment)
            else:
                warnings.warn(f"The split barcode files '{left_reads}' and '{right_reads}' do not exist or are empty.", category=UserWarning)
                # raise Exception("Split barcode files "+ left_reads +" and "+right_reads +"' do not exist or are empty.")

    def read_amplicon(self):
        infile = open(self.amplicon, "r")
        amplicon_db_dic = {}
        lines_list = infile.readlines()

        i = 0
        while i < len(lines_list):
            line = lines_list[i].strip()
            if ">" not in line:
                i += 1
                continue
            if "_" not in line:
                raise Exception(f"Please append a suffix for '{line}' and separate it with a underline ('_').")
            amplicon_name = line.replace(">","").replace("\-","").replace(" ","").upper()
            tmp_list = []
            for j in range(1,4):
                i += 1
                # print("in:"+str(i))
                if i == len(lines_list):
                    break
                if ">" in lines_list[i]:
                    print(f"The {self.amplicon} file format is not desired\n")
                    exit(1)
                tmp_str = lines_list[i].strip("\n").replace(" ","").upper()
                tmp_list.append(tmp_str)
            # print("out:"+str(i))
            if len(tmp_list) == 3:
                amplicon_db_dic[amplicon_name] = tmp_list
            i += 1
        return amplicon_db_dic
 
    def run_crisprdatacraft_cmd(self, name_list, amplicon_db_dic):
        # run_Crisprdatacraft_cmds(self, amplicon_db_dic, name_list):
        crisprdatacraft_cmds = []
        # if the reads number for each amplicon < 100, then go to next amplicon
        # scan the barcode name list
        for bc_i in range(0, len(name_list)):
            bc_name = name_list[bc_i]
            # print(bc_name)

            word = self.output_dir.split("/")
            if len(word) > 0:
                tmp_out = word[-1]
            else:
                print(f"Error: Directory '{output_dir}' does not exist. Program will now terminate.")
                raise SystemExit
            merged_fastq = self.output_dir+"/merge/"+tmp_out+"_"+bc_name+".extendedFrags.fastq"
            if not os.path.exists(merged_fastq):
                # raise Exception("merged_fastq file '"+merged_fastq+"' does not exist.")
                warnings.warn("Merged fastq file '"+ merged_fastq +"' does not exist.", category = UserWarning)
            # print('\n Processing: %d, %s'%(bc_i,bc_name))

            filter_dir = self.output_dir+"/merge/"
            filter_output_filename =_jp(filter_dir, os.path.basename(merged_fastq.replace('.fastq', '')).replace('.gz', '')+'_filtered.fastq.gz')
            # temprarily masked
            filter_fastqs(merged_fastq, filter_output_filename, min_av_read_qual = 30, min_bp_qual_or_N= 20)

            ampid_list = []
            amp_left_seq_list = []
            amp_right_seq_list = []
            amp_target_seq_list = []
            for amplicon_name in amplicon_db_dic: # The amplicon seq id and barcode id share a common string
                tmp_idx = bc_name.rfind("_")
                pre_bc_name = bc_name[:tmp_idx]

                tmp_idx = amplicon_name.rfind("_")
                pre_ap_name = amplicon_name[:tmp_idx]
                if pre_ap_name.upper() == pre_bc_name.upper():
                    ampid_list.append(amplicon_name)
                    amp_left_seq_list.append(amplicon_db_dic[amplicon_name][0])
                    amp_target_seq_list.append(amplicon_db_dic[amplicon_name][1])
                    amp_right_seq_list.append(amplicon_db_dic[amplicon_name][2])

            if len(ampid_list) == 0:
                print("Could not find corresponding genes for the barcode "+bc_name+" split sequences\n")
                continue

            if len(amp_left_seq_list) != len(amp_target_seq_list):
                raise Exception("The format of gene information for barcode "+bc_name+" is not correct\n")


            print("\nCalling base for barcode "+bc_name+"\n")
            for amp_idx in range(len(ampid_list)):
                ampid = ampid_list[amp_idx]
                amp_left_seq = amp_left_seq_list[amp_idx]
                amp_target_seq = amp_target_seq_list[amp_idx]
                amp_right_seq = amp_right_seq_list[amp_idx]

                crisprdatacraft_cmd = '/public-supool/home/gaolab/anaconda3/bin/python /public-supool/home/gaolab/SXL/OOP_Hybrid_pipeline/CRISPRDataCraft.py -fq %s --amplicon_seq %s -o %s --amplicon_name %s --quantification_window_center %d --quantification_window_size %d --left_primer_seq %s --right_primer_seq %s --n_processes_for_pooled %d --min_average_read_quality %d --min_bp_quality_or_N %d' % (filter_output_filename, amp_target_seq, self.output_dir, ampid, -3, 20, amp_left_seq, amp_right_seq, self.n_processes_for_pooled, 30, 20)

                # crisprdatacraft_cmd = 'python3 CRISPRDataCraft.py -fq %s --amplicon_seq %s -o %s --amplicon_name %s --quantification_window_center %d  --quantification_window_size %d --left_primer_seq %s --right_primer_seq %s --n_processes_for_pooled %d --min_average_read_quality %d --min_bp_quality_or_N %d' % (filter_output_filename, amp_target_seq, self.output_dir, ampid, -3, 20, amp_left_seq, amp_right_seq, self.n_processes_for_pooled, 30, 20)

                # print(crisprdatacraft_cmd)
                # sb.run(crisprdatacraft_cmd, shell=True)
                # crisprdatacraft_cmds.append(crisprdatacraft_cmd)

                prefix_str = self.output_dir.split("/")[-1]+'_'+name_list[bc_i]+'_'+ampid.split("_")[1]
                content = "#BSUB -L /bin/bash\n#BSUB -J "+prefix_str+"\n#BSUB -n 1\n#BSUB -e %J.err \n#BSUB -o %J.out\n#BSUB -q standardA\n\n"+crisprdatacraft_cmd+"\n"
                pbs_file = self.output_dir+'/'+prefix_str+".basecall.pbs"
                with open(pbs_file, 'w') as file:
                    file.write(content)

                sb.call("bsub < "+pbs_file, shell=True)


if __name__ == '__main__':
    # start time
    start_time = time.perf_counter()
    # initialize parser
    parser = argparse.ArgumentParser(description="Tools for the analysis of genome editing outcomes from deep sequencing data",formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-r1', '--fastq_r1', type=str,  help='First fastq file', default='', required=True)
    parser.add_argument('-r2', '--fastq_r2', type=str,  help='Second fastq file for paired end reads', default='', required=True)
    #parser.add_argument('-c', '--change_number', type=int,  help='Expected Maximum Number of Mutated Bases', required=True)

    parser.add_argument('-b', '--barcode', type=str,  help='Barcode file', required=True)
    parser.add_argument('-a', '--amplicon_seq', type=str,  help='Amplicon Sequences file', required=True)
    parser.add_argument('-et', '--edited_type', type=str,  help='Type of edit (ABE, CBE, STEME)', default = "STEME", required=False)

    parser.add_argument('-o', '--output', type=str,  help='Output file directory', required=True)
    parser.add_argument('-p', '--n_processes_for_pooled', type = int, default=1,  help='Multiple processes number')

    # read arguments from command line
    args = parser.parse_args()
    config = vars(args)
    print(config)

    fastq1 = config['fastq_r1']
    fastq2 = config['fastq_r2']  

    # change_number = config['change_number']
    barcode = config['barcode']
    amplicon = config['amplicon_seq']
    out_type = config['edited_type']
    output_dir = config['output']
    n_processes_for_pooled = args.n_processes_for_pooled

    print(fastq1, fastq2, barcode, amplicon, output_dir, n_processes_for_pooled)

    # check the parameter input information 
    if not os.path.exists(barcode):
        print(f"Error: The file '{barcode}' does not exist.")
        sys.exit(1)

    if not os.path.exists(amplicon):
        print(f"Error: The file '{amplicon}' does not exist.")
        sys.exit(1)

    if not os.path.exists(fastq1):
        print(f"Error: The file '{fastq1}' does not exist.")
        sys.exit(1)

    if not os.path.exists(fastq2):
        print(f"Error: The file '{fastq2}' does not exist.")
        sys.exit(1)

    if fastq1 == fastq2:
        print(f"Error: The file fastq1: '{fastq1}' and fastq2 '{fastq2}' should not be the same.")
        sys.exit(1)


    name_list = []
    bcl_list = []
    bcr_list = []
    amplicon_db_dic = {}

    my_split_barcode = SplitBarcodeCleanReads(barcode, fastq1, fastq2, amplicon, out_type, output_dir, n_processes_for_pooled)

    # read barcode
    print("read barcode.....\n")
    name_list, bcl_list, bcr_list  = my_split_barcode.read_barcode()
    '''
    for idx in range(len(name_list)):
        print(name_list[idx], bcl_list[idx], bcr_list[idx])
    exit(0)
    '''
    # read amplicon
    print("read amplicon.....\n")
    amplicon_db_dic =  my_split_barcode.read_amplicon()
    '''
    for k in amplicon_db_dic.keys():
        print(k)
        for idx in range(0,3):
            print("***"+amplicon_db_dic[k][idx]+"***")
    exit(0)
    '''

    # split barcode
    print("split barcode.....\n")
    my_split_barcode.split_barcode_clean_read(name_list, bcl_list, bcr_list)

    # merge paired reads
    print("merge paired reads.....\n")
    my_split_barcode.merge_paired_reads(name_list)

    # running basecalling command
    my_split_barcode.run_crisprdatacraft_cmd(name_list, amplicon_db_dic)

    name_list = []
    bcl_list = []
    bcr_list = []
    amplicon_db_dic = {}

    # get end time
    end_time = time.perf_counter()
    # calculating running time
    runTime = end_time - start_time
    # output running time
    print("Elapsed Time is: ", runTime, "Second\n")
