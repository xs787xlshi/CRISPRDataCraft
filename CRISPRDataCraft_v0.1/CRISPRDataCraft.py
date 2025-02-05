# -*- coding: utf-8 -*-

'''
Developed by modifying CRISPRessoCRISPResso2-2.2.14
Author: Xiaoli Shi
Date: 2024.06.17
Pipeline for the analysis of genome editing outcomes from deep sequencing data

1) Use the Monte Carlo method to perform sequence alignment
2) Quantify the indels and substitutions for each amplicon
3) Aggregate the total number of indels and substitutions
4) Generate and display statistical results summarizing the indels and substitutions
5) Visualize the data by creating plots to illustrate the distribution of indels and substitutions

(c) 2024 IGDB. All Rights Reserved.
'''


import sys
running_python3 = False
if sys.version_info > (3, 0):
    running_python3 = True

import argparse
from collections import defaultdict
from copy import deepcopy
import concurrent.futures
import multiprocessing

import errno
import gzip
import json
import zipfile
import os
import re
import subprocess as sb
import traceback
import time
import math
import io
import numpy as np

import pandas as pd
from Bio.Align import PairwiseAligner
from Bio.Align import substitution_matrices

# change to your server path
module_path = '/public-supool/home/gaolab/SXL/OOP_Hybrid_pipeline'
# local module path
# module_path = '/mnt/Data_diskE/Gene_Editing_Research_Server/OOP_Hybrid_pipeline_v0.5/server_version'
if module_path not in sys.path:
    sys.path.append(module_path)

import CRISPRessoShared
import CRISPRessoPlot
import IGDB_pairwise_aligner
import CRISPRessoCOREResources

from datetime import datetime
present = datetime.now()
import logging

logging.basicConfig(
                     format='%(levelname)-5s @ %(asctime)s:\n\t %(message)s \n',
                     datefmt='%a, %d %b %Y %H:%M:%S',
                     stream=sys.stderr,
                     filemode="w"
                     )

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

error   = logger.critical
warn    = logger.warning
debug   = logger.debug
info    = logger.info

_jp=lambda OUTPUT_DIRECTORY, filename: os.path.join(OUTPUT_DIRECTORY, filename)
get_name_from_fasta=lambda  x: os.path.basename(x).replace('.fastq', '').replace('.gz', '').replace('.fq', '')
get_name_from_bam=lambda  x: os.path.basename(x).replace('.bam', '')
_ROOT = os.path.abspath(os.path.dirname(__file__))

def get_n_reads_fastq(fastq_filename):
    p = sb.Popen(('z' if fastq_filename.endswith('.gz') else '' ) +"cat < \"%s\" | wc -l" % fastq_filename, shell=True, stdout=sb.PIPE)
    return int(float(p.communicate()[0])/4.0)

def check_file(filename):
    try:
        with open(filename): pass
    except IOError:
        files_in_curr_dir = os.listdir('.')
        if len(files_in_curr_dir) > 15:
            files_in_curr_dir = files_in_curr_dir[0:15]
            files_in_curr_dir.append("(Complete listing truncated)")
        dir_string = ""
        file_dir = os.path.dirname(filename)
        if file_dir == "":
            dir_string = ""
        elif os.path.isdir(file_dir):
            files_in_file_dir = os.listdir(file_dir)
            if len(files_in_file_dir) > 15:
                files_in_file_dir = files_in_file_dir[0:15]
                files_in_file_dir.append("(Complete listing truncated)")
            dir_string = "\nAvailable files in " + file_dir + ":\n\t" + "\n\t".join(files_in_file_dir)
        else:
            dir_string = "\nAdditionally, the folder '" + os.path.dirname(filename) + "' does not exist"

        raise BadParameterException("The specified file '"+filename + "' cannot be opened.\nAvailable files in current directory:\n\t" + "\n\t".join(files_in_curr_dir) + dir_string)


def get_allele_row(reference_name, variant_count, aln_ref_names, aln_ref_scores, variant_payload, write_detailed_allele_table):
    """
    gets a row for storing allele information in the allele table
    parameters:
    reference_name: string Reference name to write
    variant_count: number of times this allele appears
    aln_ref_names_str: '&'-joined string of references this allele aligned to
    aln_ref_scores_str: '&'-joined string of references this allele aligned to
    variant_payload: payload object (dict) with keys containing information about the allele
    write_detailed_allele_table: bool for whether to write detailed row
    returns:
        row to put into allele table
    """
    if write_detailed_allele_table:
        allele_row = {'#Reads':variant_count,
            'Aligned_Sequence': variant_payload['aln_seq'],
            'Reference_Sequence':variant_payload['aln_ref'],
            'n_inserted':variant_payload['insertion_n'],
            'n_deleted':variant_payload['deletion_n'],
            'n_mutated':variant_payload['substitution_n'],
            'Reference_Name':reference_name,
            'Read_Status':variant_payload['classification'],
            'Aligned_Reference_Names':aln_ref_names,
            'Aligned_Reference_Scores':aln_ref_scores,
            'ref_positions':variant_payload['ref_positions'],
            'all_insertion_positions': variant_payload['all_insertion_positions'], # arr with 1's where there are insertions (including those outside of include_idxs quantification window)
            'all_insertion_left_positions': variant_payload['all_insertion_left_positions'], # arr with 1's to the left of where the insertion occurs
            'insertion_positions': variant_payload['insertion_positions'], # arr with 1's where there are insertions (1bp before and 1bp after insertion) that overlap with include_idxs quantification window
            'insertion_coordinates': variant_payload['insertion_coordinates'], # one entry per insertion, tuple of (start,end)
            'insertion_sizes': variant_payload['insertion_sizes'],
            'all_deletion_positions': variant_payload['all_deletion_positions'], # arr with 1's where there are insertions
            'deletion_positions': variant_payload['deletion_positions'], # arr with 1's where there are insertions that overlap the include_idxs quantification window
            'deletion_coordinates': variant_payload['deletion_coordinates'], # one entry per deletion
            'deletion_sizes': variant_payload['deletion_sizes'], # correspond to entries in 'deletion_coordinates'
            'all_substitution_positions': variant_payload['all_substitution_positions'],
            'substitution_positions': variant_payload['substitution_positions'],
            'substitution_values': variant_payload['substitution_values']
    	}
    else:
        allele_row = {'#Reads':variant_count,
            'Aligned_Sequence': variant_payload['aln_seq'],
            'Reference_Sequence':variant_payload['aln_ref'],
            'n_inserted':variant_payload['insertion_n'],
            'n_deleted':variant_payload['deletion_n'],
            'n_mutated':variant_payload['substitution_n'],
            'Reference_Name':reference_name,
            'Read_Status':variant_payload['classification'],
            'Aligned_Reference_Names':aln_ref_names_str,
            'Aligned_Reference_Scores':aln_ref_scores_str,
            'ref_positions':variant_payload['ref_positions']
    	}
    return allele_row

def save_vector_to_file(vector, filename):
    # np.savetxt(_jp('%s.txt' %name), np.vstack([(np.arange(len(vector))+1),vector]).T, fmt=['%d','%.18e'],delimiter='\t', newline='\n', header='amplicon position\teffect',footer='', comments='# ')
    np.savetxt(filename, np.vstack([(np.arange(len(vector))+1), vector]).T, fmt=['%d', '%.18e'], delimiter='\t', newline='\n', header='amplicon position\teffect', footer='', comments='# ')

def save_count_vectors_to_file(vectors, vectorNames, refSeq, filename):
    outfile = open(filename, "w")
    outfile.write("Sequence\t"+"\t".join(list(refSeq))+"\n") # first row: reference sequence
    for vector, vectorName in zip(vectors, vectorNames):
        outfile.write(vectorName +"\t" + "\t".join([str(x) for x in vector]) + "\n") # next, vectors are printed
    outfile.close()

def get_plot_title_with_ref_name(plotTitle, ref_name):
    # if n_refs > 1:
    return (plotTitle + ": " + ref_name)
    # return plotTitle

def count_alternate_alleles(sub_base_vectors, ref_name, ref_sequence, ref_total_aln_reads):
    # create vectors with all allele frequencies -- not just the substitution (the reference allele will not be 0)
    alph = ['A', 'C', 'G', 'T', 'N']

    # count the total number of times each substitution occurs
    count_sub_base_vectors = {}
    alt_nuc_counts = {}
    for a in alph:
        alt_nuc_counts[a] = {}
        count_sub_base_vectors[a] = list(sub_base_vectors[ref_name+"_"+a])
        for b in alph:
            alt_nuc_counts[a][b] = 0

    for idx, c in enumerate(ref_sequence):
        tot_sub_at_idx = 0
        for a in alph:
            sub = sub_base_vectors[ref_name+"_" + a][idx]
            alt_nuc_counts[c][a] += sub
            tot_sub_at_idx += sub

    # df_subs = pd.DataFrame([count_sub_base_vectors["A"],count_sub_base_vectors["C"],count_sub_base_vectors["G"],count_sub_base_vectors["T"],count_sub_base_vectors["N"]])
    df_subs = pd.DataFrame([count_sub_base_vectors[a] for a in alph])
    df_subs.index = alph
    df_subs.columns = list(ref_sequence)
    return (df_subs, alt_nuc_counts)

def reverse_complement(sequence):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return ''.join(complement[nucleotide] for nucleotide in sequence[::-1])

class Edition_stat_each_amplicon(object):
    def __init__(self, args):
        # load and validate template file
        bc_seg = re.search(r'_([^_]+)\.extendedFrags', args.merged_fastq).group(1)
        # self.OUTPUT_DIRECTORY='%s/IGDB_on_%s_%s' % (args.output_dir, args.amplicon_name, bc_seg)

        tmp_idx = args.amplicon_name.rfind("_")
        pre_ap_name = args.amplicon_name[:tmp_idx]
        pos_ap_name = args.amplicon_name[tmp_idx+1:]
        self.OUTPUT_DIRECTORY='%s/IGDB_on_%s_%s_%s' % (args.output_dir, pre_ap_name, bc_seg, pos_ap_name)
        print("OUTPUT_DIRECTORY: ", self.OUTPUT_DIRECTORY)


        # keep track of all information for this run to be pickled and saved at the end of the run
        self.crispresso2_info = {'running_info': {}, 'results': {'alignment_stats': {}, 'general_plots': {}}} 
        self.crispresso2_info['running_info']['version'] = CRISPRessoShared.__version__
        self.crispresso2_info['running_info']['args'] = deepcopy(args)


        log_filename=_jp(self.OUTPUT_DIRECTORY, 'CRISPResso_RUNNING_LOG.txt')
        self.crispresso2_info['running_info']['log_filename'] = os.path.basename(log_filename)

        crispresso_cmd_to_write = ' '.join(sys.argv)
        self.crispresso2_info['running_info']['command_used'] = crispresso_cmd_to_write
        # print("crispresso_cmd_to_write ", crispresso_cmd_to_write)

        try:
            os.makedirs(self.OUTPUT_DIRECTORY)
            info('Creating Folder %s' % self.OUTPUT_DIRECTORY)
        except:
            warn('Folder %s already exists.' % self.OUTPUT_DIRECTORY)

        finally:
            logger.addHandler(logging.FileHandler(log_filename))

            with open(log_filename, 'w+') as outfile:
                outfile.write('CRISPResso version %s\n[Command used]:\n%s\n\n[Execution log]:\n' %(CRISPRessoShared.__version__, crispresso_cmd_to_write))

        self.files_to_remove = [] #these files will be deleted at the end of the run

        # SET REFERENCES TO COMPARE
        self.amplicon_seq = args.amplicon_seq
        self.amplicon_name = args.amplicon_name
        wrong_nt=CRISPRessoShared.find_wrong_nt(self.amplicon_seq)
        if wrong_nt:
            raise Exception('Reference amplicon sequence %s contains invalid characters: %s'%(self.amplicon_name, ' '.join(wrong_nt)))

        self.aln_seed_len = int(math.log2(len(args.amplicon_seq))+1)
        self.aln_seed_count = int(self.aln_seed_len/2+1)
        # build reference object
        self.refs = {}

        # merged reads
        self.N_READS_INPUT = 0
        if args.merged_fastq:
            self.N_READS_INPUT = get_n_reads_fastq(args.merged_fastq)

        if self.N_READS_INPUT == 0:
            raise Exception('The input contains 0 reads.')

        # count filtered reads
        self.N_READS_AFTER_PREPROCESSING = 0

        # INITIALIZE CACHE
        self.variantCache = {}
        self.aln_stats = {}

        self.N_TOTAL = 0
        self.N_DISCARDED = 0

        self.counts_total = 0
        self.class_counts = {}

        self.base_count_vectors = {}
        self.all_base_count_vectors = {}

        # indel mutation vectors
        this_len_amplicon = len(args.amplicon_seq)
        self.all_insertion_count_vectors            = np.zeros(this_len_amplicon) # all insertions (including quantification window bases)
        self.all_insertion_left_count_vectors       = np.zeros(this_len_amplicon) # all insertions (including quantification window bases)
        self.all_deletion_count_vectors             = np.zeros(this_len_amplicon)
        self.all_substitution_count_vectors         = np.zeros(this_len_amplicon)
        self.all_indelsub_count_vectors             = np.zeros(this_len_amplicon)


        self.insertion_count_vectors                = np.zeros(this_len_amplicon) # insertions that are in the quantification window
        self.deletion_count_vectors                 = np.zeros(this_len_amplicon)
        self.substitution_count_vectors             = np.zeros(this_len_amplicon)
        self.indelsub_count_vectors                 = np.zeros(this_len_amplicon)

        self.effective_len_dicts = {} # dict of effective lengths for all reads
        self.inserted_n_dicts = {} # dict of number of insertions for all reads: dict{len}->number of reads with that insertion length


        #length by position in amplicon
        self.insertion_length_vectors               = np.zeros(this_len_amplicon)
        self.deletion_length_vectors                = np.zeros(this_len_amplicon)

        self.alleles_list = [] # will be turned into df with rows with information for each variant (allele)

        # count times substitutions occur
        self.all_substitution_base_vectors = {}
        # these two are taken as a subset of all_substitution_base_vectors afterward to save time
        self.substitution_base_vectors = {} 
        # end of __init__

    def get_unaligned_sequences(self, fq_seq_dic):
        unaligned_fastq_seq_list = []
        unaligned_fq_seq_count_list = []
        for fastq_seq, count in fq_seq_dic.items():
            if fastq_seq not in self.variantCache or self.variantCache[fastq_seq]['best_match_score'] < 0.9 * self.refs['sequence_length']:
                unaligned_fastq_seq_list.append(fastq_seq)
                unaligned_fq_seq_count_list.append(count)
        return unaligned_fastq_seq_list, unaligned_fq_seq_count_list


    def process_fastq(self, fastq_filename, args):
        N_TOT_READS = 0
        N_CACHED_ALN = 0      # reads found in cache
        N_CACHED_NOTALN = 0   # reads found in 'not aligned' cache
        N_COMPUTED_ALN = 0    # not in cache, aligned to at least 1 sequence with min cutoff
        N_COMPUTED_NOTALN = 0 # not in cache, not aligned to any sequence with min cutoff

        if fastq_filename.endswith('.gz'):
            fastq_handle = gzip.open(fastq_filename, 'rt')
        else:
            fastq_handle = open(fastq_filename)

        # read fastq file into a temporary dictionary
        consensus_seq_dic = {}
        while(fastq_handle.readline()):
            # read through fastq in sets of 4
            fastq_seq = fastq_handle.readline().strip()
            fastq_plus = fastq_handle.readline()
            fastq_qual = fastq_handle.readline()
            consensus_seq_dic[fastq_seq] = consensus_seq_dic.get(fastq_seq, 0) + 1
        fastq_handle.close()

        # extract the edited target sequences
        fq_seq_dic = {}
        for fastq_seq, seq_count in consensus_seq_dic.items():
            # use L to split read
            seq_word = fastq_seq.split(self.refs['left_primer_seq'])
            if len(seq_word) != 2:   # if L could not split a read, pick its reversed strand
	            seq_rev = reverse_complement(fastq_seq)
	            seq_word = seq_rev.split(self.refs['left_primer_seq'])

            if len(seq_word) != 2:
                # print("stopped in the first layer")
                continue # if the reversed strand was not split yet, go to next read

            # use R split 3' end segment of the read
            target_word = seq_word[1].split(self.refs['right_primer_seq'])
            if len(target_word) != 2:
                # print("stopped in the second layer")
                continue
            edited_target = target_word[0]

            if len(edited_target) == 0:
                continue
            fq_seq_dic[edited_target] = fq_seq_dic.get(edited_target, 0) + seq_count
            N_TOT_READS+= seq_count
        consensus_seq_dic = {}
        '''
        for k, v in fq_seq_dic.items():
            print(k,v)
        exit(0)
        '''
        info("Processing reads; N_TOT_READS: %d ......"%(N_TOT_READS))

        align_retry = 10 # try 10 times to scan best alignment
        for i_retry in range(align_retry):
            unaligned_fastq_seq_list, unaligned_fq_seq_count_list = self.get_unaligned_sequences(fq_seq_dic)
            if not unaligned_fastq_seq_list:
                break
            n = len(unaligned_fastq_seq_list)
            bg = [args.amplicon_name]*n
            logging.info(f"Retry {i_retry}, unaligned sequences: {len(unaligned_fastq_seq_list)}")

            # print("i_retry ", str(i_retry))
            # print("unaligned_fastq_seq ", str(len(unaligned_fastq_seq_list)))
            # with concurrent.futures.ThreadPoolExecutor() as executor:
            with concurrent.futures.ProcessPoolExecutor() as executor:
               # submit tasks and process results
               cs = max(10, args.n_processes_for_pooled // multiprocessing.cpu_count())
               for result in executor.map(self.get_new_variant_object, bg, unaligned_fastq_seq_list, unaligned_fq_seq_count_list, chunksize=cs):
                   # print(result)
                   new_variant,fq_seq = result

                   if new_variant == 0 or new_variant['best_match_score'] <= 10:
                       pass
                   else:
                       if fq_seq in self.variantCache and new_variant['best_match_score'] <= self.variantCache[fq_seq]['best_match_score']:
                           pass
                       else:
                           self.variantCache[fq_seq] = new_variant
            # print("how many fastq_seq are not aligned ", str(len(not_aln)))
        info("Processing reads Done")

        # if the sequence can't be aligned, skip it
        for fastq_seq, seq_count in fq_seq_dic.items():
            if fastq_seq in self.variantCache:
               # print(fastq_seq, self.variantCache[fastq_seq]['best_match_score'])
               if self.variantCache[fastq_seq]['best_match_score'] > 0.5*self.refs['sequence_length']:
                   N_COMPUTED_ALN += 1
                   N_CACHED_ALN +=  seq_count -1
               else:
                   N_COMPUTED_NOTALN += 1
                   N_CACHED_NOTALN +=  seq_count-1
            else:
                N_COMPUTED_NOTALN += 1
                N_CACHED_NOTALN +=  seq_count-1

        info("Finished reads; N_TOT_READS: %d N_COMPUTED_ALN: %d N_CACHED_ALN: %d N_COMPUTED_NOTALN: %d N_CACHED_NOTALN: %d"%(N_TOT_READS, N_COMPUTED_ALN, N_CACHED_ALN, N_COMPUTED_NOTALN, N_CACHED_NOTALN))
        self.aln_stats = {"N_TOT_READS" : N_TOT_READS,
               "N_CACHED_ALN" : N_CACHED_ALN,
               "N_CACHED_NOTALN" : N_CACHED_NOTALN,
               "N_COMPUTED_ALN" : N_COMPUTED_ALN,
               "N_COMPUTED_NOTALN" : N_COMPUTED_NOTALN
               }
        # print("aln_stats ")
        # print(self.aln_stats)
        # print("self.variantCache ", self.variantCache)


    def arrange_aligned_pair(self, compr, tmp_seq1, tmp_seq2):
        s1 = '' # fastq_seq
        s2 = '' # self.refs['sequence']
        score = 0

        # deal with unaligned region before the first aligned segment
        e = compr[0]
        # print(e[0][0], e[0][1], e[1][0], e[1][1])
        '''
        # padding method1
        if e[0][0] > e[1][0]:##deletion in fastq-seq
            s1 += "-"*(e[0][0] - e[1][0])+ tmp_seq1[:e[1][0]]
            s2 += tmp_seq2[:e[0][0]]
        elif e[0][0] < e[1][0]:
            s2 += "-"*(e[1][0] - e[0][0]) + tmp_seq2[:e[0][0]]
            s1 += tmp_seq1[:e[1][0]]
        else:
            s2 += tmp_seq2[:e[0][0]]
            s1 += tmp_seq1[:e[1][0]]
        '''
        # padding method2
        if e[0][0] != 0 or e[1][0] != 0:
            if e[0][0] == e[1][0] and e[0][0] <= 2:
                s1 += tmp_seq1[:e[1][0]]
                s2 += tmp_seq2[:e[0][0]]
            else:
                s1 += "-"*e[0][0] + tmp_seq1[:e[1][0]]
                s2 += tmp_seq2[:e[0][0]]+"-"*e[1][0]

        # print("s1 before match: ", s1)
        # print("s2 before match: ", s2)
        # start from the 1st alignment segments
        for idx in range(0,len(compr)-1):
            nxt_idx = idx + 1
            e = compr[idx]
            nxt_e = compr[nxt_idx]
            s2 += tmp_seq2[e[0][0]:e[0][1] + 1]
            s1 += tmp_seq1[e[1][0]:e[1][1] + 1]
            score += e[0][1] - e[0][0] +1

            # distance between two adjacent matched pairs
            int_len1 = nxt_e[0][0] -  e[0][1] -1
            int_len2 = nxt_e[1][0] -  e[1][1] -1

            '''
            # padding method1
            if int_len1 > int_len2:
                s1 += "-"*(int_len1 - int_len2)+ tmp_seq1[e[1][1]+1:nxt_e[1][0] ]
                s2 += tmp_seq2[e[0][1]+1:nxt_e[0][0]]
            elif int_len2 > int_len1:
                s2 += "-"*(int_len2 - int_len1)+ tmp_seq2[e[0][1]+1:nxt_e[0][0]]
                s1 += tmp_seq1[e[1][1]+1:nxt_e[1][0]]
            else:
                s2 += tmp_seq2[e[0][1]+1:nxt_e[0][0]]
                s1 += tmp_seq1[e[1][1]+1:nxt_e[1][0]]
            '''
            # padding method2
            if int_len1 == int_len2 and int_len1 <= 2:
                s1 += tmp_seq1[e[1][1]+1:nxt_e[1][0]]
                s2 += tmp_seq2[e[0][1]+1:nxt_e[0][0]]
            else:
                s1 += "-"*int_len1+ tmp_seq1[e[1][1]+1:nxt_e[1][0]]
                s2 += tmp_seq2[e[0][1]+1:nxt_e[0][0]]+"-"*int_len2

            # print("middle match ", str(idx), s1)
            # print("middle match ", str(idx), s2)

        # deal with the last paired match
        if len(compr) == 1:
            nxt_e = e

        s2 += tmp_seq2[nxt_e[0][0]:nxt_e[0][1] + 1]
        s1 += tmp_seq1[nxt_e[1][0]:nxt_e[1][1] + 1]
        score += nxt_e[0][1] - nxt_e[0][0] +1

        end_len2 = len(tmp_seq2) - nxt_e[0][1] -1
        end_len1 = len(tmp_seq1) - nxt_e[1][1] -1

        '''
        # padding method1
        if end_len1 > end_len2:
            s2 += "-"*(end_len1 - end_len2)+ tmp_seq2[nxt_e[0][1]+1:len(tmp_seq2)]
            s1 += tmp_seq1[nxt_e[1][1]+1:len(tmp_seq1)]
        elif end_len2 > end_len1:
            s1 += "-"*(end_len2 - end_len1)+ tmp_seq1[nxt_e[1][1]+1:len(tmp_seq1)]
            s2 += tmp_seq2[nxt_e[0][1]+1:len(tmp_seq2)]
        else:
            s2 += tmp_seq2[nxt_e[0][1]+1:len(tmp_seq2)]
            s1 += tmp_seq1[nxt_e[1][1]+1:len(tmp_seq1)]
        '''
        # padding method2
        if end_len1 == end_len2 and end_len1 <= 2:
            s1 += tmp_seq1[nxt_e[1][1]+1:]
            s2 += tmp_seq2[nxt_e[0][1]+1:]
        else:
            s1 += "-"*end_len2+ tmp_seq1[nxt_e[1][1]+1:]
            s2 += tmp_seq2[nxt_e[0][1]+1:]+"-"*end_len1

        return s1, s2, score
        # print("s1 = ", s1, "\ns2 = ", s2, "\nscore = ", score)

    def get_new_variant_object(self, amplicon_name, fastq_seq, fastq_seq_count):
        # get alignment and score from cython
        # score = 100 * # matchedBases / length(including gaps)
        aln_strand = '+'
        aln_seed_min = 2

        fq_len = len(fastq_seq)
        ref_len = self.refs['sequence_length']
        fq_seq_comp = reverse_complement(fastq_seq)
        # n_retry: length of reference sequence
        n_retry = ref_len

        # if fq_len >= ref_len:
        #    comp = IGDB_pairwise_aligner.CmpSeq(fastq_seq, self.refs['sequence'], n_retry)
        # else:
        #    print("ref_seq , fq_seq, fq_seq_comp ", self.refs['sequence'],fastq_seq, fq_seq_comp)

        # create an CmpSeq object
        if self.refs['sequence'] == fastq_seq:
            r = [((0, fq_len - 1), (0, fq_len - 1))]
            aln_strand = "+"
        elif self.refs['sequence'] == fq_seq_comp:
            r = [((0, fq_len - 1), (0, fq_len - 1))]
            aln_strand = "-"
        else:
            comp = IGDB_pairwise_aligner.CmpSeq(self.refs['sequence'], fastq_seq, n_retry)
            r, aln_strand = comp.cmp()

        if not r:
            print("fastq_seq ", fastq_seq)
            print("refer_seq ", self.refs["sequence"])
            print("Could not find match for this fastq sequence ")
            print("\n")
            return 0, fastq_seq

        '''
        print("fastq_seq ", fastq_seq)
        print("refer_seq ", self.refs["sequence"])
        for e in r:
            print("e:", e)
            print("ref: ", self.refs['sequence'][e[0][0]:e[0][1] + 1])
            print("faq: ",fastq_seq[e[1][0]:e[1][1] + 1])
        '''
        # print('r =',r)
        compr = []
        for item in r:
            compr.append(((item[0][0], item[0][1]), (item[1][0], item[1][1])))
        # compr = [((item[1][0], item[1][1]), (item[0][0], item[0][1])) if fq_len < ref_len else ((item[0][0], item[0][1]), (item[1][0], item[1][1])) for item in r]

        ############################################################################
        # call arrange_aligned_pair to arrange the alignment pair (consensus)
        tmp_seq2 = self.refs['sequence']
        if aln_strand == '-':
            tmp_seq1 = fq_seq_comp
        else:
            tmp_seq1 = fastq_seq
        s1, s2, score = self.arrange_aligned_pair(compr, tmp_seq1, tmp_seq2)
        '''
        print("s1 ", s1)
        print("s2 ", s2)
        print("score ", score)
        for e in compr:
            print(e)
        exit(0)
        '''
        # the 'min_aln_score' is calculated using only the changes in 'include_idxs'
        # start to get variations
        if score > 10:
            new_variant = {}
            new_variant['count'] = fastq_seq_count
            new_variant['aln_ref_names'] = amplicon_name
            new_variant['best_match_score'] = score
            class_names = []

            payload = CRISPRessoCOREResources.find_indels_substitutions(s1, s2, self.refs['include_idxs'])

            '''
            print("payload")
            for k, v in payload.items():
                print(k, v)
            exit(0)
            '''
            payload['ref_name'] = amplicon_name
 
            # if there is an insertion/deletion/substitution in the quantification window, the read is modified.
            is_modified = False
            if payload['deletion_n'] > 0:
                is_modified = True
            elif payload['insertion_n'] > 0:
                is_modified = True
            elif payload['substitution_n'] > 0:
                is_modified = True

            if is_modified:
                class_names.append(amplicon_name+"_MODIFIED")
                payload['classification'] = 'MODIFIED'
            else:
                class_names.append(amplicon_name+"_UNMODIFIED")
                payload['classification'] = 'UNMODIFIED'

            payload['aln_seq'] = s1
            payload['aln_ref'] = s2
            payload['aln_strand'] = aln_strand

            new_variant['variant_'+amplicon_name] = payload
            new_variant['class_name'] = "&".join(class_names)
        else:
            new_variant = {}
            new_variant['count'] = fastq_seq_count
            new_variant['best_match_score'] = score
            # return a new variant with the best match score set to 0, while retaining the scores of insufficient alignments
            return new_variant, fastq_seq
        return new_variant, fastq_seq


    def refObj_for_alignment(self, args):
        this_seq = self.amplicon_seq.strip().upper()
        this_seq_length = len(this_seq)

        include_idxs = []
        for idx in range(0, this_seq_length):
            include_idxs.append(idx)
        refObj = { 'sequence': this_seq,
                   'sequence_length': this_seq_length,
                   'include_idxs': include_idxs,
                   'left_primer_seq': args.left_primer_seq,
                   'right_primer_seq': args.right_primer_seq,
                   'aln_genome': None,
                   'aln_chr': None,
                   'aln_start': None,
                   'aln_end': None,
                   'aln_strand': None
            }
        # refs: dict with info for all refs
        self.refs = refObj


    def reads_amplicon_alignment(self, args):
        # reads quality control
        min_av_quality = None
        if args.min_average_read_quality > 0:
            min_av_quality = args.min_average_read_quality

        min_bp_quality_or_N = None
        if args.min_bp_quality_or_N > 0:
            min_bp_quality_or_N = args.min_bp_quality_or_N

        # print(args.merged_fastq)

        # count filtered reads
        self.N_READS_AFTER_PREPROCESSING=get_n_reads_fastq(args.merged_fastq)
        if self.N_READS_AFTER_PREPROCESSING == 0:
            raise Exception('No reads in input or no reads survived the average or single bp quality filtering.')

        info('Aligning sequences...')
        # INITIALIZE CACHE
        # put empty sequence into cache
        cache_fastq_seq = ''
        self.variantCache[cache_fastq_seq] = {}
        self.variantCache[cache_fastq_seq]['count'] = 0

        self.process_fastq(args.merged_fastq, args)
        '''
        print("self.aln_stats")
        for k, v in self.aln_stats.items():
        print(k, v)
        print("self.variantCache")
        for k, v in self.variantCache.items():
        print(k, v)
        end_time = time.perf_counter()
        runtime = end_time - start_time
        print(f"Time of program runing: {runtime}second")
        exit(0)
        '''
        info('Sequences Alignment Done!')

    def analyze_alignments(self, args):
        # ANALYZE ALIGNMENTS
        # class_counts = {} # number of reads in each class e.g. "ref1_UNMODIFIED" -> 50
        # alleles_list = [] # will be turned into df with rows with information for each variant (allele)

        # substitution_base_vectors = {} # these two are taken as a subset of all_substitution_base_vectors afterward to save time
        # base_count_vectors = {} # calculated after as well

        deleted_n_dicts = {}
        substituted_n_dicts = {}

        # initialize data structures for each ref
        this_len_amplicon = self.refs['sequence_length']

        self.counts_total                    = 0
        counts_modified                      = 0
        counts_unmodified                    = 0
        counts_discarded                     = 0

        counts_insertion                     = 0
        counts_deletion                      = 0
        counts_substitution                  = 0

        counts_only_insertion                = 0
        counts_only_deletion                 = 0
        counts_only_substitution             = 0

        counts_insertion_and_deletion        = 0
        counts_insertion_and_substitution    = 0
        counts_deletion_and_substitution     = 0
        counts_insertion_and_deletion_and_substitution   = 0

        # for each reference, the following are computed individually
        # count times substitutions occur
        # all_substitution_base_vectors = {}
        self.all_substitution_base_vectors        [args.amplicon_name+"_A" ] = np.zeros(this_len_amplicon)
        self.all_substitution_base_vectors        [args.amplicon_name+"_C" ] = np.zeros(this_len_amplicon)
        self.all_substitution_base_vectors        [args.amplicon_name+"_G" ] = np.zeros(this_len_amplicon)
        self.all_substitution_base_vectors        [args.amplicon_name+"_T" ] = np.zeros(this_len_amplicon)
        self.all_substitution_base_vectors        [args.amplicon_name+"_N" ] = np.zeros(this_len_amplicon)

        # all_base_count_vectors = {} # number of times each base is seen
        self.all_base_count_vectors               [args.amplicon_name+"_A" ] = np.zeros(this_len_amplicon)
        self.all_base_count_vectors               [args.amplicon_name+"_C" ] = np.zeros(this_len_amplicon)
        self.all_base_count_vectors               [args.amplicon_name+"_G" ] = np.zeros(this_len_amplicon)
        self.all_base_count_vectors               [args.amplicon_name+"_T" ] = np.zeros(this_len_amplicon)
        self.all_base_count_vectors               [args.amplicon_name+"_N" ] = np.zeros(this_len_amplicon)
        self.all_base_count_vectors               [args.amplicon_name+"_-" ] = np.zeros(this_len_amplicon)


        self.inserted_n_dicts                  = defaultdict(int)
        self.effective_len_dicts               = defaultdict(int)
        deleted_n_dicts                        = defaultdict(int)
        substituted_n_dicts                    = defaultdict(int)

        # take care of empty seqs
        cache_fastq_seq = ''
        self.variantCache[cache_fastq_seq]['count'] = 0

        # iterate through variants one by one
        for variant in self.variantCache:
            # skip variant if there is none observed
            variant_count = self.variantCache[variant]['count']
            if (variant_count == 0):
                continue

            # check to see if this sequence's reverse complement is in the variant
            rc_variant = reverse_complement(variant)
            if rc_variant in self.variantCache and self.variantCache[rc_variant]['count'] > 0:
                variant_count += self.variantCache[rc_variant]['count']
                self.variantCache[rc_variant]['count'] = 0
                self.variantCache[variant]['count'] = variant_count
            self.N_TOTAL += variant_count

            aln_ref_names = self.variantCache[variant]['aln_ref_names']
            aln_ref_scores = self.variantCache[variant]['best_match_score']

            class_name = self.variantCache[variant]['class_name'] #for classifying read e.g. 'HDR_MODIFIED' for pie chart
            if class_name not in self.class_counts:
                self.class_counts[class_name] = 0
            self.class_counts[class_name]+=variant_count

            # iterate through payloads -- if a read aligned equally-well to two references, it could have more than one payload
            variant_payload = self.variantCache[variant]["variant_"+aln_ref_names]
            self.counts_total += variant_count
            if variant_payload['classification'] == 'MODIFIED':
                counts_modified += variant_count
            else:
                counts_unmodified += variant_count

            write_detailed_allele_table = True
            allele_row = get_allele_row(aln_ref_names, variant_count, aln_ref_names, aln_ref_scores, variant_payload, write_detailed_allele_table)
            self.alleles_list.append(allele_row)

            # how long is this alignment (insertions increase length, deletions decrease length)
            this_effective_len= self.refs['sequence_length'] 

            # deal with insertions
            this_has_insertions = False
            self.all_insertion_count_vectors[variant_payload['all_insertion_positions']]+=variant_count
            self.all_insertion_left_count_vectors[variant_payload['all_insertion_left_positions']]+=variant_count

            self.inserted_n_dicts[variant_payload['insertion_n']] += variant_count
            self.insertion_count_vectors[variant_payload['insertion_positions']]+=variant_count
            this_effective_len = this_effective_len + variant_payload['insertion_n']
            if variant_payload['insertion_n'] > 0:
                counts_insertion += variant_count
                this_has_insertions = True

            # deal with deletions
            this_has_deletions = False
            self.all_deletion_count_vectors[variant_payload['all_deletion_positions']]+=variant_count
            deleted_n_dicts[variant_payload['deletion_n']] += variant_count
            self.deletion_count_vectors[variant_payload['deletion_positions']]+=variant_count
            this_effective_len = this_effective_len - variant_payload['deletion_n']
            if variant_payload['deletion_n'] > 0:
                counts_deletion += variant_count
                this_has_deletions = True

            self.effective_len_dicts[this_effective_len] += variant_count

            # deal with substitutions
            this_has_substitutions = False
            self.all_substitution_count_vectors[variant_payload['all_substitution_positions']] += variant_count

            substituted_n_dicts[variant_payload['substitution_n']] += variant_count
            self.substitution_count_vectors[variant_payload['substitution_positions']] += variant_count
            if variant_payload['substitution_n'] > 0:
                counts_substitution += variant_count
                this_has_substitutions = True

            nucs = ['A', 'T', 'C', 'G', 'N']
            for nuc in nucs:
                isNuc = [n == nuc for n in variant_payload['all_substitution_values']]
                if(np.sum(isNuc) > 0):
                    locs = np.array(variant_payload['all_substitution_positions'])[isNuc]
                    self.all_substitution_base_vectors[aln_ref_names + "_" + nuc ][locs] += variant_count

            # start to count "and"
            if this_has_deletions:
                if this_has_insertions:
                    if this_has_substitutions:
                        counts_insertion_and_deletion_and_substitution += variant_count
                    else:
                        counts_insertion_and_deletion += variant_count
                else:
                    if this_has_substitutions:
                        counts_deletion_and_substitution += variant_count
                    else:
                        counts_only_deletion += variant_count
            else: # no deletions
                if this_has_insertions:
                    if this_has_substitutions:
                        counts_insertion_and_substitution += variant_count
                    else:
                        counts_only_insertion += variant_count
                else:
                    if this_has_substitutions:
                        counts_only_substitution += variant_count

            # set all_base_count_vectors
            aln_seq = variant_payload['aln_seq']
            ref_pos = variant_payload['ref_positions']

            for i in range(len(aln_seq)):
                if ref_pos[i] < 0:
                    continue
                nuc = aln_seq[i]
                self.all_base_count_vectors[aln_ref_names + "_" + nuc][ref_pos[i]] += variant_count

            if this_has_insertions or this_has_deletions or this_has_substitutions: #only count modified reads
                insertion_coordinates = variant_payload['insertion_coordinates']
                insertion_sizes = variant_payload['insertion_sizes']
                all_insertion_positions = variant_payload['all_insertion_positions']
                all_insertion_left_positions = variant_payload['all_insertion_left_positions']
                insertion_positions = variant_payload['insertion_positions']
                deletion_coordinates = variant_payload['deletion_coordinates']
                deletion_sizes = variant_payload['deletion_sizes']
                all_deletion_positions = variant_payload['all_deletion_positions']
                deletion_positions = variant_payload['deletion_positions']
                all_substitution_positions = variant_payload['all_substitution_positions']
                substitution_positions = variant_payload['substitution_positions']

                length_modified_positions_exons=[]
                current_read_exons_modified = False
                current_read_spliced_modified = False

                for idx_ins, (ins_start, ins_end) in enumerate(insertion_coordinates):
                    self.insertion_length_vectors[ins_start]+=(insertion_sizes[idx_ins]*variant_count)
                    self.insertion_length_vectors[ins_end]+=(insertion_sizes[idx_ins]*variant_count)

                for idx_del, (del_start, del_end) in enumerate(deletion_coordinates):
                    self.deletion_length_vectors[list(range(del_start, del_end))] += (deletion_sizes[idx_del]*variant_count)
        # done iterating through variantCache objects

        this_include_idx = self.refs['include_idxs']
        aln_ref_names = args.amplicon_name
        self.substitution_base_vectors      [aln_ref_names+"_A" ] = [self.all_substitution_base_vectors[aln_ref_names+"_A"][x] for x in this_include_idx]
        self.substitution_base_vectors      [aln_ref_names+"_C" ] = [self.all_substitution_base_vectors[aln_ref_names+"_C"][x] for x in this_include_idx]
        self.substitution_base_vectors      [aln_ref_names+"_G" ] = [self.all_substitution_base_vectors[aln_ref_names+"_G"][x] for x in this_include_idx]
        self.substitution_base_vectors      [aln_ref_names+"_T" ] = [self.all_substitution_base_vectors[aln_ref_names+"_T"][x] for x in this_include_idx]
        self.substitution_base_vectors      [aln_ref_names+"_N" ] = [self.all_substitution_base_vectors[aln_ref_names+"_N"][x] for x in this_include_idx]

        self.base_count_vectors             [aln_ref_names+"_A" ] = [self.all_base_count_vectors[aln_ref_names+"_A"][x] for x in this_include_idx]
        self.base_count_vectors             [aln_ref_names+"_C" ] = [self.all_base_count_vectors[aln_ref_names+"_C"][x] for x in this_include_idx]
        self.base_count_vectors             [aln_ref_names+"_G" ] = [self.all_base_count_vectors[aln_ref_names+"_G"][x] for x in this_include_idx]
        self.base_count_vectors             [aln_ref_names+"_T" ] = [self.all_base_count_vectors[aln_ref_names+"_T"][x] for x in this_include_idx]
        self.base_count_vectors             [aln_ref_names+"_N" ] = [self.all_base_count_vectors[aln_ref_names+"_N"][x] for x in this_include_idx]
        self.base_count_vectors             [aln_ref_names+"_-" ] = [self.all_base_count_vectors[aln_ref_names+"_-"][x] for x in this_include_idx]

        self.all_indelsub_count_vectors = self.all_insertion_count_vectors + self.all_deletion_count_vectors + self.all_substitution_count_vectors
        indelsub_count_vectors = self.insertion_count_vectors + self.deletion_count_vectors + self.substitution_count_vectors
 
        # order class_counts
        class_counts_order = []
        for class_count_name in self.class_counts:
            # print(" class_count_name: ",  class_count_name)
            # for idx, ref_name in enumerate(ref_names):
            class_counts_order.append(class_count_name)

        if self.N_TOTAL == 0:
            raise Exception('Error: No alignments were found')

        # create alleles table
        info('Calculating allele frequencies...')


        # set up allele table
        df_alleles = pd.DataFrame(self.alleles_list)
        # df_alleles['%Reads']=df_alleles['#Reads']/df_alleles['#Reads'].sum()*100 # sum of #reads will be >= self.N_TOTAL because an allele appears once for each reference it aligns to
        df_alleles['%Reads']=df_alleles['#Reads']/self.N_TOTAL*100
        df_alleles[['n_deleted', 'n_inserted', 'n_mutated']] = df_alleles[['n_deleted', 'n_inserted', 'n_mutated']].astype(int)

        df_alleles.sort_values(by='#Reads', ascending=False, inplace=True)


        all_insertion_pct_vectors = {} # all insertions/total (including quantification window bases)
        all_deletion_pct_vectors = {}
        all_substitution_pct_vectors = {}
        all_indelsub_pct_vectors = {}

        insertion_pct_vectors = {} # insertions that are in the quantification window
        deletion_pct_vectors = {}
        substitution_pct_vectors = {}
        indelsub_pct_vectors = {}


        # for ref_name in ref_names:
        # save these values in the ref object-- we need to print them later

        ref_len = self.refs['sequence_length']
        min_cut=ref_len/2
        max_cut=ref_len/2
        xmin, xmax=-min_cut, +max_cut

        self.refs['min_cut'] = min_cut
        self.refs['max_cut'] = max_cut

        max_mut=max(15, max(substituted_n_dicts.keys() or [0])) # the or bit is for when there are no keys
        max_ins=max(15, max(self.inserted_n_dicts.keys() or [0]))
        max_del=max(15, max(deleted_n_dicts.keys() or [0]))

        x_bins_mut = np.arange(max_mut+1)
        y_values_mut = np.array([substituted_n_dicts[x] for x in x_bins_mut])
        x_bins_ins = np.arange(max_ins+1)
        y_values_ins = np.array([self.inserted_n_dicts[x] for x in x_bins_ins])
        x_bins_del = np.arange(max_del+1)
        y_values_del = np.array([deleted_n_dicts[x] for x in x_bins_del])

        self.refs['y_values_mut'] = y_values_mut
        self.refs['x_bins_mut'] = x_bins_mut
        self.refs['y_values_ins'] = y_values_ins
        self.refs['x_bins_ins'] = x_bins_ins
        self.refs['y_values_del'] = y_values_del
        self.refs['x_bins_del'] = x_bins_del

        min_eff_len = min(ref_len-15, min(self.effective_len_dicts.keys() or [0]))
        max_eff_len = max(ref_len+15, max(self.effective_len_dicts.keys() or [0]))
        hlengths = np.arange(min_eff_len, max_eff_len+1)

        hdensity = np.array([self.effective_len_dicts[x] for x in hlengths])
        hlengths = hlengths-ref_len
        center_index=np.where(hlengths==0)[0][0]

        self.refs['hdensity'] = hdensity
        self.refs['hlengths'] = hlengths
        self.refs['center_index'] = center_index

        count_tot = self.counts_total
        if count_tot > 0:
            # normalize effect vectors
            all_insertion_pct_vectors = 100 * self.all_insertion_count_vectors / count_tot
            all_deletion_pct_vectors = 100 * self.all_deletion_count_vectors / count_tot
            all_substitution_pct_vectors = 100 * self.all_substitution_count_vectors / count_tot
            all_indelsub_pct_vectors = 100 * self.all_indelsub_count_vectors / count_tot

            insertion_pct_vectors = 100 * self.insertion_count_vectors / count_tot
            deletion_pct_vectors = 100 * self.deletion_count_vectors / count_tot
            substitution_pct_vectors = 100 * self.substitution_count_vectors / count_tot
            indelsub_pct_vectors = 100 * self.indelsub_count_vectors / count_tot

            # insertion_pct_vectors_noncoding[ref_name] = 100 * insertion_count_vectors_noncoding[ref_name] / count_tot
            # deletion_pct_vectors_noncoding[ref_name] = 100 * deletion_count_vectors_noncoding[ref_name] / count_tot
            # substitution_pct_vectors_noncoding[ref_name] = 100 * substitution_count_vectors_noncoding[ref_name] / count_tot

            deletion_length_vectors_sum = self.deletion_length_vectors
            self.deletion_length_vectors = np.zeros(ref_len)
            delMask = self.deletion_count_vectors > 0
            self.deletion_length_vectors[delMask] = deletion_length_vectors_sum[delMask]/self.deletion_count_vectors[delMask]

            insertion_length_vectors_sum = self.insertion_length_vectors
            self.insertion_length_vectors = np.zeros(ref_len)
            insMask = self.insertion_count_vectors > 0
            self.insertion_length_vectors[insMask] = insertion_length_vectors_sum[insMask]/self.insertion_count_vectors[insMask]

        else: # no reads for this ref
            this_len_amplicon = self.refs['sequence_length']
            all_insertion_pct_vectors = np.zeros(this_len_amplicon)
            all_deletion_pct_vectors  = np.zeros(this_len_amplicon)
            all_substitution_pct_vectors = np.zeros(this_len_amplicon)
            all_indelsub_pct_vectors = np.zeros(this_len_amplicon)

            insertion_pct_vectors = np.zeros(this_len_amplicon)
            deletion_pct_vectors = np.zeros(this_len_amplicon)
            substitution_pct_vectors = np.zeros(this_len_amplicon)
            indelsub_pct_vectors = np.zeros(this_len_amplicon)

            # insertion_pct_vectors_noncoding[ref_name] = np.zeros(this_len_amplicon)
            # deletion_pct_vectors_noncoding[ref_name] = np.zeros(this_len_amplicon)
            # substitution_pct_vectors_noncoding[ref_name] = np.zeros(this_len_amplicon)

            self.deletion_length_vectors = np.zeros(ref_len)
            self.insertion_length_vectors = np.zeros(ref_len)

        info('Done!')
        self.crispresso2_info['results']['ref_name'] = aln_ref_names
        self.crispresso2_info['results']['refs'] = self.refs

        info('Saving processed data...')

        # write alleles table
        # crispresso1Cols = ["Aligned_Sequence","Reference_Sequence","NHEJ","UNMODIFIED","HDR","n_deleted","n_inserted","n_mutated","#Reads","%Reads"]
        # df_alleles.loc[:,crispresso1Cols].to_csv(_jp('Alleles_frequency_table.txt'),sep='\t',header=True,index=None)
        # crispresso2Cols = ["Aligned_Sequence","Reference_Sequence","Reference_Name","Read_Status","n_deleted","n_inserted","n_mutated","#Reads","%Reads"]
        # crispresso2Cols = ["Aligned_Sequence","Reference_Sequence","Reference_Name","Read_Status","n_deleted","n_inserted","n_mutated","#Reads","%Reads","Aligned_Reference_Names","Aligned_Reference_Scores"]
        # crispresso2Cols = ["Read_Sequence","Amplicon_Sequence","Amplicon_Name","Read_Status","n_deleted","n_inserted","n_mutated","#Reads","%Reads"]
        crispresso2Cols = ["Aligned_Sequence", "Reference_Sequence", "Reference_Name", "Read_Status", "n_deleted", "n_inserted", "n_mutated", "#Reads", "%Reads"]

        # allele_frequency_table_filename = 'Alleles_frequency_table.txt'
        allele_frequency_table_filename = 'Alleles_frequency_table.csv'
        allele_frequency_table_fileLoc = _jp(self.OUTPUT_DIRECTORY, allele_frequency_table_filename)

        allele_frequency_table_zip_filename = _jp(self.OUTPUT_DIRECTORY, 'Alleles_frequency_table.zip')

        # df_alleles.loc[:, crispresso2Cols].to_csv(allele_frequency_table_fileLoc, sep='\t', header=True, index=None)
        df_alleles.loc[:, crispresso2Cols].to_csv(allele_frequency_table_fileLoc, sep=',', header=True, index=None)
        with zipfile.ZipFile(allele_frequency_table_zip_filename, 'w', zipfile.ZIP_DEFLATED, allowZip64=True) as myzip:
            myzip.write(allele_frequency_table_fileLoc, allele_frequency_table_filename)
        os.remove(allele_frequency_table_fileLoc)
        self.crispresso2_info['running_info']['allele_frequency_table_filename'] = os.path.basename(allele_frequency_table_filename) # filename is the name of the file in zip
        self.crispresso2_info['running_info']['allele_frequency_table_zip_filename'] = os.path.basename(allele_frequency_table_zip_filename)
        quant_of_editing_freq_filename =_jp(self.OUTPUT_DIRECTORY, 'CRISPResso_quantification_of_editing_frequency.txt')

        with open(quant_of_editing_freq_filename, 'w+') as outfile:
            outfile.write('Amplicon\tUnmodified%\tModified%\tReads_in_input\tReads_aligned_all_amplicons\tReads_aligned\tUnmodified\tModified\tDiscarded\tInsertions\tDeletions\tSubstitutions\tOnly Insertions\tOnly Deletions\tOnly Substitutions\tInsertions and Deletions\tInsertions and Substitutions\tDeletions and Substitutions\tInsertions Deletions and Substitutions\n')
            # for ref_name in ref_names:
            n_aligned = self.counts_total
            n_unmod = counts_unmodified
            n_mod = counts_modified
            n_discarded = counts_discarded

            n_insertion = counts_insertion
            n_deletion = counts_deletion
            n_substitution = counts_substitution
            n_only_insertion = counts_only_insertion
            n_only_deletion = counts_only_deletion
            n_only_substitution = counts_only_substitution
            n_insertion_and_deletion = counts_insertion_and_deletion
            n_insertion_and_substitution = counts_insertion_and_substitution
            n_deletion_and_substitution = counts_deletion_and_substitution
            n_insertion_and_deletion_and_substitution = counts_insertion_and_deletion_and_substitution

            unmod_pct = "NA"
            mod_pct = "NA"
            if n_aligned > 0:
                unmod_pct = round(100*n_unmod/float(n_aligned), 8)
                mod_pct = round(100*n_mod/float(n_aligned), 8)

            vals = [aln_ref_names]
            vals.extend([str(x) for x in [unmod_pct, mod_pct, self.N_READS_INPUT, self.N_TOTAL, n_aligned, n_unmod, n_mod, n_discarded, n_insertion, n_deletion, n_substitution, n_only_insertion, n_only_deletion, n_only_substitution, n_insertion_and_deletion, n_insertion_and_substitution, n_deletion_and_substitution, n_insertion_and_deletion_and_substitution]])
            outfile.write("\t".join(vals) + "\n")

        self.crispresso2_info['running_info']['quant_of_editing_freq_filename'] = os.path.basename(quant_of_editing_freq_filename)

        mapping_stats_filename = _jp(self.OUTPUT_DIRECTORY, 'CRISPResso_mapping_statistics.txt')
        with open(mapping_stats_filename, 'w+') as outfile:
            outfile.write('READS IN INPUTS\tREADS AFTER PREPROCESSING\tREADS ALIGNED\tN_COMPUTED_ALN\tN_CACHED_ALN\tN_COMPUTED_NOTALN\tN_CACHED_NOTALN\n')
            outfile.write("\t".join([str(x) for x in[self.N_READS_INPUT, self.N_READS_AFTER_PREPROCESSING, self.N_TOTAL, self.aln_stats['N_COMPUTED_ALN'], self.aln_stats['N_CACHED_ALN'], self.aln_stats['N_COMPUTED_NOTALN'], self.aln_stats['N_CACHED_NOTALN']]]) + "\n")
        self.crispresso2_info['running_info']['alignment_stats'] = self.aln_stats
        self.crispresso2_info['running_info']['mapping_stats_filename'] = os.path.basename(mapping_stats_filename)

        self.crispresso2_info['results']['alignment_stats']['insertion_pct_vectors'] = insertion_pct_vectors
        self.crispresso2_info['results']['alignment_stats']['deletion_pct_vectors'] = insertion_pct_vectors
        self.crispresso2_info['results']['alignment_stats']['substitution_pct_vectors'] = insertion_pct_vectors
        self.crispresso2_info['results']['alignment_stats']['indelsub_pct_vectors'] = insertion_pct_vectors

        # set unique plot name to appear as prefix to files for each reference
        seen_ref_names = {}

        # only show reference name in filenames if more than one reference
        ref_plot_name =  aln_ref_names
        ref_plot_name = ""
        seen_ref_names[ref_plot_name] = 1
        self.refs['ref_plot_name'] = ref_plot_name

        ref_plot_name = self.refs['ref_plot_name']

        # n_this_category = self.counts_total
        # if n_this_category < 1:
        #    continue

        ins_pct_vector_filename = _jp(self.OUTPUT_DIRECTORY, ref_plot_name+'Effect_vector_insertion.txt')
        save_vector_to_file(insertion_pct_vectors, ins_pct_vector_filename)
        self.crispresso2_info['results']['refs']['insertion_pct_vector_filename'] = os.path.basename(ins_pct_vector_filename)

        del_pct_vector_filename = _jp(self.OUTPUT_DIRECTORY, ref_plot_name+'Effect_vector_deletion.txt')
        save_vector_to_file(deletion_pct_vectors, del_pct_vector_filename)
        self.crispresso2_info['results']['refs']['deletion_pct_vector_filename'] = os.path.basename(del_pct_vector_filename)

        sub_pct_vector_filename = _jp(self.OUTPUT_DIRECTORY, ref_plot_name+'Effect_vector_substitution.txt')
        save_vector_to_file(substitution_pct_vectors, sub_pct_vector_filename)
        self.crispresso2_info['results']['refs']['substitution_pct_vector_filename'] = os.path.basename(sub_pct_vector_filename)

        indelsub_pct_vector_filename = _jp(self.OUTPUT_DIRECTORY, ref_plot_name+'Effect_vector_combined.txt')
        save_vector_to_file(indelsub_pct_vectors, indelsub_pct_vector_filename)
        self.crispresso2_info['results']['refs']['combined_pct_vector_filename'] = os.path.basename(indelsub_pct_vector_filename)

        # save mods in quantification window
        quant_window_mod_count_filename = _jp(self.OUTPUT_DIRECTORY, 'Quantification_window_modification_count_vectors.txt')
        save_count_vectors_to_file([self.insertion_count_vectors,
                    self.deletion_count_vectors,
                    self.substitution_count_vectors,
                    self.indelsub_count_vectors,
                    [self.counts_total]*self.refs['sequence_length']],
                    ['Insertions', 'Deletions', 'Substitutions', 'All_modifications', 'Total'],
                        self.refs['sequence'], quant_window_mod_count_filename)
        self.crispresso2_info['results']['refs']['quant_window_mod_count_filename'] = os.path.basename(quant_window_mod_count_filename)

        # save all mods
        mod_count_filename = _jp(self.OUTPUT_DIRECTORY, ref_plot_name+'Modification_count_vectors.txt')
        save_count_vectors_to_file([self.all_insertion_count_vectors,
                    self.all_insertion_left_count_vectors,
                    self.all_deletion_count_vectors,
                    self.all_substitution_count_vectors,
                    self.all_indelsub_count_vectors,
                    [self.counts_total]*self.refs['sequence_length']],
                    ['Insertions', 'Insertions_Left', 'Deletions', 'Substitutions', 'All_modifications', 'Total'],
                    self.refs['sequence'], mod_count_filename)
        self.crispresso2_info['results']['refs']['mod_count_filename'] = os.path.basename(mod_count_filename)
        self.crispresso2_info['results']['refs']['mod_count_filename_caption'] = "A tab-separated file showing the number of modifications for each position in the amplicon. " \
            "The first row shows the amplicon sequence, and successive rows show the number of reads with insertions (row 2), insertions_left (row 3), deletions (row 4), substitutions (row 5) and the sum of all modifications (row 6)." \
            "Additionally, the last row shows the number of reads aligned. If an insertion occurs between bases 5 and 6, the insertions vector will be incremented at bases 5 and 6. " \
            "However, the insertions_left vector will only be incremented at base 5 so the sum of the insertions_left row represents an accurate count of the number of insertions, " \
            "whereas the sum of the insertions row will yield twice the number of insertions."


        hdensity = self.refs['hdensity']
        hlengths = self.refs['hlengths']
        center_index = self.refs['center_index']

        y_values_mut = self.refs['y_values_mut']
        x_bins_mut = self.refs['x_bins_mut']
        y_values_ins = self.refs['y_values_ins']
        x_bins_ins = self.refs['x_bins_ins']
        y_values_del = self.refs['y_values_del']
        x_bins_del = self.refs['x_bins_del']

        # if not args.suppress_plots:
        indel_histogram_file = _jp(self.OUTPUT_DIRECTORY, ref_plot_name+'Indel_histogram.txt')
        pd.DataFrame({'indel_size':pd.Series(hlengths, dtype='int'),'fq':hdensity})[['indel_size', 'fq']].to_csv(indel_histogram_file, index=None, sep='\t')
        self.crispresso2_info['results']['refs']['indel_histogram_filename'] = os.path.basename(indel_histogram_file)
        self.crispresso2_info['results']['refs']['indel_histogram_filename_caption'] = "A tab-separated text file that shows a histogram of the length of indels (both insertions and deletions) in the " + aln_ref_names +" sequence in the quantification window. " \
            "Indels outside of the quantification window are not included. The indel_size column shows the number of substitutions, and the fq column shows the number of reads having an indel of that length."


        insertion_histogram_file = _jp(self.OUTPUT_DIRECTORY, ref_plot_name+'Insertion_histogram.txt')
        pd.DataFrame({'ins_size':x_bins_ins,'fq':y_values_ins})[["ins_size", "fq"]].to_csv(insertion_histogram_file, index=None, sep='\t')
        self.crispresso2_info['results']['refs']['insertion_histogram_filename'] = os.path.basename(insertion_histogram_file)
        self.crispresso2_info['results']['refs']['insertion_histogram_filename_caption'] = "A tab-separated text file that shows a histogram of the number of insertions in the " + aln_ref_names +" sequence in the quantification window. " \
            "Insertions outside of the quantification window are not included. The ins_size column shows the number of insertions, and the fq column shows the number of reads having that number of insertions."


        deletion_histogram_file = _jp(self.OUTPUT_DIRECTORY, ref_plot_name+'Deletion_histogram.txt')
        pd.DataFrame({'del_size':-1*x_bins_del,'fq':y_values_del})[['del_size', 'fq']].to_csv(deletion_histogram_file, index=None, sep='\t')
        self.crispresso2_info['results']['refs']['deletion_histogram_filename'] = os.path.basename(deletion_histogram_file)
        self.crispresso2_info['results']['refs']['deletion_histogram_filename_caption'] = "A tab-separated text file that shows a histogram of the number of deletions in the " + aln_ref_names +" sequence in the quantification window. " \
            "Deletions outside of the quantification window are not included. The del_size column shows the number of deletions, and the fq column shows the number of reads having that number of deletions."

        substitution_histogram_file = _jp(self.OUTPUT_DIRECTORY, ref_plot_name+'Substitution_histogram.txt')
        pd.DataFrame({'sub_count':x_bins_mut,'fq':y_values_mut})[['sub_count', 'fq']].to_csv(substitution_histogram_file, index=None, sep='\t')
        self.crispresso2_info['results']['refs']['substitution_histogram_filename'] = os.path.basename(substitution_histogram_file)
        self.crispresso2_info['results']['refs']['substitution_histogram_filename_caption'] = "A tab-separated text file that shows a histogram of the number of substitutions in the " + aln_ref_names +" sequence in the quantification window. " \
            "Substitutions outside of the quantification window are not included. The sub_size column shows the number of substitutions, and the fq column shows the number of reads having that number of substitutions."


    def making_alignment_plot_figure1(self, args):
        # FIGURE 1: Alignment
        n_processes = args.n_processes_for_pooled
        plot_pool = concurrent.futures.ProcessPoolExecutor(n_processes)
        plot_results = []
        ref_name = args.amplicon_name
        save_png = True
        ref_plot_name = ref_name+'_'

        plot_1a_root = _jp(self.OUTPUT_DIRECTORY, "1a.Read_barplot")
        plot_1a_input = {
            'N_READS_INPUT': self.N_READS_INPUT,
            'N_READS_AFTER_PREPROCESSING': self.N_READS_AFTER_PREPROCESSING,
            'N_TOTAL': self.N_TOTAL,
            'plot_root': plot_1a_root,
            'save_png': save_png
        }
        if n_processes > 1:
            plot_results.append(plot_pool.submit(
                CRISPRessoPlot.plot_read_barplot,
                **plot_1a_input,
            ))
        else:
            CRISPRessoPlot.plot_read_barplot(**plot_1a_input)
        # plot 1a and 1b using CRISPResso_mapping_statistics.txt file
        mapping_stats_filename = _jp(self.OUTPUT_DIRECTORY, 'CRISPResso_mapping_statistics.txt')
        self.crispresso2_info['results']['general_plots']['plot_1a_root'] = os.path.basename(plot_1a_root)
        self.crispresso2_info['results']['general_plots']['plot_1a_caption'] = "Figure 1a: The number of reads in input fastqs, after preprocessing, and after alignment to amplicons."
        self.crispresso2_info['results']['general_plots']['plot_1a_data'] = [('Mapping statistics', os.path.basename(mapping_stats_filename))]


        # plot 1b using CRISPResso_quantification_of_editing_frequency.txt file
        class_counts_order = []
        for class_count_name in self.class_counts:
            class_counts_order.append(class_count_name)

        quant_of_editing_freq_filename =_jp(self.OUTPUT_DIRECTORY, 'CRISPResso_quantification_of_editing_frequency.txt')
        plot_1b_root = _jp(self.OUTPUT_DIRECTORY, "1b.Alignment_pie_chart")
        plot_1c_root = _jp(self.OUTPUT_DIRECTORY, '1c.Alignment_barplot')
        plot_1bc_input = {
            'class_counts_order': class_counts_order,
            'class_counts': self.class_counts,
            'ref_names': ref_name,
            'expected_hdr_amplicon_seq': "",
            'N_TOTAL': self.N_TOTAL,
            'piechart_plot_root': plot_1b_root,
            'barplot_plot_root': plot_1c_root,
            'save_png': save_png
        }
        self.crispresso2_info['results']['general_plots']['plot_1b_root'] = os.path.basename(plot_1b_root)
        self.crispresso2_info['results']['general_plots']['plot_1b_caption'] = "Figure 1b: Alignment and editing frequency of reads as determined by the percentage and number of sequence reads showing unmodified and modified alleles."

        self.crispresso2_info['results']['general_plots']['plot_1b_data'] = [('Quantification of editing', os.path.basename(quant_of_editing_freq_filename))]

        self.crispresso2_info['results']['general_plots']['plot_1c_root'] = os.path.basename(plot_1c_root)
        self.crispresso2_info['results']['general_plots']['plot_1c_caption'] = "Figure 1c: Alignment and editing frequency of reads as determined by the percentage and number of sequence reads showing unmodified and modified alleles."
        self.crispresso2_info['results']['general_plots']['plot_1c_data'] = [('Quantification of editing', os.path.basename(quant_of_editing_freq_filename))]
        if n_processes > 1:
            plot_results.append(plot_pool.submit(
                CRISPRessoPlot.plot_class_piechart_and_barplot,
                **plot_1bc_input,
            ))
        else:
            CRISPRessoPlot.plot_class_piechart_and_barplot(
                **plot_1bc_input,
            )
        # to test, run: plot_pool.apply_async(CRISPRessoPlot.plot_class_piechart_and_barplot, kwds=plot_1bc_input).get()


    def making_nuc_pct_plot_figure2(self, args):
        # Figure 2 nucleotide percentage
        n_processes = args.n_processes_for_pooled
        plot_pool = concurrent.futures.ProcessPoolExecutor(n_processes)
        plot_results = []
        ref_name = args.amplicon_name
        save_png = True

        # hold mod pct dfs for each amplicon
        mod_pct_dfs = {}
        ref_len = self.refs['sequence_length']
        ref_seq = self.refs['sequence']
        min_cut = self.refs['min_cut']
        max_cut = self.refs['max_cut']
        hdensity = self.refs['hdensity']
        hlengths = self.refs['hlengths']
        center_index = self.refs['center_index']
        include_idxs_list = self.refs['include_idxs']
        quantification_window_ref_seq = [list(ref_seq)[x] for x in include_idxs_list]

        tot_aln_reads = self.counts_total
        n_this_category = self.counts_total

        # only show reference name in filenames if more than one reference
        # ref_plot_name = self.refs['ref_plot_name']
        ref_plot_name = ref_name+'_'

        if n_this_category < 1:
            info("The total aligned reads are 0\n")
            sys.exit(0)

        # plot quilt for this amplicon  (if not crispresso1 mode)
        # nucleotide counts
        df_nuc_freq = pd.DataFrame([self.base_count_vectors[ref_name+"_A"], self.base_count_vectors[ref_name+"_C"], self.base_count_vectors[ref_name+"_G"], self.base_count_vectors[ref_name+"_T"], self.base_count_vectors[ref_name+"_N"], self.base_count_vectors[ref_name+'_-']])
        df_nuc_freq.index = ['A', 'C', 'G', 'T', 'N', '-']
        df_nuc_freq.columns = quantification_window_ref_seq
        # print table showing nuc frequencies (sum to total alleles) (in quantification window)
        quant_window_nuc_freq_filename = _jp(self.OUTPUT_DIRECTORY, ref_plot_name + 'Quantification_window_nucleotide_frequency_table.txt')
        df_nuc_freq.to_csv(quant_window_nuc_freq_filename, sep='\t', header=True, index=True)
        self.crispresso2_info['results']['refs']['quant_window_nuc_freq_filename'] = os.path.basename(quant_window_nuc_freq_filename)

        df_nuc_pct = df_nuc_freq.divide(tot_aln_reads)
        quant_window_nuc_pct_filename = _jp(self.OUTPUT_DIRECTORY, ref_plot_name + 'Quantification_window_nucleotide_percentage_table.txt')
        df_nuc_pct.to_csv(quant_window_nuc_pct_filename, sep='\t', header=True, index=True)
        self.crispresso2_info['results']['refs']['quant_window_nuc_pct_filename'] = os.path.basename(quant_window_nuc_pct_filename)

        df_nuc_freq_all = pd.DataFrame([self.all_base_count_vectors[ref_name+"_A"], self.all_base_count_vectors[ref_name+"_C"], self.all_base_count_vectors[ref_name+"_G"], self.all_base_count_vectors[ref_name+"_T"], self.all_base_count_vectors[ref_name+"_N"], self.all_base_count_vectors[ref_name+'_-']])
        df_nuc_freq_all.index = ['A', 'C', 'G', 'T', 'N', '-']
        df_nuc_freq_all.columns = list(ref_seq)
        # print table showing nuc frequencies (sum to total alleles) (in entire region)
        nuc_freq_filename = _jp(self.OUTPUT_DIRECTORY, ref_plot_name + 'Nucleotide_frequency_table.txt')
        df_nuc_freq_all.to_csv(nuc_freq_filename, sep='\t', header=True, index=True)
        self.crispresso2_info['results']['refs']['nuc_freq_filename'] = os.path.basename(nuc_freq_filename)

        df_nuc_pct_all = df_nuc_freq_all.divide(tot_aln_reads)
        nuc_pct_filename = _jp(self.OUTPUT_DIRECTORY, ref_plot_name + 'Nucleotide_percentage_table.txt')
        df_nuc_pct_all.to_csv(nuc_pct_filename, sep='\t', header=True, index=True)
        self.crispresso2_info['results']['refs']['nuc_pct_filename'] = os.path.basename(nuc_pct_filename)

        mod_pcts = []

        tot = float(self.counts_total)
        mod_pcts.append(np.concatenate((['Insertions'], np.array(self.all_insertion_count_vectors).astype(np.float64)/tot)))
        mod_pcts.append(np.concatenate((['Insertions_Left'], np.array(self.all_insertion_left_count_vectors).astype(np.float64)/tot)))
        mod_pcts.append(np.concatenate((['Deletions'], np.array(self.all_deletion_count_vectors).astype(np.float64)/tot)))
        mod_pcts.append(np.concatenate((['Substitutions'], np.array(self.all_substitution_count_vectors).astype(np.float64)/tot)))
        mod_pcts.append(np.concatenate((['All_modifications'], np.array(self.all_indelsub_count_vectors).astype(np.float64)/tot)))
        mod_pcts.append(np.concatenate((['Total'], [self.counts_total]*self.refs['sequence_length'])))

        # print("mod_pcts \n", mod_pcts)
        colnames = ['Modification']+list(ref_seq)
        # modification_percentage_summary_df = pd.DataFrame(mod_pcts, columns=colnames).apply(pd.to_numeric, errors='ignore')
        modification_percentage_summary_df = pd.DataFrame(mod_pcts, columns=colnames)

        # print("modification_percentage_summary_df \n", modification_percentage_summary_df)
        nuc_df_for_plot = df_nuc_pct_all.reset_index().rename(columns={'index':'Nucleotide'})
        nuc_df_for_plot.insert(0, 'Batch', ref_name) #this function was designed for plottin batch... so just add a column in there to make it happy
        mod_df_for_plot = modification_percentage_summary_df.copy()
        mod_df_for_plot.insert(0, 'Batch', ref_name)

        plot_root = _jp(self.OUTPUT_DIRECTORY, '2a.'+ref_plot_name + 'Nucleotide_percentage_quilt')
        plot_2a_input = {
            'nuc_pct_df': nuc_df_for_plot,
            'mod_pct_df': mod_df_for_plot,
            'fig_filename_root': plot_root,
            # 'save_also_png': save_png,
            # 'sgRNA_intervals': sgRNA_intervals,
            # 'sgRNA_names': 'sgRNA_name',
            # 'sgRNA_mismatches': sgRNA_mismatches,
            'quantification_window_idxs': include_idxs_list,
        }
        if n_processes > 1:
            plot_results.append(plot_pool.submit(
                CRISPRessoPlot.plot_nucleotide_quilt,
                **plot_2a_input,
            ))
        else:
            CRISPRessoPlot.plot_nucleotide_quilt(**plot_2a_input)

        self.crispresso2_info['results']['refs']['plot_2a_root'] = os.path.basename(plot_root)
        self.crispresso2_info['results']['refs']['plot_2a_caption'] = "Figure 2a: Nucleotide distribution across amplicon. At each base in the reference amplicon, the percentage of each base as observed in sequencing reads is shown (A = green; C = orange; G = yellow; T = purple). Black bars show the percentage of reads for which that base was deleted. Brown bars between bases show the percentage of reads having an insertion at that position."
        self.crispresso2_info['results']['refs']['plot_2a_data'] = [('Nucleotide frequency table', os.path.basename(nuc_freq_filename))]

        self.crispresso2_info['results']['refs']['plot_2b_roots'] = []
        self.crispresso2_info['results']['refs']['plot_2b_captions'] = []
        self.crispresso2_info['results']['refs']['plot_2b_datas'] = []
        for i in range(self.refs['sequence_length']):
            # get nucleotide columns to print for this sgRNA
            plot_window_size = self.refs['sequence_length']
            # print("plot_window_size ", plot_window_size)
            # print("sgRNA_names ", sgRNA_names)
            sel_cols = [0, 1]
            plot_half_window = max(1, plot_window_size)
            new_sel_cols_start = 2 ##new_sel_cols_start = max(2, cut_point-plot_half_window+1)##modified by Xiaoli 20220809
            new_sel_cols_end = ref_len ##new_sel_cols_end = min(ref_len, cut_point+plot_half_window+1) ##modified by Xiaoli 20220809
            sel_cols.extend(list(range(new_sel_cols_start+2, new_sel_cols_end+2)))
            # print("This is the third: ", str(new_sel_cols_start), str(new_sel_cols_end))
            # get new intervals
            new_sgRNA_intervals = []
            # add annotations for each sgRNA (to be plotted on this sgRNA's plot)
            # for (int_start, int_end) in self.refs['sgRNA_intervals']:
            #    new_sgRNA_intervals += [(int_start - new_sel_cols_start, int_end - new_sel_cols_start)]

            new_include_idx = []
            for x in include_idxs_list:
                new_include_idx += [x - new_sel_cols_start]
            plot_root = _jp(self.OUTPUT_DIRECTORY, '2b.'+ref_plot_name + 'Nucleotide_percentage_quilt_around')

            plot_2b_input = {
                'nuc_pct_df': nuc_df_for_plot.iloc[:, sel_cols],
                'mod_pct_df': mod_df_for_plot.iloc[:, sel_cols],
                'fig_filename_root': plot_root,
                'save_also_png': save_png,
                # 'sgRNA_intervals': new_sgRNA_intervals,
                # 'sgRNA_names': sgRNA_names,
                # 'sgRNA_mismatches': sgRNA_mismatches,
                'quantification_window_idxs': new_include_idx,
            }
            if n_processes > 1:
                plot_results.append(plot_pool.submit(
                    CRISPRessoPlot.plot_nucleotide_quilt, **plot_2b_input,
                ))
            else:
                CRISPRessoPlot.plot_nucleotide_quilt(
                    **plot_2b_input,
                )
            self.crispresso2_info['results']['refs']['plot_2b_roots'].append(os.path.basename(plot_root))
            self.crispresso2_info['results']['refs']['plot_2b_captions'].append('Figure 2b: Nucleotide distribution inside the target segement .')
            self.crispresso2_info['results']['refs']['plot_2b_datas'].append([('Nucleotide frequency in quantification window', os.path.basename(quant_window_nuc_freq_filename))])

    def making_read_length_plot_figure3(self, args):
        # figure3 visualize effective lengths of reads aligning to this amplicon
        if self.counts_total < 1:
            info("No mapped reads to the reference\n")
            sys.exit(1)

        plot_results = []
        n_processes = args.n_processes_for_pooled
        plot_pool = concurrent.futures.ProcessPoolExecutor(n_processes)
        save_png = True

        ref_name = args.amplicon_name
        ref_plot_name = ref_name+'_'
        # print("ref_name ", ref_name)

        hdensity =  self.refs['hdensity']
        hlengths = self.refs['hlengths']
        center_index = self.refs['center_index']

        xmin = min(hlengths)
        xmax = max(hlengths)
        # get 99% cutoff
        plot_histogram_outliers = True
        if not plot_histogram_outliers:
            sum_cutoff = .99 * hdensity.sum()
            sum_so_far = 0
            for idx, val in enumerate(hlengths):
                sum_so_far += hdensity[idx]
                if sum_so_far > sum_cutoff:
                    xmax = val
                    break
            sum_so_far = 0
            for idx, val in enumerate(hlengths[::-1]):
                sum_so_far += hdensity[idx]
                if sum_so_far > sum_cutoff:
                    xmin = val
                    break
        xmin = min(xmin, -15)
        xmax = max(xmax, 15)

        plot_root = _jp(self.OUTPUT_DIRECTORY, '3a.' + ref_plot_name + 'Indel_size_distribution',)
        plot_3a_input = {
            'hdensity': hdensity,
            'hlengths': hlengths,
            'center_index': center_index,
            'n_this_category': self.counts_total,
            'xmin': xmin,
            'xmax': xmax,
            'title': get_plot_title_with_ref_name(
                'Indel size distribution', ref_name,
            ),
            'plot_root': plot_root,
            'save_also_png': save_png,
        }
        if n_processes > 1:
            plot_results.append(plot_pool.submit(
                CRISPRessoPlot.plot_indel_size_distribution,
                **plot_3a_input,
            ))
        else:
            CRISPRessoPlot.plot_indel_size_distribution(**plot_3a_input)
        clipped_string = ""
        if xmax < max(hlengths):
            clipped_string += " (Maximum " + str(int(max(hlengths))) + " not shown)"
        if xmin > min(hlengths):
            clipped_string += " (Minimum " + str(int(min(hlengths))) + " not shown)"
        if clipped_string != "":
            clipped_string = " Note that histograms are clipped to show 99% of the data. To show all data, run using the parameter '--plot_histogram_outliers'. " + clipped_string

        self.crispresso2_info['results']['refs']['plot_3a_root'] = os.path.basename(plot_root)
        self.crispresso2_info['results']['refs']['plot_3a_caption'] = "Figure 3a: Frequency distribution of alleles with indels (blue) and without indels (red)." + clipped_string
        self.crispresso2_info['results']['refs']['plot_3a_data'] = [('Indel histogram', os.path.basename(self.crispresso2_info['results']['refs']['indel_histogram_filename']))]

        #########################################################################################
        # (3b) a graph of frequency of deletions and insertions of various sizes (deletions could be consider as negative numbers and insertions as positive);

        y_values_del = self.refs['y_values_del']
        y_values_ins = self.refs['y_values_ins']
        x_bins_ins  = self.refs['x_bins_ins']
        x_bins_del =  self.refs['x_bins_del']
        x_bins_mut =  self.refs['x_bins_mut']

        xmax_ins = max(x_bins_ins)
        if not plot_histogram_outliers:
            sum_cutoff = 0.99 * hdensity.sum()
            sum_so_far = 0
            for idx, val in enumerate(y_values_ins):
                sum_so_far += x_bins_ins[idx]
                if sum_so_far > sum_cutoff:
                    xmax_ins = val
                    break
        xmax_ins = max(15, xmax_ins)

        clipped_string = ""
        if xmax_ins < max(x_bins_ins):
            clipped_string += " (Insertion maximum " + str(int(max(x_bins_ins))) + " not shown)"
        xmax_del = max(x_bins_del)
        if not plot_histogram_outliers:
            sum_cutoff = .99 * hdensity.sum()
            sum_so_far = 0
            for idx, val in enumerate(y_values_del):
                sum_so_far += x_bins_del[idx]
                if sum_so_far > sum_cutoff:
                    xmax_del = val
                    break
        xmax_del = max(15, xmax_del)

        if xmax_del < max(x_bins_del):
            clipped_string += " (Deletion minimum -" + str(int(max(x_bins_del))) + " not shown)"
        xmax_mut = max(x_bins_mut)
        if not plot_histogram_outliers:
            sum_cutoff = .99 * hdensity.sum()
            sum_so_far = 0
            for idx, val in enumerate(y_values_mut):
                sum_so_far += x_bins_mut[idx]
                if sum_so_far > sum_cutoff:
                    xmax_mut = val
                    break
        xmax_mut = max(15, xmax_mut)
        if xmax_mut < max(x_bins_mut):
            clipped_string += " (Mutation maximum " + str(int(max(x_bins_mut))) + " not shown)"

        plot_root =  _jp(self.OUTPUT_DIRECTORY, '3b.' + ref_plot_name + 'Insertion_deletion_substitutions_size_hist',)
        plot_3b_input = {
            'ref': self.refs,
            'counts_total': self.counts_total,
            'plot_path': plot_root,
            'plot_titles': {
                'ins': get_plot_title_with_ref_name(
                    'Insertions', ref_name,
                ),
                'del': get_plot_title_with_ref_name(
                    'Deletions', ref_name,
                ),
                'mut': get_plot_title_with_ref_name(
                    'Substitutions', ref_name,
                ),
            },
            'xmax_del': xmax_del,
            'xmax_ins': xmax_ins,
            'xmax_mut': xmax_mut,
            'save_also_png': save_png,
        }
        if n_processes > 1:
            plot_results.append(plot_pool.submit(
                CRISPRessoPlot.plot_frequency_deletions_insertions,
                **plot_3b_input,
            ))
        else:
            CRISPRessoPlot.plot_frequency_deletions_insertions(
                **plot_3b_input,
            )

        if clipped_string != "":
            clipped_string = " Note that histograms are clipped to show 99% of the data. To show all data, run using the parameter '--plot_histogram_outliers'. " + clipped_string

        self.crispresso2_info['results']['refs']['plot_3b_root'] = os.path.basename(plot_root)
        self.crispresso2_info['results']['refs']['plot_3b_caption'] = "Figure 3b: Left panel, frequency distribution of sequence modifications that increase read length with respect to the reference amplicon, classified as insertions (positive indel size). Middle panel, frequency distribution of sequence modifications that reduce read length with respect to the reference amplicon, classified as deletions (negative indel size). Right panel, frequency distribution of sequence modifications that do not alter read length with respect to the reference amplicon, which are classified as substitutions (number of substituted positions shown)." + clipped_string
        self.crispresso2_info['results']['refs']['plot_3b_data'] = [('Insertions frequency', self.crispresso2_info['results']['refs']['insertion_histogram_filename']), ('Deletions Frequency', self.crispresso2_info['results']['refs']['deletion_histogram_filename']), ('Substitutions Frequency', self.crispresso2_info['results']['refs']['substitution_histogram_filename'])]

    def making_modification_frequency_plot_figure4(self, args):
        # figure 4 another graph with the frequency that each nucleotide within the amplicon was modified in any way (perhaps would consider insertion as modification of the flanking nucleotides);
        # indels location Plots

        n_this_category_modified = 0
        n_this_category = self.counts_total


        plot_results = []
        n_processes = args.n_processes_for_pooled
        plot_pool = concurrent.futures.ProcessPoolExecutor(n_processes)
        save_png = True

        ref_name = args.amplicon_name
        ref_plot_name = ref_name
        ref_len = len(args.amplicon_seq)

        # cut_points = self.refs['sgRNA_cut_points']
        # plot_cut_points = self.refs['sgRNA_plot_cut_points']
        # sgRNA_intervals = self.refs['sgRNA_intervals']
        # include_idxs_list = self.refs['include_idxs']

        mod_count_filename = _jp(self.OUTPUT_DIRECTORY, 'Modification_count_vectors.txt')

        modifiedName = ref_name + "_MODIFIED"
        if modifiedName in self.class_counts:
            n_this_category_modified = self.class_counts[modifiedName]

        y_max = max(self.all_indelsub_count_vectors) * 1.1
        plot_root = _jp(self.OUTPUT_DIRECTORY, '4a.' + ref_plot_name + 'Combined_insertion_deletion_substitution_locations',)

        plot_4a_input = {
            'all_indelsub_count_vectors': self.all_indelsub_count_vectors,
            'include_idxs_list': self.refs['include_idxs'],
            'cut_points': '',
            'plot_cut_points': '',
            'sgRNA_intervals': '',
            'n_total':self.N_TOTAL,
            'n_this_category': n_this_category,
            'ref_name': ref_name,
            'num_refs': 1,
            'ref_len': ref_len,
            'y_max': y_max,
            'plot_titles': {
                'combined': get_plot_title_with_ref_name(
                    'Combined Insertions/Deletions/Substitutions', ref_name,
                ),
                'main': get_plot_title_with_ref_name(
                    'Mutation position distribution', ref_name,
                ),
            },
            'plot_root': plot_root,
            'save_also_png': save_png,
        }
        if n_processes > 1:
            plot_results.append(plot_pool.submit(
                CRISPRessoPlot.plot_amplicon_modifications, **plot_4a_input,
            ))
        else:
            CRISPRessoPlot.plot_amplicon_modifications(**plot_4a_input)
        self.crispresso2_info['results']['refs']['plot_4a_root'] = os.path.basename(plot_root)
        self.crispresso2_info['results']['refs']['plot_4a_caption'] = "Figure 4a: Combined frequency of any modification across the amplicon. Modifications outside of the quantification window are also shown."
        self.crispresso2_info['results']['refs']['plot_4a_data'] = []

        plot_root = _jp(self.OUTPUT_DIRECTORY,'4b.' + ref_plot_name + 'Insertion_deletion_substitution_locations',)
        plot_4b_input = {
            'include_idxs_list': self.refs['include_idxs'],
            'all_insertion_count_vectors': self.all_insertion_count_vectors,
            'all_deletion_count_vectors': self.all_deletion_count_vectors,
            'all_substitution_count_vectors': self.all_substitution_count_vectors,
            'sgRNA_intervals': '',
            'ref_len': ref_len,
            'ref_name': ref_name,
            'num_refs': 1,
            'n_total': self.N_TOTAL,
            'n_this_category': n_this_category,
            'cut_points': '',
            'plot_cut_points': '',
            'y_max': y_max,
            'plot_title': get_plot_title_with_ref_name(
                'Mutation position distribution', ref_name,
            ),
            'plot_root': plot_root,
            'save_also_png': save_png,
        }
        if n_processes > 1:
            plot_results.append(plot_pool.submit(
                CRISPRessoPlot.plot_modification_frequency,
                **plot_4b_input,
            ))
        else:
            CRISPRessoPlot.plot_modification_frequency(**plot_4b_input)
        self.crispresso2_info['results']['refs']['plot_4b_root'] = os.path.basename(plot_root)
        self.crispresso2_info['results']['refs']['plot_4b_caption'] = "Figure 4b: Frequency of insertions (red), deletions (purple), and substitutions (green) across the entire amplicon, including modifications outside of the quantification window."
        self.crispresso2_info['results']['refs']['plot_4b_data'] = [('Modification frequency', os.path.basename(mod_count_filename))]

        quant_window_mod_count_filename = _jp(self.OUTPUT_DIRECTORY, 'Quantification_window_modification_count_vectors.txt')
        plot_root = _jp(self.OUTPUT_DIRECTORY, '4c.' + ref_plot_name + 'Quantification_window_insertion_deletion_substitution_locations',)
        plot_4c_input = {
            'insertion_count_vectors': self.insertion_count_vectors,
            'deletion_count_vectors': self.deletion_count_vectors,
            'substitution_count_vectors': self.substitution_count_vectors,
            'include_idxs_list': self.refs['include_idxs'],
            'cut_points': '',
            'plot_cut_points': '',
            'sgRNA_intervals': '',
            'ref_len': ref_len,
            'num_refs': 1,
            'n_total': self.N_TOTAL,
            'n_this_category': n_this_category,
            'plot_title': get_plot_title_with_ref_name(
                'Mutation position distribution', ref_name,
            ),
            'ref_name': ref_name,
            'plot_root': plot_root,
            'save_also_png': save_png,
        }
        if n_processes > 1:
            plot_results.append(plot_pool.submit(
                CRISPRessoPlot.plot_quantification_window_locations,
                **plot_4c_input,
            ))
        else:
            CRISPRessoPlot.plot_quantification_window_locations(
                **plot_4c_input,
            )
        self.crispresso2_info['results']['refs']['plot_4c_root'] = os.path.basename(plot_root)
        self.crispresso2_info['results']['refs']['plot_4c_caption'] = "Figure 4c: Frequency of insertions (red), deletions (purple), and substitutions (green) across the entire amplicon, considering only modifications that overlap with the quantification window."
        self.crispresso2_info['results']['refs']['plot_4c_data'] = [('Modification frequency in quantification window', os.path.basename(quant_window_mod_count_filename))]

        # position dependent indels plot
        plot_root = _jp(self.OUTPUT_DIRECTORY, '4d.' + ref_plot_name + 'Position_dependent_average_indel_size',)
        plot_4d_input = {
            'insertion_length_vectors': self.insertion_length_vectors,
            'deletion_length_vectors': self.deletion_length_vectors,
            'cut_points': '',
            'plot_cut_points': '',
            'ref_len': ref_len,
            'plot_titles': {
                'ins': get_plot_title_with_ref_name(
                    'Position dependent insertion size', ref_name,
                ),
                'del': get_plot_title_with_ref_name(
                    'Position dependent deletion size', ref_name,
                ),
            },
            'plot_root': plot_root,
            'save_also_png': save_png,
        }
        if n_processes > 1:
            #print("here we are")
            plot_results.append(plot_pool.submit(
                CRISPRessoPlot.plot_position_dependent_indels,
                **plot_4d_input,
            ))
        else:
            CRISPRessoPlot.plot_position_dependent_indels(
                **plot_4d_input,
            )
        self.crispresso2_info['results']['refs']['plot_4d_root'] = os.path.basename(plot_root)
        self.crispresso2_info['results']['refs']['plot_4d_caption'] = "Figure 4d: Position dependent insertion size(left) and deletion size (right), including only modifications that overlap with the quantification window."
        self.crispresso2_info['results']['refs']['plot_4d_data'] = []

    def making_allele_distribution_plot_figure9(self, args):
        # new plots alleles around cut_sites
        # sgRNA_sequences = self.refs['sgRNA_sequences']
        # sgRNA_cut_points = self.refs['sgRNA_cut_points']
        # sgRNA_plot_cut_points = self.refs['sgRNA_plot_cut_points']
        # sgRNA_intervals = self.refs['sgRNA_intervals']
        # sgRNA_names = 'sgRNA_name'
        # sgRNA_plot_idxs = self.refs['sgRNA_plot_idxs']
        # sgRNA_mismatches = self.refs['sgRNA_mismatches']
        ref_len = len(args.amplicon_seq)
        ref_name = args.amplicon_name
        ref_plot_name = ref_name

        plot_results = []
        n_processes = args.n_processes_for_pooled
        plot_pool = concurrent.futures.ProcessPoolExecutor(n_processes)
        save_png = True

        self.crispresso2_info['results']['refs']['plot_9_roots'] = []
        self.crispresso2_info['results']['refs']['plot_9_captions'] = []
        self.crispresso2_info['results']['refs']['plot_9_datas'] = []
        self.crispresso2_info['results']['refs']['allele_frequency_files'] = []

        plot_window_size = len(args.amplicon_seq)

        df_alleles = pd.DataFrame(self.alleles_list)
        df_alleles['%Reads']=df_alleles['#Reads']/self.N_TOTAL*100
        df_alleles[['n_deleted', 'n_inserted', 'n_mutated']] = df_alleles[['n_deleted', 'n_inserted', 'n_mutated']].astype(int)

        df_alleles.sort_values(by='#Reads', ascending=False, inplace=True)

        # use df_alleles to replace df_alleles_around_cut 
        # df_alleles_around_cut=CRISPRessoShared.get_dataframe_around_cut(df_alleles.loc[df_alleles['Reference_Name'] == ref_name], -3, plot_window_size)
        count_total = self.counts_total

        # write alleles table to file
        # allele_filename = _jp(self.OUTPUT_DIRECTORY, 'Alleles_frequency_table_around_'+sgRNA_label+'.txt')
        allele_filename = _jp(self.OUTPUT_DIRECTORY, 'Alleles_frequency_table.txt')
        df_alleles.to_csv(allele_filename, sep='\t', header=True)
        self.crispresso2_info['results']['refs']['allele_frequency_files'].append(os.path.basename(allele_filename))
        # ref_seq_around_cut=refs[ref_name]['sequence'][cut_point-plot_half_window+1:cut_point+plot_half_window+1]
        ref_seq_around_cut=self.refs['sequence'][0:ref_len]
        fig_filename_root = _jp(self.OUTPUT_DIRECTORY, '9.'+ref_plot_name+'Alleles_frequency_table')
        min_frequency_alleles_around_cut_to_plot = 0.2
        n_good = df_alleles[df_alleles['%Reads']>=min_frequency_alleles_around_cut_to_plot].shape[0]

        if n_good > 0:
            df_to_plot = df_alleles
            df_to_plot = df_alleles.groupby(['Aligned_Sequence', 'Reference_Sequence']).sum().reset_index().set_index('Aligned_Sequence')
            df_to_plot.sort_values(by='%Reads', inplace=True, ascending=False)
            max_rows_alleles_around_cut_to_plot = 50
            plot_9_input = {
                'reference_seq': ref_seq_around_cut,
                'df_alleles': df_to_plot,
                'fig_filename_root': fig_filename_root,
                'MIN_FREQUENCY': min_frequency_alleles_around_cut_to_plot,
                'MAX_N_ROWS': max_rows_alleles_around_cut_to_plot,
                'SAVE_ALSO_PNG': save_png,
                # 'plot_cut_point': plot_cut_point,
                'annotate_wildtype_allele': '',
            }
            '''
            if n_processes > 1:
                plot_results.append(plot_pool.submit(
                    CRISPRessoPlot.plot_alleles_table,
                    **plot_9_input,
                ))
            else:
            '''
            CRISPRessoPlot.plot_alleles_table(
                **plot_9_input,
            )
            self.crispresso2_info['results']['refs']['plot_9_roots'].append(os.path.basename(fig_filename_root))
            self.crispresso2_info['results']['refs']['plot_9_captions'].append("Figure 9: Visualization of the distribution of identified alleles around the cleavage site for the sgRNA_legend. Nucleotides are indicated by unique colors (A = green; C = red; G = yellow; T = purple). Substitutions are shown in bold font. Red rectangles highlight inserted sequences. Horizontal dashed lines indicate deleted sequences. The vertical dashed line indicates the predicted cleavage site.")
            self.crispresso2_info['results']['refs']['plot_9_datas'].append([('Allele frequency table', os.path.basename(allele_filename))])

    def making_nucleotide_substitution_plot_figure10(self, args):
        ref_seq = args.amplicon_seq
        ref_len = len(args.amplicon_seq)
        ref_name = args.amplicon_name
        ref_plot_name = ref_name

        plot_results = []
        n_processes = args.n_processes_for_pooled
        plot_pool = concurrent.futures.ProcessPoolExecutor(n_processes)
        save_png = True

        tot_aln_reads = self.counts_total

        fig_filename_root= _jp(self.OUTPUT_DIRECTORY, '10a.'+ref_plot_name+'Substitution_frequencies_at_each_bp')
        plot_10a_input = {
            'ref_len': ref_len,
            'ref_seq': ref_seq,
            'ref_name': ref_name,
            'ref_count': tot_aln_reads,
            'all_substitution_base_vectors': self.all_substitution_base_vectors,
            'plot_title': get_plot_title_with_ref_name(
                'Substitution frequency', ref_name,
            ),
            'fig_filename_root': fig_filename_root,
            'save_also_png': save_png,
            'quantification_window_idxs': self.refs['include_idxs'],
        }
        if n_processes > 1:
            plot_results.append(plot_pool.submit(
                CRISPRessoPlot.plot_subs_across_ref,
                **plot_10a_input,
            ))
        else:
            CRISPRessoPlot.plot_subs_across_ref(
                **plot_10a_input,
            )
        self.crispresso2_info['results']['refs']['plot_10a_root'] = os.path.basename(fig_filename_root)
        self.crispresso2_info['results']['refs']['plot_10a_caption'] = "Figure 10a: Substitution frequencies across the amplicon."

        df_sub_freq_all, alt_nuc_counts_all = count_alternate_alleles(
            sub_base_vectors = self.all_substitution_base_vectors,
            ref_name = ref_name,
            ref_sequence = ref_seq,
            ref_total_aln_reads = tot_aln_reads
        )


        if 'nuc_freq_filename' in self.crispresso2_info['results']['refs']:
            nuc_freq_filename = self.crispresso2_info['results']['refs']['nuc_freq_filename']
            self.crispresso2_info['results']['refs']['plot_10a_data'] = [('Nucleotide frequencies', os.path.basename(nuc_freq_filename))]

        # plot all substitution rates in entire region
        fig_filename_root = _jp(self.OUTPUT_DIRECTORY, '10b.'+ref_plot_name+'Substitution_frequency_barplot')
        nuc_freq_filename = _jp(self.OUTPUT_DIRECTORY, ref_plot_name + 'Nucleotide_frequency_table.txt')
        plot_10b_input = {
            'alt_nuc_counts': alt_nuc_counts_all,
            'plot_title': get_plot_title_with_ref_name(
                'Substitution frequency\nin entire amplicon', ref_name,
            ),
            'fig_filename_root': fig_filename_root,
            'save_also_png': save_png
        }
        if n_processes > 1:
            plot_results.append(plot_pool.submit(
                CRISPRessoPlot.plot_sub_freqs,
                **plot_10b_input,
            ))
        else:
            CRISPRessoPlot.plot_sub_freqs(
                **plot_10b_input,
            )
        self.crispresso2_info['results']['refs']['plot_10b_root'] = os.path.basename(fig_filename_root)
        self.crispresso2_info['results']['refs']['plot_10b_caption'] = "Figure 10b: Substitution frequencies across the amplicon."
        self.crispresso2_info['results']['refs']['plot_10b_data'] = [('Nucleotide frequencies', os.path.basename(nuc_freq_filename))]

        # plot all substitution rates in quantification_window
        fig_filename_root = _jp(self.OUTPUT_DIRECTORY, '10c.'+ref_plot_name+'Substitution_frequency_barplot_in_quantification_window')

        include_idxs_list = self.refs['include_idxs']
        quantification_window_ref_seq = [list(ref_seq)[x] for x in include_idxs_list]

        df_sub_freq, alt_nuc_counts = count_alternate_alleles(
            sub_base_vectors = self.substitution_base_vectors,
            ref_name = ref_name,
            ref_sequence = quantification_window_ref_seq,
            ref_total_aln_reads = tot_aln_reads
        )
        plot_10c_input = {
            'alt_nuc_counts': alt_nuc_counts,
            'plot_title': get_plot_title_with_ref_name('Substitution frequency\nin quantification window', ref_name),
            'fig_filename_root': fig_filename_root,
            'save_also_png': save_png
        }
        if n_processes > 1:
            plot_results.append(plot_pool.submit(
                CRISPRessoPlot.plot_sub_freqs,
                **plot_10c_input,
            ))
        else:
            CRISPRessoPlot.plot_sub_freqs(
                **plot_10c_input,
            )

        quant_window_sub_freq_filename =_jp(self.OUTPUT_DIRECTORY, ref_plot_name + 'Quantification_window_substitution_frequency_table.txt')
        self.crispresso2_info['results']['refs']['plot_10c_root'] = os.path.basename(fig_filename_root)
        self.crispresso2_info['results']['refs']['plot_10c_caption'] = "Figure 10c: Substitution frequencies in the quantification window"
        self.crispresso2_info['results']['refs']['plot_10c_data'] = [('Nucleotide frequencies in quantification window', os.path.basename(quant_window_sub_freq_filename))]

        # guide-specific base editor plots figure 10d, e, f and g
        self.crispresso2_info['results']['refs']['plot_10d_roots'] = []
        self.crispresso2_info['results']['refs']['plot_10d_captions'] = []
        self.crispresso2_info['results']['refs']['plot_10d_datas'] = []

        self.crispresso2_info['results']['refs']['plot_10e_roots'] = []
        self.crispresso2_info['results']['refs']['plot_10e_captions'] = []
        self.crispresso2_info['results']['refs']['plot_10e_datas'] = []

        self.crispresso2_info['results']['refs']['plot_10f_roots'] = []
        self.crispresso2_info['results']['refs']['plot_10f_captions'] = []
        self.crispresso2_info['results']['refs']['plot_10f_datas'] = []

        self.crispresso2_info['results']['refs']['plot_10g_roots'] = []
        self.crispresso2_info['results']['refs']['plot_10g_captions'] = []
        self.crispresso2_info['results']['refs']['plot_10g_datas'] = []

        plot_half_window = 70
        # cut_point = self.refs['sgRNA_cut_points'][0]
        # print("cut_points ", self.refs['sgRNA_cut_points'])

        # left_end = max(1, cut_point-plot_half_window+1)
        # right_end = min(cut_point+plot_half_window+1, ref_len)
        # ref_seq_around_cut=self.refs['sequence'][left_end:right_end]
        ref_seq_around_cut=self.refs['sequence']

        df_nuc_freq_all = pd.DataFrame([self.all_base_count_vectors[ref_name+"_A"], self.all_base_count_vectors[ref_name+"_C"], self.all_base_count_vectors[ref_name+"_G"], self.all_base_count_vectors[ref_name+"_T"], self.all_base_count_vectors[ref_name+"_N"], self.all_base_count_vectors[ref_name+'_-']])
        df_nuc_freq_all.index = ['A', 'C', 'G', 'T', 'N', '-']
        df_nuc_freq_all.columns = list(ref_seq)

        plot_idxs = self.refs['include_idxs']
        df_nuc_pct_all = df_nuc_freq_all.divide(tot_aln_reads)

        plot_ref_seq = ref_seq_around_cut
        plot_nuc_pcts = df_nuc_pct_all.iloc[:, plot_idxs]
        plot_nuc_freqs = df_nuc_freq_all.iloc[:, plot_idxs]

        # get computation window in plotted region
        is_window = np.zeros(ref_len)
        for ind in include_idxs_list:
            is_window[ind] = 1
        plot_is_window = np.zeros(len(plot_idxs)) #binary array whether base should be plotted
        plot_quant_window_idxs = []
        for ind, loc in enumerate(plot_idxs):
            plot_is_window[ind] = is_window[loc]
            if is_window[loc]:
                plot_quant_window_idxs.append(ind-2)

        # sgRNA = self.refs['sgRNA_orig_sequences'][0]
        # sgRNA_label = "sgRNA_"+ sgRNA # for file names
        # sgRNA_legend = "sgRNA " + sgRNA # for legends
        from_nuc_indices = [pos for pos, char in enumerate(list(plot_nuc_pcts.columns.values)) if char == args.conversion_nuc_from]
        just_sel_nuc_pcts = plot_nuc_pcts.iloc[:, from_nuc_indices].copy() #only nucleotides targeted by base editing
        # just_sel_nuc_pcts.columns = [char + str(pos+1) for pos,char in enumerate(list(just_sel_nuc_pcts.columns.values))]
        just_sel_nuc_pcts.columns = [args.conversion_nuc_from + str(pos+1) for pos in from_nuc_indices]
        just_sel_nuc_freqs = plot_nuc_freqs.iloc[:, from_nuc_indices].copy()
        just_sel_nuc_freqs.columns = [args.conversion_nuc_from + str(pos+1) for pos in from_nuc_indices]

        quant_window_sel_nuc_pct_filename = _jp(self.OUTPUT_DIRECTORY, ref_plot_name + 'Selected_nucleotide_percentage_table_around_'+sgRNA_label+'.txt')
        just_sel_nuc_pcts.to_csv(quant_window_sel_nuc_pct_filename, sep='\t', header=True, index=True)
        # not storing the name because it is unique to this sgRNA
        # crispresso2_info['quant_window_sel_nuc_pct_filename'] = os.path.basename(quant_window_sel_nuc_pct_filename)


        quant_window_sel_nuc_freq_filename = _jp(self.OUTPUT_DIRECTORY, ref_plot_name + 'Selected_nucleotide_frequency_table_around_'+sgRNA_label+'.txt')
        just_sel_nuc_freqs.to_csv(quant_window_sel_nuc_freq_filename, sep='\t', header=True, index=True)
        # not storing the name because it is unique to this sgRNA
        # crispresso2_info['quant_window_sel_nuc_freq_filename'] = os.path.basename(quant_window_sel_nuc_freq_filename)
 
        # print table showing all nuc frequencies (sum to total alleles) (in entire region)

        #                CRISPRessoPlot.plot_nuc_freqs(
        #                    df_nuc_freq = df_nuc_freq,
        #                    tot_aln_reads = tot_aln_reads,
        #                    plot_title = get_plot_title_with_ref_name('Nucleotide Frequencies',ref_name),
        #                    fig_filename_root = _jp('14a.'+ref_name+'.nucleotide_frequency'),
        #                    save_also_png = save_png
        #                    )

        fig_filename_root = _jp(self.OUTPUT_DIRECTORY,
            '10d.' + ref_plot_name + 'Log2_nucleotide_frequency_around_' + sgRNA_label,
        )
        plot_10d_input = {
            'df_nuc_freq': plot_nuc_freqs,
            'tot_aln_reads': tot_aln_reads,
            'plot_title': get_plot_title_with_ref_name(
                'Log2 Nucleotide Frequencies Around the ' + sgRNA_legend,
                ref_name,
            ),
            'fig_filename_root': fig_filename_root,
            'save_also_png': save_png,
            'quantification_window_idxs': plot_quant_window_idxs,
        }
        if n_processes > 1:
            plot_results.append(plot_pool.submit(
                CRISPRessoPlot.plot_log_nuc_freqs,
                **plot_10d_input,
            ))
        else:
            CRISPRessoPlot.plot_log_nuc_freqs(
                **plot_10d_input,
            )
        self.crispresso2_info['results']['refs']['plot_10d_roots'].append(os.path.basename(fig_filename_root))
        self.crispresso2_info['results']['refs']['plot_10d_captions'].append("Figure 10d: Log2 nucleotide frequencies for each position in the plotting window around the " + sgRNA_legend + ". The quantification window is outlined by the dotted box.")
        self.crispresso2_info['results']['refs']['plot_10d_datas'].append([])

        fig_filename_root = _jp(self.OUTPUT_DIRECTORY,
            '10e.' + ref_plot_name + 'Selected_conversion_at_' + args.conversion_nuc_from + 's_around_' + sgRNA_label,
        )
        plot_10e_input = {
            'df_subs': plot_nuc_pcts,
            'ref_name': ref_name,
            'ref_sequence': plot_ref_seq,
            'plot_title': get_plot_title_with_ref_name(
                'Substitution Frequencies at ' + args.conversion_nuc_from + 's around the ' + sgRNA_legend,
                ref_name,
            ),
            'conversion_nuc_from': args.conversion_nuc_from,
            'fig_filename_root': fig_filename_root,
            'save_also_png': save_png,
        }
        if n_processes > 1:
            plot_results.append(plot_pool.submit(
                CRISPRessoPlot.plot_conversion_at_sel_nucs,
                **plot_10e_input,
            ))
        else:
            CRISPRessoPlot.plot_conversion_at_sel_nucs(
                **plot_10e_input,
            )
        self.crispresso2_info['results']['refs']['plot_10e_roots'].append(os.path.basename(fig_filename_root))
        self.crispresso2_info['results']['refs']['plot_10e_captions'].append("Figure 10e: Proportion of each base at each nucleotide targeted by base editors in the plotting window around the " + sgRNA_legend + ". The number of each target base is annotated on the reference sequence at the bottom of the plot.")
        self.crispresso2_info['results']['refs']['plot_10e_datas'].append([('Nucleotide frequencies at ' + args.conversion_nuc_from + 's', os.path.basename(quant_window_sel_nuc_freq_filename))])

        fig_filename_root = _jp(self.OUTPUT_DIRECTORY,'10f.'+ref_plot_name+'Selected_conversion_no_ref_at_'+args.conversion_nuc_from+'s_around_'+sgRNA_label)
        plot_10f_input = {
            'df_subs': plot_nuc_pcts,
            'ref_name': ref_name,
            'ref_sequence': plot_ref_seq,
            'plot_title': get_plot_title_with_ref_name('Substitution Frequencies at '+args.conversion_nuc_from+'s around the ' + sgRNA_legend, ref_name),
            'conversion_nuc_from': args.conversion_nuc_from,
            'fig_filename_root': fig_filename_root,
            'save_also_png': save_png,
        }
        if n_processes > 1:
            plot_results.append(plot_pool.submit(
                CRISPRessoPlot.plot_conversion_at_sel_nucs_not_include_ref,
                **plot_10f_input,
            ))
        else:
            CRISPRessoPlot.plot_conversion_at_sel_nucs_not_include_ref(
                **plot_10f_input,
            )
        self.crispresso2_info['results']['refs']['plot_10f_roots'].append(os.path.basename(fig_filename_root))
        self.crispresso2_info['results']['refs']['plot_10f_captions'].append("Figure 10f: Non-reference base proportions. For target nucleotides in the plotting window, this plot shows the proportion of non-reference (non-"+args.conversion_nuc_from + ") bases as a percentage of all non-reference sequences. The number of each target base is annotated on the reference sequence at the bottom of the plot.")
        self.crispresso2_info['results']['refs']['plot_10f_datas'].append([('Nucleotide frequencies at ' + args.conversion_nuc_from + 's', os.path.basename(quant_window_sel_nuc_freq_filename))])

        fig_filename_root = _jp(self.OUTPUT_DIRECTORY, '10g.'+ref_plot_name+'Selected_conversion_no_ref_scaled_at_'+args.conversion_nuc_from+'s_around_'+sgRNA_label)
        plot_10g_input = {
            'df_subs': plot_nuc_pcts,
            'ref_name': ref_name,
            'ref_sequence': plot_ref_seq,
            'plot_title': get_plot_title_with_ref_name(
                'Substitution Frequencies at ' + args.conversion_nuc_from + 's around the ' + sgRNA_legend,
                ref_name,
            ),
            'conversion_nuc_from': args.conversion_nuc_from,
            'fig_filename_root': fig_filename_root,
            'save_also_png': save_png
        }
        if n_processes > 1:
            plot_results.append(plot_pool.submit(
                CRISPRessoPlot.plot_conversion_at_sel_nucs_not_include_ref_scaled,
                **plot_10g_input,
            ))
        else:
            CRISPRessoPlot.plot_conversion_at_sel_nucs_not_include_ref_scaled(
                **plot_10g_input,
            )
        self.crispresso2_info['results']['refs']['plot_10g_roots'].append(os.path.basename(fig_filename_root))
        self.crispresso2_info['results']['refs']['plot_10g_captions'].append("Figure 10g: Non-reference base counts. For target nucleotides in the plotting window, this plot shows the number of non-reference (non-" + args.conversion_nuc_from + ") bases. The number of each target base is annotated on the reference sequence at the bottom of the plot.")
        self.crispresso2_info['results']['refs']['plot_10g_datas'].append([('Nucleotide frequencies at ' + args.conversion_nuc_from +'s', os.path.basename(quant_window_sel_nuc_freq_filename))])
        # END GUIDE SPECIFIC PLOTS


if __name__ == '__main__':
    # start time
    start_time = time.time()

    # initialize parser
    parser = argparse.ArgumentParser(description="Compare the edited reads to the given amplicon sequences and find out variations.",formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # define command line arguments
    parser.add_argument('-fq', '--merged_fastq', type=str,  help='Fastq file containing merged reads which were assigned to corresponding amplicons', default='', required=True)
    parser.add_argument('--amplicon_name', type=str,  help='The name of amplicons', default='', required=True)
    parser.add_argument('-o','--output_dir', type=str,  help='The directory of output files.', default='', required=True)


    parser.add_argument('--amplicon_seq', type=str,  help='Amplicon Sequence', required=True)
    # parser.add_argument('-et', '--edited_type', type=str,  help='Type of edition (ABE, CBE, STEME)', required=False)
    parser.add_argument('-p', '--n_processes_for_pooled', type = int, default=1,  help='Multiple processes number')

    parser.add_argument('-wc', '--quantification_window_center', type=int, help="Center of quantification window to use within respect to the 3' end of the provided sgRNA sequence. Remember that the sgRNA sequence must be entered without the PAM. For cleaving nucleases, this is the predicted cleavage position. The default is -3 and is suitable for the Cas9 system. For alternate nucleases, other cleavage offsets may be appropriate, for example, if using Cpf1 this parameter would be set to 1. For base editors, this could be set to -17 to only include mutations near the 5' end of the sgRNA. Multiple quantification window centers (corresponding to each guide specified by --guide_seq) can be specified with a comma-separated list.", default = -3)

    parser.add_argument('--quantification_window_size', type=int, help="Defines the size (in bp) of the quantification window extending from the position specified by the \"--quantification_window_center\" parameter in relation to the provided guide RNA sequence(s) (--sgRNA). Mutations within this number of bp from the quantification window center are used in classifying reads as modified or unmodified. A value of 0 disables this window and indels in the entire amplicon are considered. Default is 20, 20bp on each side of the cleavage position for a total length of 40bp. Multiple quantification window sizes (corresponding to each guide specified by --guide_seq) can be specified with a comma-separated list", default = 20)

    parser.add_argument('--left_primer_seq', type = str, help='Sequences of left primer of the amplicon sequence')
    parser.add_argument('--right_primer_seq', type = str, help='Sequences of right primer of the amplicon sequence')

    parser.add_argument('--min_average_read_quality', type=int, help='Minimum average quality score (phred33) to keep a read', default=30)

    parser.add_argument('--min_bp_quality_or_N', type=int, help='Bases with a quality score (phred33) less than this value will be set to "N"', default=20)

    # read arguments from command line
    args = parser.parse_args()
    print(args)

    if args.merged_fastq:
        check_file(args.merged_fastq)
    else:
        arg_parser.print_help()
        raise CRISPRessoShared.BadParameterException('Please provide input data for analysis e.g. using the --fastq_r1 parameter.')

    if args.amplicon_seq is None:
        arg_parser.print_help()
        raise Exception('Please provide an amplicon sequence for analysis using the --amplicon_seq parameter.')

    # normalize name and remove not allowed characters
    if not args.amplicon_name:
        arg_parser.print_help()
        raise Exception('Please provide an amplicon name for analysis using the --amplicon_seq parameter.')

    my_edition_stat_each_amplicon = Edition_stat_each_amplicon(args)##initializing a new class
    my_edition_stat_each_amplicon.refObj_for_alignment(args)
    my_edition_stat_each_amplicon.reads_amplicon_alignment(args)
    my_edition_stat_each_amplicon.analyze_alignments(args)

    my_edition_stat_each_amplicon.making_alignment_plot_figure1(args)
    my_edition_stat_each_amplicon.making_nuc_pct_plot_figure2(args)
    # my_edition_stat_each_amplicon.making_read_length_plot_figure3(args)
    # my_edition_stat_each_amplicon.making_modification_frequency_plot_figure4(args)
    # my_edition_stat_each_amplicon.making_allele_distribution_plot_figure9(args)
    # my_edition_stat_each_amplicon.making_nucleotide_substitution_plot_figure10(args)

    # get end time
    end_time = time.time()
    # calculating running time
    runTime = end_time - start_time
    # output running time
    print(f"Elapsed Time is: {runTime} Seconds\n")


