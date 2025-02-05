'''
CRISPResso2 - Kendell Clement and Luca Pinello 2020
Software pipeline for the analysis of genome editing outcomes from deep sequencing data
(c) 2020 The General Hospital Corporation. All Rights Reserved.
'''

import argparse
import datetime
import errno
import gzip
import json
import numpy as np
import os
import pandas as pd
import re
import string
import shutil
import signal
import subprocess as sb
import unicodedata

__version__ = "2.2.9"

def find_wrong_nt(sequence):
    return list(set(sequence.upper()).difference({'A', 'T', 'C', 'G', 'N'}))

def reverse_complement(seq):
    nt_complement=dict({'A':'T','C':'G','G':'C','T':'A','N':'N','_':'_','-':'-'})
    return "".join([nt_complement[c] for c in seq.upper()[-1::-1]])

def set_guide_array(vals, guides, property_name):
    """
    creates an int array of vals of the same length as guide, filling in missing values with the first item in vals
    vals: comma-separated string of values
    guides: list of guides
    property_name: property name, useful for creating warning to user

    returns: list of int vals
    """
    #if only one is given, set it for all guides
    vals_array = str(vals).split(",")
    ret_array = [int(vals_array[0])]*len(guides)
    #len(vals_array) is always one -- don't freak out if it's longer than 0 guides
    if len(vals_array) > 1 and len(vals_array) > len(guides):
        raise Exception("More %s values were given than guides. Guides: %d %ss: %d"%(property_name, len(guides), property_name, len(vals_array)))

    if len(guides) == 0:
        return []

    for idx, val in enumerate(vals_array):
        if val != '':
            ret_array[idx] = int(val)
    return ret_array

#def get_amplicon_info_for_guides(ref_seq,guides,guide_mismatches,guide_names,quantification_window_centers,quantification_window_sizes,quantification_window_coordinates,exclude_bp_from_left,exclude_bp_from_right,plot_window_size,guide_plot_cut_points,discard_guide_positions_overhanging_amplicon_edge=False):
def get_amplicon_info_for_guides(ref_seq,guides,guide_mismatches,quantification_window_centers,quantification_window_sizes,exclude_bp_from_left,exclude_bp_from_right,guide_plot_cut_points):
    """
    gets cut site and other info for a reference sequence and a given list of guides

    input:
    ref_seq : reference sequence
    guides : a list of guide sequences
    guide_mismatches : a list of positions where a guide may have mismatches (for flexiguides)
    guide_names : a list of names for each guide
    quantification_window_centers : a list of positions where quantification is centered for each guide
    quantification_window_sizes : a list of lengths of quantification windows extending from quantification_window_center for each guide
    quantification_window_coordinates: if given, these override quantification_window_center and quantification_window_size for setting quantification window. These are specific for this amplicon.
    exclude_bp_from_left : these bp are excluded from the quantification window
    exclude_bp_from_right : these bp are excluded from the quantification window
    #plot_window_size : length of window extending from quantification_window_center to plot
    guide_plot_cut_points : whether or not to add cut point to plot (prime editing flaps don't have cut points)
    #discard_guide_positions_overhanging_amplicon_edge : if True, for guides that align to multiple positions, guide positions will be discarded if plotting around those regions would included bp that extend beyond the end of the amplicon.

    returns:
    this_sgRNA_sequences : list of sgRNAs that are in this amplicon
    this_sgRNA_intervals : indices of each guide
    this_sgRNA_cut_points : cut points for each guide (defined by quantification_window_center)
    this_sgRNA_plot_cut_points : whether or not a cut point is plotted
    this_sgRNA_plot_idxs : list of indices to be plotted for each sgRNA
    this_sgRNA_mismatches: list of mismatches between the guide and the amplicon
    this_sgRNA_names : list of names for each sgRNA (to disambiguate in case a sequence aligns to multiple positions)
    this_include_idxs : list of indices to be included in quantification
    this_exclude_idxs : list of indices to be excluded from quantification
    """
    ref_seq_length = len(ref_seq)

    this_sgRNA_sequences = []
    this_sgRNA_intervals = []
    this_sgRNA_cut_points = []
    this_sgRNA_plot_cut_points = []
    this_sgRNA_plot_idxs=[]
    this_sgRNA_mismatches = []
    this_sgRNA_names = []
    this_include_idxs=[]
    this_exclude_idxs=[]

    if exclude_bp_from_left:
       this_exclude_idxs+=range(exclude_bp_from_left)

    if exclude_bp_from_right:
       this_exclude_idxs+=range(ref_seq_length)[-exclude_bp_from_right:]

    #window_around_cut=max(1, plot_window_size)

    seen_cut_points = {} #keep track of cut points in case 2 gudes cut at same position (so they can get different names)
    #seen_guide_names = {} #keep track of guide names (so we don't assign a guide the same name as another guide)
    for guide_idx, current_guide_seq in enumerate(guides):
        if current_guide_seq == '':
            continue
        offset_fw=quantification_window_centers[guide_idx]+len(current_guide_seq)-1
        offset_rc=(-quantification_window_centers[guide_idx])-1

        #.. run once with findall to get number of matches
        fw_matches = re.findall(current_guide_seq, ref_seq, flags=re.IGNORECASE)
        rv_matches = re.findall(reverse_complement(current_guide_seq), ref_seq, flags=re.IGNORECASE)
        match_count = len(fw_matches) + len(rv_matches)

        #and now create the iter which will keep track of the locations of matches
        fw_matches = re.finditer(current_guide_seq, ref_seq, flags=re.IGNORECASE)
        rv_matches = re.finditer(reverse_complement(current_guide_seq), ref_seq, flags=re.IGNORECASE)
        '''
        for fw in fw_matches:
            print(fw.start())
        for rv in rv_matches:
            print(rv.start())
        '''
        #for every match, append:
        # this_sgRNA_cut_points, this_sgRNA_intervals,this_sgRNA_mismatches,this_sgRNA_names,this_sgRNA_sequences,include_idxs
        for m in fw_matches:
            cut_p = m.start() + offset_fw ##convert the cut positon coordination from guide RNA coordination to amplicon seq coordination
            this_sgRNA_cut_points.append(cut_p)
            this_sgRNA_plot_cut_points.append(guide_plot_cut_points[guide_idx])
            this_sgRNA_intervals.append((m.start(), m.start() + len(current_guide_seq)-1))
            this_sgRNA_mismatches.append(guide_mismatches[guide_idx])

            if quantification_window_sizes[guide_idx] > 0:
                st=max(0, cut_p-quantification_window_sizes[guide_idx]+1)
                en=min(ref_seq_length-1, cut_p+quantification_window_sizes[guide_idx]+1)
                this_include_idxs.extend(range(st, en))

            '''
            this_sgRNA_name = guide_names[guide_idx]
            if match_count == 1:
                this_sgRNA_names.append(this_sgRNA_name)
            else:
                if this_sgRNA_name == "":
                    this_sgRNA_name = current_guide_seq.upper()
                this_potential_name = this_sgRNA_name + "_" + str(m.start())
                curr_guide_idx = 1
                while this_potential_name in seen_guide_names:
                    curr_guide_idx += 1
                    this_potential_name = this_sgRNA_name + "_" + str(m.start()) + "_" + curr_guide_idx
                this_sgRNA_names.append(this_potential_name)
            '''
            this_sgRNA_sequences.append(current_guide_seq.upper())

        for m in rv_matches:
            cut_p = m.start() + offset_rc #cut position
            this_sgRNA_cut_points.append(cut_p)
            this_sgRNA_plot_cut_points.append(guide_plot_cut_points[guide_idx])
            this_sgRNA_intervals.append((m.start(), m.start() + len(current_guide_seq)-1))
            this_sgRNA_mismatches.append([len(current_guide_seq)-(x+1) for x in guide_mismatches[guide_idx]])

            if quantification_window_sizes[guide_idx] > 0:
                st=max(0, cut_p-quantification_window_sizes[guide_idx]+1)
                en=min(ref_seq_length-1, cut_p+quantification_window_sizes[guide_idx]+1)
                this_include_idxs.extend(range(st, en))

            #this_sgRNA_name = guide_names[guide_idx]
            '''
            if match_count == 1:
                this_sgRNA_names.append(this_sgRNA_name)
            else:
                if this_sgRNA_name == "":
                    this_sgRNA_name = current_guide_seq.upper()
                this_potential_name = this_sgRNA_name + "_" + str(m.start())
                curr_guide_idx = 1
                while this_potential_name in seen_guide_names:
                    curr_guide_idx += 1
                    this_potential_name = this_sgRNA_name + "_" + str(m.start()) + "_" + curr_guide_idx
                this_sgRNA_names.append(this_potential_name)
            '''
            this_sgRNA_sequences.append(current_guide_seq.upper())


    #create mask of positions in which to include/exclude indels for the quantification window
    #first, if exact coordinates have been given, set those
    given_include_idxs = []
    '''
    if quantification_window_coordinates is not None:
        coordinate_include_idxs = []
        theseCoords = str(quantification_window_coordinates).split("_")
        for coord in theseCoords:
            coordRE = re.match(r'^(\d+)-(\d+)$', coord)
            if coordRE:
                start = int(coordRE.group(1))
                end = int(coordRE.group(2)) + 1
                if end > ref_seq_length:
                    raise NTException("End coordinate " + str(end) + " for '" + str(coord) + "' in '" + str(theseCoords) + "' is longer than the sequence length ("+str(ref_seq_length)+")")
                coordinate_include_idxs.extend(range(start, end))
            else:
                raise NTException("Cannot parse analysis window coordinate '" + str(coord) + "' in '" + str(theseCoords) + "'. Coordinates must be given in the form start-end e.g. 5-10 . Please check the --analysis_window_coordinate parameter.")
        given_include_idxs = coordinate_include_idxs
        this_include_idxs = coordinate_include_idxs
    '''
    #otherwise, take quantification centers + windows calculated for each guide above
    #elif this_sgRNA_cut_points and len(this_include_idxs) > 0:
    if this_sgRNA_cut_points and len(this_include_idxs) > 0:
        given_include_idxs = this_include_idxs
    else:
        this_include_idxs=range(ref_seq_length)

    #flatten the arrays to avoid errors with old numpy library
    this_include_idxs=np.ravel(this_include_idxs)
    this_exclude_idxs=np.ravel(this_exclude_idxs)
    given_include_idxs=np.ravel(given_include_idxs)
    pre_exclude_include_idxs = this_include_idxs.copy()
    this_include_idxs=set(np.setdiff1d(this_include_idxs, this_exclude_idxs))##remove the idxs present in this_exclude_idxs from the this_include_idxs
    '''
    print("this_include_idxs ", str(this_include_idxs))
    print("this_exclude_idxs ", str(this_exclude_idxs))
    print("given_include_idxs ", str(given_include_idxs))
    print("pre_exclude_include_idxs ", str(pre_exclude_include_idxs))
    exit(0)
    '''
    if len(np.setdiff1d(given_include_idxs, list(this_include_idxs))) > 0:
        raise BadParameterException('The quantification window has been partially exluded by the --exclude_bp_from_left or --exclude_bp_from_right parameters. Given: ' + str(given_include_idxs) + ' Pre: ' + str(pre_exclude_include_idxs) + ' Post: ' + str(this_include_idxs))
    if len(this_include_idxs) == 0:
        if len(pre_exclude_include_idxs) > 0:
            raise BadParameterException('The quantification window around the sgRNA is excluded. Please decrease the exclude_bp_from_right and exclude_bp_from_left parameters.')
        else:
            raise BadParameterException('The entire sequence has been excluded. Please enter a longer amplicon, or decrease the exclude_bp_from_right and exclude_bp_from_left parameters')

    '''
    if this_sgRNA_cut_points and plot_window_size>0:
        for cut_p in this_sgRNA_cut_points:
            if cut_p - window_around_cut + 1 < 0:
                raise BadParameterException('Offset around cut would extend to the left of the amplicon. Please decrease plot_window_size parameter. Cut point: ' + str(cut_p) + ' window: ' + str(window_around_cut) + ' reference: ' + str(ref_seq_length))
            if cut_p + window_around_cut > ref_seq_length-1:
                raise BadParameterException('Offset around cut would be greater than reference sequence length. Please decrease plot_window_size parameter. Cut point: ' + str(cut_p) + ' window: ' + str(window_around_cut) + ' reference: ' + str(ref_seq_length))
            st=max(0, cut_p-window_around_cut+1) 
            en=min(ref_seq_length-1, cut_p+window_around_cut+1) 
            this_sgRNA_plot_idxs.append(sorted(list(range(st, en))))
    else:
       this_sgRNA_plot_idxs.append(range(ref_seq_length))
    '''
    st=  exclude_bp_from_left 
    en = ref_seq_length - exclude_bp_from_right
    this_sgRNA_plot_idxs.append(sorted(list(range(st, en))))

    this_include_idxs = np.sort(list(this_include_idxs))
    this_exclude_idxs = np.sort(list(this_exclude_idxs))

    #return this_sgRNA_sequences, this_sgRNA_intervals, this_sgRNA_cut_points, this_sgRNA_plot_cut_points, this_sgRNA_plot_idxs, this_sgRNA_mismatches, this_sgRNA_names, this_include_idxs, this_exclude_idxs
    return this_sgRNA_sequences, this_sgRNA_intervals, this_sgRNA_cut_points, this_sgRNA_plot_cut_points, this_sgRNA_plot_idxs, this_sgRNA_mismatches, this_include_idxs, this_exclude_idxs

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


def slugify(value): #adapted from the Django project

    value = unicodedata.normalize('NFKD', value).encode('ascii', 'ignore')
    value = re.sub(rb'[^\w\s-]', b'_', value).strip()
    value = re.sub(rb'[-\s]+', b'-', value)

    return value.decode('utf-8')


def find_indels_substitutions(read_seq_al, ref_seq_al, _include_indx):
    #ref_positions holds the indices for which positions map back to the original reference
    # for example,
    #     1 2 3 4 5 6 7 8
    # ref A A T T G G C C
    #
    # and for a given alignment
    # ref A A T T - G G C C
    # aln A - T T T G G C C
    #     1 2 3 4-4 5 6 7 8 <ref positions. Note that the negative values/indices represent places that don't map back to the original reference
    ref_positions=[]
    all_substitution_positions=[]##reference position where read has substitution
    substitution_positions=[] ##reference position in _include_indx where read has substitution
    all_substitution_values=[] ## The subtituted nucleotide in read sequence
    substitution_values=[] ##The substituted nucleotide in read sequence 

    all_deletion_positions = []
    deletion_positions = []
    deletion_coordinates = []
    deletion_sizes = []
    start_deletion = -1  # the -1 value indicates that there currently isn't a deletion

    all_insertion_positions = []
    all_insertion_left_positions = []
    insertion_positions = []
    insertion_coordinates = []
    insertion_sizes = []
    start_insertion = -1  # the -1 value indicates that there currently isn't an insertion

    seq_len = len(ref_seq_al)
    include_indx_set = set(_include_indx)
    nucSet = set(['A', 'T', 'C', 'G', 'N'])
    idx = 0
    current_insertion_size = 0
    for idx_c, c in enumerate(ref_seq_al):
        if c != '-':##no indel in reference
            ref_positions.append(idx)
            if ref_seq_al[idx_c]!=read_seq_al[idx_c] and read_seq_al[idx_c] != '-' and read_seq_al[idx_c] != 'N':
                all_substitution_positions.append(idx)
                all_substitution_values.append(read_seq_al[idx_c])
                if idx in _include_indx:
                    substitution_positions.append(idx)
                    substitution_values.append(read_seq_al[idx_c])
            if start_insertion != -1:  # this is the end of an insertion
                all_insertion_left_positions.append(start_insertion)
                all_insertion_positions.append(start_insertion)
                all_insertion_positions.append(idx)
                if start_insertion in include_indx_set and idx in include_indx_set:
                    insertion_coordinates.append((start_insertion, idx))
                    insertion_positions.append(start_insertion)
                    insertion_positions.append(idx)
                    insertion_sizes.append(current_insertion_size)
                start_insertion = -1
            current_insertion_size = 0
            idx += 1
        else:  # the current ref position is -
            if idx == 0:
                ref_positions.append(-1)
            else:
                ref_positions.append(-idx)
            if idx > 0 and start_insertion == -1:  # this is the first index of an insertion
                start_insertion = idx - 1
            current_insertion_size += 1

        if read_seq_al[idx_c] == '-' and start_deletion == -1:  # this is the first part of a deletion
            if idx_c - 1 > 0:
                start_deletion = ref_positions[idx_c]
            else:
                start_deletion = 0
        elif read_seq_al[idx_c] != '-' and start_deletion != -1:  # this is the end of a deletion
            end_deletion = ref_positions[idx_c]
            all_deletion_positions.extend(range(start_deletion, end_deletion))
            if include_indx_set.intersection(range(start_deletion, end_deletion)):
                deletion_positions.extend(range(start_deletion, end_deletion))
                deletion_coordinates.append((start_deletion, end_deletion))
                deletion_sizes.append(end_deletion - start_deletion)
            start_deletion = -1

    if start_deletion != -1:
        end_deletion = ref_positions[seq_len - 1]
        all_deletion_positions.extend(range(start_deletion, end_deletion))
        if include_indx_set.intersection(range(start_deletion, end_deletion)):
            deletion_positions.extend(range(start_deletion, end_deletion))
            deletion_coordinates.append((start_deletion, end_deletion))
            deletion_sizes.append(end_deletion - start_deletion)

    substitution_n = len(substitution_positions)
    deletion_n = sum(deletion_sizes)
    insertion_n = sum(insertion_sizes)

    return {
        'all_insertion_positions': all_insertion_positions,
        'all_insertion_left_positions': all_insertion_left_positions,
        'insertion_positions': insertion_positions,
        'insertion_coordinates': insertion_coordinates,
        'insertion_sizes': insertion_sizes,
        'insertion_n': insertion_n,

        'all_deletion_positions': all_deletion_positions,
        'deletion_positions': deletion_positions,
        'deletion_coordinates': deletion_coordinates,
        'deletion_sizes': deletion_sizes,
        'deletion_n': deletion_n,

        'all_substitution_positions': all_substitution_positions,
        'substitution_positions': substitution_positions,
        'all_substitution_values': np.array(all_substitution_values),
        'substitution_values': np.array(substitution_values),
        'substitution_n': substitution_n,

        'ref_positions': ref_positions,
    }

##def get_row_around_cut(row, cut_point, offset):##Changed by Xiaoli 20220809
def get_row_around_cut(row, cut_point, offset, ref_len):
    cut_idx=row['ref_positions'].index(cut_point)
    ##Changed by Xiaoli 20220809
    #return row['Aligned_Sequence'][cut_idx-offset+1:cut_idx+offset+1], row['Reference_Sequence'][cut_idx-offset+1:cut_idx+offset+1], row['Read_Status']=='UNMODIFIED', row['n_deleted'], row['n_inserted'], row['n_mutated'], row['#Reads'], row['%Reads']
    return row['Aligned_Sequence'][0:ref_len], row['Reference_Sequence'][0:ref_len], row['Read_Status']=='UNMODIFIED', row['n_deleted'], row['n_inserted'], row['n_mutated'], row['#Reads'], row['%Reads']



def get_dataframe_around_cut(df_alleles, cut_point,offset,collapse_by_sequence=True):
    if df_alleles.shape[0] == 0:
        return df_alleles
    ref1 = df_alleles['Reference_Sequence'].iloc[0]
    ref1 = ref1.replace('-', '')
    #if (cut_point + offset + 1 > len(ref1)):
    #    raise(BadParameterException('The plotting window cannot extend past the end of the amplicon. Amplicon length is ' + str(len(ref1)) + ' but plot extends to ' + str(cut_point+offset+1)))

    ##Changed by Xiaoli 20220809
    #df_alleles_around_cut=pd.DataFrame(list(df_alleles.apply(lambda row: get_row_around_cut(row, cut_point, offset), axis=1).values),columns=['Aligned_Sequence', 'Reference_Sequence', 'Unedited', 'n_deleted', 'n_inserted', 'n_mutated', '#Reads', '%Reads'])
    df_alleles_around_cut=pd.DataFrame(list(df_alleles.apply(lambda row: get_row_around_cut(row, cut_point, offset, len(ref1)), axis=1).values),columns=['Aligned_Sequence', 'Reference_Sequence', 'Unedited', 'n_deleted', 'n_inserted', 'n_mutated', '#Reads', '%Reads'])

    df_alleles_around_cut=df_alleles_around_cut.groupby(['Aligned_Sequence', 'Reference_Sequence', 'Unedited', 'n_deleted', 'n_inserted', 'n_mutated']).sum().reset_index().set_index('Aligned_Sequence')

    df_alleles_around_cut.sort_values(by='%Reads', inplace=True, ascending=False)
    df_alleles_around_cut['Unedited']=df_alleles_around_cut['Unedited']>0
    return df_alleles_around_cut
