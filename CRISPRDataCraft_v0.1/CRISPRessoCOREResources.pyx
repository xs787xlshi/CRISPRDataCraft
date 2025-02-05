import cython
import numpy as np
cimport numpy as np
import re

cdef extern from "stdlib.h":
    ctypedef unsigned int size_t
    size_t strlen(char* s)


cdef extern from "Python.h":
    ctypedef void PyObject
    int _PyBytes_Resize(PyObject **, size_t)
    char * PyBytes_AS_STRING(PyObject *)


re_find_indels = re.compile("(-*-)")

@cython.boundscheck(False)
@cython.nonecheck(False)
@cython.wraparound(False)
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
    all_substitution_positions=[]
    substitution_positions=[]
    all_substitution_values=[]
    substitution_values=[]

    all_deletion_positions = []
    deletion_positions = []
    deletion_coordinates = []
    deletion_sizes = []
    cdef int start_deletion = -1  # the -1 value indicates that there currently isn't a deletion

    all_insertion_positions = []
    all_insertion_left_positions = []
    insertion_positions = []
    insertion_coordinates = []
    insertion_sizes = []
    cdef int start_insertion = -1  # the -1 value indicates that there currently isn't an insertion

    cdef size_t seq_len = len(ref_seq_al)
    include_indx_set = set(_include_indx)
    nucSet = set(['A', 'T', 'C', 'G', 'N'])
    cdef int idx = 0
    cdef int idx_c
    cdef int current_insertion_size = 0
    for idx_c, c in enumerate(ref_seq_al):
        #print("ref_seq_al ", str(idx_c), str(c))
        if c != '-':
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
            #print("start a deletion ", str(start_deletion))
        elif read_seq_al[idx_c] != '-' and start_deletion != -1:  # this is the end of a deletion
            end_deletion = abs(ref_positions[idx_c])
            all_deletion_positions.extend(range(start_deletion, end_deletion))
            '''
            print("end_deletion ", end_deletion)
            print("all_deletion_positions ", all_deletion_positions)
            print("out of include_indx_set ", str(start_deletion), str(end_deletion))
            '''
            if include_indx_set.intersection(range(start_deletion, end_deletion)):
                #print("include_indx_set ", str(start_deletion), str(end_deletion))
                deletion_positions.extend(range(start_deletion, end_deletion))
                deletion_coordinates.append((start_deletion, end_deletion))
                deletion_sizes.append(end_deletion - start_deletion)
            start_deletion = -1
            #print("end a deletion")
            #print(deletion_positions)

    if start_deletion != -1:
        end_deletion = ref_positions[seq_len - 1]
        all_deletion_positions.extend(range(start_deletion, end_deletion))
        if include_indx_set.intersection(range(start_deletion, end_deletion)):
            deletion_positions.extend(range(start_deletion, end_deletion))
            deletion_coordinates.append((start_deletion, end_deletion))
            deletion_sizes.append(end_deletion - start_deletion)

    cdef size_t substitution_n = len(substitution_positions)
    cdef size_t deletion_n = sum(deletion_sizes)
    cdef size_t insertion_n = sum(insertion_sizes)
    #print('all_deletion_positions: ', str(all_deletion_positions),'deletion_positions: ', str(deletion_positions),'deletion_coordinates: ', str(deletion_coordinates), 'deletion_sizes: ', str(deletion_sizes),'deletion_n:', str(deletion_n))

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


def calculate_homology(a, b):
    cdef char *al = a
    cdef char *bl = b
    cdef size_t l = strlen(al)
    cdef float score = 0.0

    for i in range(l):
        if al[i] == bl[i]:
            score+=1
    return score/l
