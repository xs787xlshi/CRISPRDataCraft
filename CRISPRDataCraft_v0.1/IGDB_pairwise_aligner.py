#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
IGDB_pairwise_aligner.py
Author: Xiaoli Shi
Date: 2024.10.14
Short reads alignment
(c) 2024 IGDB. All Rights Reserved.
'''

import random
import math

def repeat(v, n):
    for i in range(n):
        yield v


class CmpSeq:
    # comp_table = {'a': 'c', 'c': 'a', 'g': 't', 't': 'g'}
    comp_table = {'A': 'C', 'C': 'A', 'G': 'T', 'T': 'G', 'N': 'N'}
    n_changes = 5

    def __init__(self, ref_seq, m_seq, n_retry):
        self.ref_seq = ref_seq
        self.ref_seq_comp = self.comp_seq(ref_seq)
        self.m_seq = m_seq
        self.n_retry = n_retry
        self.cmp_result = None
        self.compl = "+"

    def comp_seq(self, seq):
        return ''.join([self.comp_table[x] for x in reversed(seq)])

    def cmp_aux(self, seq1, seq2):
        n1 = len(seq1)
        n2 = len(seq2)
        n = min(n1, n2)
        m = int(math.log2(n))
        # m = max(m, n // self.n_changes)

        idx = -1
        for i in range(max(m, self.n_retry)):
            # o = random.randint(0, n - m)
            o = random.randint(0, n2 - m)
            t = seq2[o:o + m]
            # print(t)
            idx = seq1.find(t)
            # print(idx)
            if idx >= 0:
                # print(i)
                break
        if idx < 0:
            return None
        s1 = idx
        s2 = o
        while s1 > 0 and s2 > 0: # start expanding to left
            if seq1[s1 - 1] != seq2[s2 - 1]:
                break
            s1 -= 1
            s2 -= 1
        e1 = idx + m - 1
        e2 = o + m - 1
        while e1 < (n1 - 1) and e2 < (n2 - 1): # start expanding to right
            if seq1[e1 + 1] != seq2[e2 + 1]:
                break
            e1 += 1
            e2 += 1
        # print(s1, e1, s2, e2)
        if s1 >= 2 and s2 >= 2: # there are at least 3 unaligned nt on the left side of seq1 and seq2
            rl = self.cmp_aux(seq1[0:s1], seq2[0:s2])
        else: # end of left
            if s1 == 2 and s2 == 2:
                if seq1[1] == seq2[0]:
                    rl = [((1, 1), (0, 0))]
                elif seq1[0] == seq2[1]:
                    rl = [((0, 0), (1, 1))]
                elif seq1[0] == seq2[0]:
                    rl = [((0, 0), (0, 0))]
                else:
                    rl = None
            elif s1 == 2 and s2 == 1 or s1 == 1 and s2 == 2:
                if seq1[0] == seq2[0]:
                    rl = [((0, 0), (0, 0))]
                else:
                    rl = None
            else:
                rl = None

        # start to deal with the remaining sequence at right
        nrem1 = n1 - e1 - 1
        nrem2 = n2 - e2 - 1
        if nrem1 >= 2 and nrem2 >= 2: # if the length of remaining seq > 3, call cmp_aux
            rr = self.cmp_aux(seq1[e1 + 1:], seq2[e2 + 1:])
        else: # end of right
            if nrem1 == 2 and nrem2 == 2:
                if seq1[n1 - 1] == seq2[n2 - 2]:
                    rr = [((nrem1 - 1, nrem1 - 1), (nrem2 - 2, nrem2 - 2))]
                elif seq1[n1 - 2] == seq2[n2 - 1]:
                    rr = [((nrem1 - 2, nrem1 - 2), (nrem2 - 1, nrem2 - 1))]
                elif seq1[n1 - 1] == seq2[n2 - 1]:
                    rr = [((nrem1 - 1, nrem1 - 1), (nrem2 - 1, nrem2 - 1))]
                else:
                    rr = None
            # right end ---- this bug has been corrected in version 0.5
            elif nrem1 == 2 and nrem2 == 1 or nrem1 == 1 and nrem2 == 2:
                if seq1[n1 - 1] == seq2[n2 - 1]:
                    rr = [((nrem1 - 1, nrem1 - 1), (nrem2 - 1, nrem2 - 1))]
                else:
                    rr = None
            else:
                rr = None
            pass
        ret = []
        if rl:
            ret.extend(rl)
        ret.append(((s1, e1), (s2, e2)))
        if rr:
            ret.extend([((x[0][0] + e1 + 1, x[0][1] + e1 + 1),
                         (x[1][0] + e2 + 1, x[1][1] + e2 + 1)) for x in rr])
        return ret

    def cmp_v1(self):
        ret = self.cmp_aux(self.ref_seq, self.m_seq)
        if not ret:
            ret = self.cmp_aux(self.ref_seq_comp, self.m_seq)
        self.cmp_result = ret
        return ret

    def cmp(self):
        ret = None
        for i in range(self.n_retry):
            ret = self.cmp_aux(self.ref_seq, self.m_seq)
            if ret:
                n = sum(e[0][1] - e[0][0] + 1 for e in ret)
                if n >= min(len(self.ref_seq), len(self.m_seq)) - self.n_changes:
                    break
        if not ret:
            for i in range(self.n_retry):
                ret = self.cmp_aux(self.ref_seq_comp, self.m_seq)
                if ret:
                    n = sum(e[0][1] - e[0][0] + 1 for e in ret)
                    if n >= min(len(self.ref_seq), len(self.m_seq)) - self.n_changes:
                        self.compl = "-"
                        break
        self.cmp_result = ret
        return ret, self.compl

    def __str__(self):
        if self.cmp_result is None:
            self.cmp()
        if not self.cmp_result:
            return ''
        result = self.cmp_result

        p1 = 0
        p2 = 0
        for i in range(len(self.cmp_result)):
            pass

    @classmethod
    def gen_seq(cls, length, n_mut, n_inser, n_del):
        keys = list(cls.comp_table.keys())
        seq = []
        for i in range(length):
            i = random.randint(0, len(keys) - 1)
            seq.append(keys[i])
        ref_seq = ''.join(seq)

        for i in range(n_mut):
            o = random.randint(0, len(seq) - 1)
            seq[o] = cls.comp_table[seq[o]]
        for i in range(n_inser):
            o = random.randint(0, len(seq) - 1)
            seq.insert(o, keys[o % len(keys)])
        for i in range(n_del):
            o = random.randint(0, len(seq) - 1)
            seq[o:o + 1] = []
        seq = ''.join(seq)

        return (ref_seq, seq)

    @classmethod
    def test_run(cls, n_mut = 4, n_inser = 4, n_del = 4, verbose=True, idx=0):
        cls.n_changes = n_mut + n_inser + n_del
        s1, s2 = CmpSeq.gen_seq(300, n_mut, n_inser, n_del)
        if verbose:
            print(s1)
            print(s2)

        cmp = CmpSeq(s1, s2, 5)
        r = cmp.cmp()
        if verbose:
            print(r)
            e0 = 0
            e1 = 0
            for e in r:
                print('<' + s1[e0:e[0][0]] + '>')
                print('<' + s2[e1:e[1][0]] + '>')
                print(s1[e[0][0]:e[0][1] + 1])
                print(s2[e[1][0]:e[1][1] + 1])
                e0 = e[0][1] + 1
                e1 = e[1][1] + 1

            nr = 0
            for e in r:
                nr += e[0][1] - e[0][0] + 1
            print(nr)
        return idx

    @classmethod
    def test_run_repeat(cls, n_mut, n_inser, n_del, n, verbose=True):
        for i in range(n):
            cls.test_run(n_mut, n_inser, n_del, verbose)

    @classmethod
    def test_run_thread_pool(cls, n_mut, n_inser, n_del, n, verbose=True):
        import concurrent.futures
        with concurrent.futures.ThreadPoolExecutor() as executor:
            futures = {
                executor.submit(cls.test_run, n_mut, n_inser, n_del, verbose): i
                for i in range(n)
            }
            for f in concurrent.futures.as_completed(futures):
                # print(f'Task {futures[f]} finished')
                pass

    @classmethod
    def test_run_thread_pool_2(cls, n_mut, n_inser, n_del, n, verbose=True):
        import concurrent.futures
        import multiprocessing
        with concurrent.futures.ThreadPoolExecutor() as executor:
            args = [n_mut, n_inser, n_del, verbose] * n
            cs = max(10, n // multiprocessing.cpu_count())
            # print(f'{cs=}')
            for i, f in zip(
                    range(n),
                    # executor.map(cls.test_run, args, range(n), chunksize=cs)):
                    executor.map(cls.test_run, repeat(n_mut, n), repeat(n_inser, n), repeat(n_del, n), repeat(verbose, n), range(n), chunksize=cs)):
                # print(f'Task {i} {f} finished')
                pass

    @classmethod
    def test_run_process_pool(cls, n_mut, n_inser, n_del, n, verbose=True):
        import concurrent.futures
        with concurrent.futures.ProcessPoolExecutor() as executor:
            futures = {
                executor.submit(cls.test_run, n_mut, n_inser, n_del, verbose): i
                for i in range(n)
            }
            for f in concurrent.futures.as_completed(futures):
                # print(f'Task {futures[f]} finished')
                pass

    @classmethod
    def test_run_process_pool_2(cls, n_mut, n_inser, n_del, n, verbose=True):
        import concurrent.futures
        import multiprocessing
        with concurrent.futures.ProcessPoolExecutor() as executor:
            args = [n_mut, n_inser, n_del, verbose] * n
            cs = max(10, n // multiprocessing.cpu_count())
            # print(f'{cs=}')
            for i, f in zip(
                    range(n),
                    # executor.map(cls.test_run, args, range(n), chunksize=cs)):
                    executor.map(cls.test_run, repeat(n_mut, n), repeat(n_inser, n), repeat(n_del, n), repeat(verbose, n), range(n), chunksize=cs)):
                # print(f'Task {i} {f} finished')
                pass


if __name__ == '__main__':
    CmpSeq.test_run(5, 5, 5)
    '''
    s1, s2 = CmpSeq.gen_seq(300, 4, 4, 4)
    print(s1)
    print(s2)

    cmp = CmpSeq(s1, s2, 5)
    r = cmp.cmp()
    print(r)
    for e in r:
        print(s1[e[0][0]:e[0][1] + 1])
        print(s2[e[1][0]:e[1][1] + 1])
    '''
