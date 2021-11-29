#!/usr/bin/env python

import os
import sys
import glob
import string
import argparse

def read_a3m(fn):
    a3m = []
    with open(fn) as fp:
        for line in fp:
            if line.startswith(">"):
                seq = [line.strip(), ""]
                a3m.append(seq)
            else:
                seq[1] += line.strip()
    return a3m

def read_sto(fn):
    name_s = []
    seq_s = []
    with open(fn) as fp:
        n_chunk = -1
        for line in fp:
            if line.startswith("#"):
                continue
            if line.startswith("//"):
                break
            if line.strip() == '':
                i = -1
                continue
            i += 1
            if i == 0: n_chunk += 1
            #
            x = line.strip().split()
            name = x[0] ; seq = x[1]
            if n_chunk == 0:
                name_s.append(name)
                seq_s.append("")
            seq_s[i] += seq
    #
    for i in range(len(seq_s)):
        seq_s[i] = seq_s[i].replace(".","")
    #
    a3m = []
    seqq = seq_s[0]
    for seqt in seq_s:
        a3m_seq = []
        for q,t in zip(seqq,seqt):
            if q != '-':
                a3m_seq.append(t)
            elif t != '-':
                a3m_seq.append(t.lower())
        a3m.append(''.join(a3m_seq))

    out = []
    for i,name in enumerate(name_s):
        out.append([name, a3m[i]])
    return out

def read_msa(msa_home):
    msa_s = []
    for fn in glob.glob("%s/*"%msa_home):
        if fn.endswith("a3m"):
            msa_s.extend(read_a3m(fn))
        elif fn.endswith("sto"):
            msa_s.extend(read_sto(fn))
    return msa_s

def main():
    arg = argparse.ArgumentParser(prog='build_oligomer_msa')
    arg.add_argument('-o', '--output', dest='output', default=None)
    arg.add_argument('-m', '--msa', dest='msa_home_s', required=True, nargs='*')
    arg.add_argument('-n', dest='n_mol', required=True, nargs='*', type=int)
    #
    if len(sys.argv) == 1:
        arg.print_help()
        return
    arg = arg.parse_args()
    if len(arg.msa_home_s) != len(arg.n_mol):
        sys.exit("ERROR: the numbers of msa_home and n_mol do not match!")
    #
    msa_s = []
    for n,msa_home in zip(arg.n_mol, arg.msa_home_s):
        msa_home = os.path.abspath(msa_home)
        if os.path.isdir(msa_home) and msa_home.split("/")[-1] != 'msas':
            msa_home += '/msas'
        sys.stderr.write("Reading ... %s\n"%msa_home)
        if os.path.isdir(msa_home):
            msa = read_msa(msa_home)
        else:
            msa = read_a3m(msa_home)
        if n > 1:
            msa = [[name,seq*n] for name,seq in msa]
        msa_s.append(msa)
    #
    out = []
    out.append(">Query\n")
    #
    seq = [] ; l_seq = []
    for msa in msa_s:
        seq.append(msa[0][1])
        l_seq.append(len(msa[0][1]))
    out.append("".join(seq) + '\n')
    #
    for i,msa in enumerate(msa_s):
        l_nt = sum(l_seq[:i])
        l_ct = sum(l_seq[i+1:])
        #
        for k,(name,seq) in enumerate(msa):
            if k == 0:
                continue
            out.append('>%s'%name + '\n')
            out.append('-'*l_nt + seq + '-'*l_ct + '\n')
    #
    if arg.output is not None:
        with open(arg.output, 'wt') as fout:
            fout.writelines(out)
    else:
        sys.stdout.writelines(out)

if __name__ == '__main__':
    main()

