#!/usr/bin/env python

import os
import sys
import json
import pickle
import pathlib

def read_rank(fn):
    with open(fn) as fp:
        rank_s = json.load(fp)
    return rank_s

def read_residue_error(pkl_fn):
    with open(pkl_fn, 'rb') as fp:
        pkl = pickle.load(fp)
    return pkl['plddt']

def apply_error(out_fn, pdb_fn, plddt_global, residue_error):
    wrt = []
    wrt.append("REMARK  plddt_global    %6.2f\n"%plddt_global)
    with open(pdb_fn, 'r') as fp:
        for line in fp:
            if not line.startswith("ATOM"):
                wrt.append(line)
                continue
            resNo = int(line[22:26])
            i_res = resNo-1
            wrt.append("%s%6.2f%s"%(line[:60], residue_error[i_res], line[66:]))
    with open(out_fn, 'wt') as fout:
        fout.writelines(wrt)

def main():
    if len(sys.argv) == 1:
        run_home = pathlib.Path('.').resolve()
    else:
        run_home = pathlib.Path(sys.argv[1]).resolve()
    #
    ranking_json = run_home / "ranking_debug.json"
    if not ranking_json.exists():
        sys.exit("file does not exist: %s"%ranking_json)
    rank_s = read_rank(ranking_json)
    #
    for i,rank in enumerate(rank_s['order']):
        plddt_global = rank_s['plddts'][rank]
        #
        pkl_fn = run_home / ("result_%s.pkl"%rank)
        if not pkl_fn.exists():
            sys.exit("file does not exist: %s"%pkl_fn)
        residue_error = read_residue_error(pkl_fn)
        #
        pdb_fn = run_home / ("ranked_%d.pdb"%i)
        out_fn = run_home / ("model_%d.pdb"%(i+1))
        apply_error(out_fn, pdb_fn, plddt_global, residue_error)

if __name__ == '__main__':
    main()

