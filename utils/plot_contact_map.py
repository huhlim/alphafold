#!/usr/bin/env python

import os
import sys
import tqdm
import path
import numpy as np
import scipy.interpolate
import scipy.integrate
import pickle
import warnings
import matplotlib
import matplotlib.pyplot as plt
import argparse

warnings.filterwarnings("ignore")

font={"size": 14}
matplotlib.rc("font", **font)

D_MAX = 20.0
D_CNT =  8.0
d2c = lambda d: 1. - 1./(1.+np.exp(-2*(d-D_CNT-2)))

def read_pdb(pdb_fn):
    pdb = [] ; resNo_s = []
    chain_break_s = [0]
    with pdb_fn.open() as fp:
        resNo_prev_chain = 0
        for line in fp:
            if not line.startswith("ATOM"): 
                if line.startswith("TER"):
                    chain_break_s.append(resNo)
                continue
            #
            resName = line[17:20]
            atmName = line[12:16].strip()
            #
            if resName == 'GLY' and atmName != 'CA':
                continue
            elif resName != 'GLY' and atmName != 'CB':
                continue
            #
            resNo = int(line[22:26]) + chain_break_s[-1]
            r = np.array([line[30:38], line[38:46], line[46:54]], dtype=float)
            pdb.append(r)
            resNo_s.append(resNo)
    return np.array(pdb), np.array(resNo_s, dtype=int), chain_break_s

def get_contact_from_PDB(l_seq, pdb_fn):
    pdb, resNo_s, chain_break_s = read_pdb(pdb_fn)
    #
    status = np.ones(l_seq, dtype=bool)
    R = np.zeros((l_seq, 3), dtype=float)
    k = -1
    for resNo in range(l_seq):
        if np.isin(resNo+1, resNo_s):
            k += 1 ; R[resNo] = pdb[k]
        else:
            status[resNo] = False
    #
    dR = R[:,None] - R[None,:]
    dij = np.sqrt(np.sum(dR**2, axis=-1))
    return d2c(dij), chain_break_s

def run(pkl_fn, pdb_fn0, plot_mode):
    png_fn = '%s.png'%pkl_fn.name()
    npy_fn = '%s.contact_prob.npy'%pkl_fn.name()
    if os.path.exists(npy_fn):
        contact = np.load(npy_fn)
        l_seq = contact.shape[0]
    else:
        with pkl_fn.open("rb") as fp:
            X = pickle.load(fp)['distogram']
        dist_bin = X['bin_edges']
        dist = (dist_bin[1:] + dist_bin[:-1]) / 2
        logits_exp = np.exp(X['logits'])
        isinf = np.where(np.isinf(logits_exp))
        prob = logits_exp / logits_exp.sum(axis=-1)[...,None]
        prob[isinf[0],isinf[1]] = 0.
        prob[isinf] = 1.
        del logits_exp
        del X
        #
        l_seq = prob.shape[0]
        contact = np.zeros((l_seq, l_seq), dtype=np.float)
        n_pair = l_seq*(l_seq-1) // 2
        progress = tqdm.tqdm(total=n_pair)
        for i in range(l_seq):
            for j in range(i+1, l_seq):
                p0 = prob[i,j,0]
                p1 = prob[i,j,1:-1][dist < 8.0].sum()
                d0 = dist[dist<8.0][-1]
                func = scipy.interpolate.interp1d(dist, prob[i,j][1:-1], kind='cubic')
                p2 = scipy.integrate.quad(func, d0, 8.0, limit=10)[0]
                contact[i,j] = p0 + p1 + p2
                progress.update(1)
        contact += contact.T
        np.save(npy_fn, contact)
        progress.close()
    #
    if plot_mode > 0:
        model_number = pkl_fn.name().split("_")[-1]
        if pdb_fn0 is None:
            pdb_fn = path.Path('relaxed_model_%s.pdb'%model_number)
            pdb_fn_alt = path.Path('unrelaxed_model_%s.pdb'%model_number)
        else:
            pdb_fn = pdb_fn0
        if pdb_fn.status():
            pdb, chain_break_s = get_contact_from_PDB(l_seq, pdb_fn)
        elif pdb_fn_alt.status() and (pdb_fn0 is None):
            pdb, chain_break_s = get_contact_from_PDB(l_seq, pdb_fn_alt)
        else:
            sys.exit("cannot find PDB file, %s"%pdb_fn)
    #
    options = {"vmin": 0.2, "vmax": 0.5, "cmap": plt.get_cmap("hot_r"), "marker": 's'}
    #
    fig, ax = plt.subplots(figsize=(8.4, 7.2))
    if plot_mode == 0:
        X,Y = np.where(contact >= options['vmin'])
        ax.scatter(X+1, Y+1, s=1, alpha=1.0, c=contact[X,Y], **options)
    elif plot_mode == 1:
        X,Y = np.where(contact >= options['vmin'])
        lower_triangular = (X > Y)
        X = X[lower_triangular]
        Y = Y[lower_triangular]
        ax.scatter(X+1, Y+1, s=1, alpha=1.0, c=contact[X,Y], **options)
        #
        X,Y = np.where(pdb >= options['vmin'])
        upper_triangular = (X < Y)
        X = X[upper_triangular]
        Y = Y[upper_triangular]
        ax.scatter(X+1, Y+1, s=1, alpha=1.0, c=contact[X,Y], **options)
    elif plot_mode == 2:
        X,Y = np.where(contact >= options['vmin'])
        ax.scatter(X+1, Y+1, s=16, alpha=0.2, c='grey', marker='s')
        X,Y = np.where(pdb >= options['vmin'])
        ax.scatter(X+1, Y+1, s=1, alpha=1.0, c=pdb[X,Y], **options)
    #
    xylim = np.array([0.5, l_seq+0.5])
    ax.plot(xylim, xylim, 'k-')
    #
    ax.grid(True, linestyle='--', color='grey')
    ax.set_xlim(xylim)
    ax.set_ylim(xylim)
    ax.set_aspect("equal")
    #
    if plot_mode > 0:
        for chain_break in chain_break_s[1:-1]:
            ax.plot(xylim, np.full(2, chain_break), 'k-')
            ax.plot(np.full(2, chain_break), xylim, 'k-')
        if len(chain_break_s) > 2:
            ticks = ax.get_xticks()
            tick_space = ticks[1] - ticks[0]
            ticks = [] ; ticklabels = []
            for chain_start, chain_end in zip(chain_break_s[:-1], chain_break_s[1:]):
                tick = np.arange(chain_start, chain_end, tick_space)[1:]
                ticklabel = ['%d'%(resSeq-chain_start) for resSeq in tick]
                ticks.extend(tick.tolist())
                ticklabels.extend(ticklabel)
            #
            ax.set_xticks(ticks)
            ax.set_xticklabels(ticklabels)
            ax.set_yticks(ticks)
            ax.set_yticklabels(ticklabels)
    #
    fig.tight_layout()
    plt.savefig(png_fn)
    plt.close("all")

def main():
    arg = argparse.ArgumentParser(prog='plot_contact_map')
    arg.add_argument("--pkl", dest='pkl_fns', required=True, nargs='*')
    arg.add_argument("--pdb", dest='pdb_fns', nargs='*', default=[])
    arg.add_argument('--mode', dest='plot_mode', default=2, type=int, \
            help='plot_mode: 0 (predicted map only); 1 (triangular); 2 (overlap, default)')
    #
    if len(sys.argv) == 1:
        return arg.print_help()
    arg = arg.parse_args()
    #
    if len(arg.pdb_fns) == 1:
        arg.pdb_fns = [arg.pdb_fns[0] for _ in arg.pkl_fns]
    arg.pdb_fns = [path.Path(fn) for fn in arg.pdb_fns]
    #
    arg.pkl_fns = [path.Path(fn) for fn in arg.pkl_fns]
    #
    if len(arg.pdb_fns) == 0:
        for pkl_fn in arg.pkl_fns:
            run(pkl_fn, None, arg.plot_mode)
    else:
        for pkl_fn, pdb_fn in zip(arg.pkl_fns, arg.pdb_fns):
            run(pkl_fn, pdb_fn, arg.plot_mode)

if __name__ == '__main__':
    main()
