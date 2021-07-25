#!/usr/bin/env python

import sys
import string

REMOVE_INSERTION = str.maketrans("", "", string.ascii_lowercase)

def read_a3m(fn):
    a3m = []
    with open(fn) as fp:
        for line in fp:
            if line.startswith(">"):
                seq = [line.strip(), ""]
                a3m.append(seq)
            else:
                seq[1] += line.strip().translate(REMOVE_INSERTION)
    return a3m

def main():
    a3m_fn = sys.argv[1]
    n_mer = int(sys.argv[2])
    #
    a3m_in = read_a3m(a3m_fn)
    l_seq_per_chain = len(a3m_in[0][1])
    #
    a3m_out = []
    for k in range(n_mer):
        l_nt = l_seq_per_chain * k
        l_ct = l_seq_per_chain * (n_mer - k -1)
        for name,seq in a3m_in:
            a3m_out.append(name)
            a3m_out.append('-'*l_nt + seq + '-'*l_ct)
    #
    for line in a3m_out:
        sys.stdout.write("%s\n"%line)

if __name__ == '__main__':
    main()

