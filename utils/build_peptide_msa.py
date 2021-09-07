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
    pept_fn = sys.argv[2]
    #
    a3m_in = read_a3m(a3m_fn)
    pept_seq = read_a3m(pept_fn)
    l_pept_seq = len(pept_seq[0][1])
    #
    a3m_out = []
    for i,(name,seq) in enumerate(a3m_in):
        a3m_out.append(name)
        if i == 0:
            a3m_out.append(seq + pept_seq[0][1])
        else:
            a3m_out.append(seq + '-'*l_pept_seq)
    #
    for line in a3m_out:
        sys.stdout.write("%s\n"%line)



if __name__ == '__main__':
    main()

