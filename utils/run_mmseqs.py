#!/usr/bin/env python

import os
import sys
import time

import tempfile
import tarfile
import requests
from typing import Tuple, List
import argparse
import pathlib

import logging
logger = logging.getLogger(__name__)

HOST_URL = "https://a3m.mmseqs.com"

PROXY_SERVER = 'http://markov.bch.msu.edu:9999'
PROXY = {
        "http":  PROXY_SERVER,
        "https": PROXY_SERVER,
        "ftp":   PROXY_SERVER,
        }

def run_mmseqs2(x, path, use_env=True, use_filter=True,
                use_templates=False, use_pairing=False,
                ) -> Tuple[List[str], List[str]]:
  submission_endpoint = "ticket/pair" if use_pairing else "ticket/msa"

  def submit(seqs, mode, N=101):
    n, query = N, ""
    for seq in seqs:
      query += f">{n}\n{seq}\n"
      n += 1

    res = requests.post(f'{HOST_URL}/{submission_endpoint}', data={'q':query,'mode': mode}, proxies=PROXY)
    try:
      out = res.json()
    except ValueError:
      logger.error(f"Server didn't reply with json: {res.text}")
      out = {"status":"ERROR"}
    return out

  def status(ID):
    res = requests.get(f'{HOST_URL}/ticket/{ID}', proxies=PROXY)
    try:
      out = res.json()
    except ValueError:
      logger.error(f"Server didn't reply with json: {res.text}")
      out = {"status":"ERROR"}
    return out

  def download(ID, path):
    res = requests.get(f'{HOST_URL}/result/download/{ID}', proxies=PROXY)
    with open(path,"wb") as out: out.write(res.content)

  # process input x
  seqs = [x] if isinstance(x, str) else x

  # setup mode
  if use_filter:
    mode = "env" if use_env else "all"
  else:
    mode = "env-nofilter" if use_env else "nofilter"

  if use_pairing:
    mode = ""
    use_templates = False
    use_env = False

  # call mmseqs2 api
  tar_gz_file = f'{path}/out.tar.gz'
  N,REDO = 101,True

  # deduplicate and keep track of order
  seqs_unique = []
  #TODO this might be slow for large sets
  [seqs_unique.append(x) for x in seqs if x not in seqs_unique]
  Ms = [N + seqs_unique.index(seq) for seq in seqs]
  # lets do it!
  if not os.path.isfile(tar_gz_file):
    TIME_ESTIMATE = 150 * len(seqs_unique)
    while REDO:
      # Resubmit job until it goes through
      out = submit(seqs_unique, mode, N)
      while out["status"] in ["UNKNOWN", "RATELIMIT"]:
        sleep_time = 5 
        # resubmit
        time.sleep(sleep_time)
        out = submit(seqs_unique, mode, N)

      if out["status"] == "ERROR":
        raise Exception(f'MMseqs2 API is giving errors. Please confirm your input is a valid protein sequence. If error persists, please try again an hour later.')

      if out["status"] == "MAINTENANCE":
        raise Exception(f'MMseqs2 API is undergoing maintenance. Please try again in a few minutes.')

      # wait for job to finish
      ID,TIME = out["id"],0
      while out["status"] in ["UNKNOWN","RUNNING","PENDING"]:
        t = 5
        time.sleep(t)
        out = status(ID)
        if out["status"] == "RUNNING":
          TIME += t

      if out["status"] == "COMPLETE":
        REDO = False

      if out["status"] == "ERROR":
        REDO = False
        raise Exception(f'MMseqs2 API is giving errors. Please confirm your input is a valid protein sequence. If error persists, please try again an hour later.')

    # Download results
    download(ID, tar_gz_file)

  # prep list of a3m files
  if use_pairing:
    a3m_files = [f"{path}/pair.a3m"]
  else:
    a3m_files = [f"{path}/uniref.a3m"]
    if use_env: a3m_files.append(f"{path}/bfd.mgnify30.metaeuk30.smag30.a3m")

  # extract a3m files
  if any(not os.path.isfile(a3m_file) for a3m_file in a3m_files):
    with tarfile.open(tar_gz_file) as tar_gz:
      tar_gz.extractall(path)

  # templates
  if use_templates:
    templates = {}
    for line in open(f"{path}/pdb70.m8","r"):
      p = line.rstrip().split()
      M,pdb,qid,e_value = p[0],p[1],p[2],p[10]
      M = int(M)
      if M not in templates: templates[M] = []
      templates[M].append(pdb)

    template_paths = {}
    for k,TMPL in templates.items():
      TMPL_PATH = f"{path}/templates_{k}"
      if not os.path.isdir(TMPL_PATH):
        os.mkdir(TMPL_PATH)
        TMPL_LINE = ",".join(TMPL[:20])
        os.system(f"curl -s --proxy {PROXY_SERVER} https://a3m-templates.mmseqs.com/template/{TMPL_LINE} | tar xzf - -C {TMPL_PATH}/")
        os.system(f"cp {TMPL_PATH}/pdb70_a3m.ffindex {TMPL_PATH}/pdb70_cs219.ffindex")
        os.system(f"touch {TMPL_PATH}/pdb70_cs219.ffdata")
      template_paths[k] = TMPL_PATH

  # gather a3m lines
  a3m_lines = {}
  for a3m_file in a3m_files:
    update_M,M = True,None
    for line in open(a3m_file,"r"):
      if len(line) > 0:
        if "\x00" in line:
          line = line.replace("\x00","")
          update_M = True
        if line.startswith(">") and update_M:
          M = int(line[1:].rstrip())
          update_M = False
          if M not in a3m_lines: a3m_lines[M] = []
        a3m_lines[M].append(line)

  # return results

  a3m_lines = ["".join(a3m_lines[n]) for n in Ms]

  if use_templates:
    template_paths_ = []
    for n in Ms:
      if n not in template_paths:
        template_paths_.append(None)
        print(f"{n-N}\tno_templates_found")
      else:
        template_paths_.append(template_paths[n])
    template_paths = template_paths_
  else:
    template_paths = [None for _ in a3m_lines]
  return path, a3m_lines, template_paths

def main():
    arg = argparse.ArgumentParser(prog='run_mmseqs')
    arg.add_argument(dest='fa_fn_s', nargs='+')
    arg.add_argument('--env', dest='use_env', default=True, action='store_true')
    arg.add_argument('--noenv', dest='use_env', default=True, action='store_false')
    arg.add_argument('--templates', dest='use_templates', default=False, action='store_true')
    arg.add_argument('--notemplates', dest='use_templates', default=False, action='store_false')
    arg.add_argument('--pair', dest='use_pairing', default=False, action='store_true')
    arg.add_argument('--nopair', dest='use_pairing', default=False, action='store_false')
    if len(sys.argv) == 1:
        return arg.print_help()
    arg = arg.parse_args()

    def read_fasta(fn):
        seq_s = []
        with open(fn) as fp:
            for line in fp:
                if line.startswith(">"):
                    seq = [] ; seq_s.append(seq)
                else:
                    seq.append(line.strip())
        return [''.join(seq) for seq in seq_s if len(seq) > 0]
    #
    name_s = []
    seq_s = []
    for fa_fn in arg.fa_fn_s:
        seq = read_fasta(fa_fn)
        n_seq = len(seq)
        name = pathlib.Path(fa_fn).stem
        name_s.extend([name for _ in range(n_seq)])
        seq_s.extend(seq)
    #
    tmpdir = tempfile.TemporaryDirectory(prefix='mmseqs.')
    path, a3m_s, template_paths = run_mmseqs2(seq_s, tmpdir.name, \
            use_env=arg.use_env, \
            use_templates=arg.use_templates, \
            use_pairing=arg.use_pairing)
    #
    for name, a3m, template_path in zip(name_s, a3m_s, template_paths):
        with open(f"{name}.mmseqs.a3m", 'wt') as fout:
            fout.write(a3m)
        if template_path is not None:
            template_dir = f"{name}.mmseqs.pdb70"
            if os.path.exists(template_dir):
                os.system(f"rm -rf {template_dir}")
            os.system(f"cp -rf {template_path} {template_dir}")

if __name__ == '__main__':
    main()
