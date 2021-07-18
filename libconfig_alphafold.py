#!/usr/bin/env python

import os

jackhmmer_binary_path = "PATH_TO_JACKHMMER"
hhblits_binary_path = "PATH_TO_HHBLITS"
hhsearch_binary_path = "PATH_TO_HHSEARCH"
kalign_binary_path = "PATH_TO_KALIGN"

DOWNLOAD_DIR = 'PATH_TO_DATABASE_BASE'

data_dir = DOWNLOAD_DIR

# Path to the Uniref90 database for use by JackHMMER.
uniref90_database_path = os.path.join(
    DOWNLOAD_DIR, 'uniref90', 'uniref90.fasta')

# Path to the MGnify database for use by JackHMMER.
mgnify_database_path = os.path.join(
    DOWNLOAD_DIR, 'mgnify', 'mgy_clusters.fa')

# Path to the BFD database for use by HHblits.
bfd_database_path = os.path.join(
    DOWNLOAD_DIR, 'bfd',
    'bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt')

# Path to the Uniclust30 database for use by HHblits.
uniclust30_database_path = os.path.join(
    DOWNLOAD_DIR, 'uniclust30', 'uc30')
    #DOWNLOAD_DIR, 'uniclust30', 'uniclust30_2018_08', 'uniclust30_2018_08')

# Path to the PDB70 database for use by HHsearch.
pdb70_database_path = os.path.join(DOWNLOAD_DIR, 'pdb70', 'pdb70')

# Path to a directory with template mmCIF structures, each named <pdb_id>.cif')
template_mmcif_dir = os.path.join(DOWNLOAD_DIR, 'pdb_mmcif', 'mmcif_files')

# Path to a file mapping obsolete PDB IDs to their replacements.
obsolete_pdbs_path = os.path.join(DOWNLOAD_DIR, 'pdb_mmcif', 'obsolete.dat')

max_template_date = '2099-12-31'

os.environ['NVIDIA_VISIBLE_DEVICES'] = os.getenv("CUDA_VISIBLE_DEVICES", "")
os.environ['TF_FORCE_UNIFIED_MEMORY'] = '1'
os.environ['XLA_PYTHON_CLIENT_MEM_FRACTION'] = '4.0'
if os.getenv("CUDA_VISIBLE_DEVICES", "") == "":
    os.environ['JAX_PLATFORM_NAME'] = 'cpu'

model_names = [
    'model_1',
    'model_2',
    'model_3',
    'model_4',
    'model_5',
]
