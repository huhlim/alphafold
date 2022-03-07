#!/bin/bash

CUDA=11.5

ENVNAME=af2

conda create --name $ENVNAME -y -c conda-forge -c nvidia \
    python=3.7 \
    openmm=7.5.1 \
    cudatoolkit=$CUDA \
    pdbfixer \
    pip \
    absl-py=0.13.0  \
    cudnn \
    libcusolver 

source activate $ENVNAME

pip install --upgrade pip 
pip3 install -r requirements.txt 
# this need to be updated 
pip3 install --upgrade "jax[cuda11_cudnn82]" jaxlib -f \
      https://storage.googleapis.com/jax-releases/jax_releases.html

openmm_patch=$(readlink -f docker/openmm.patch)
cd $CONDA_PREFIX/lib/python3.7/site-packages
patch -p0 < $openmm_patch
cd -

if [[ ! -e alphafold/common/stereo_chemical_props.txt ]]; then
    wget -q -P alphafold/common/ \
        https://git.scicore.unibas.ch/schwede/openstructure/-/raw/7102c63615b64735c4941278d92b554ec94415f8/modules/mol/alg/src/stereo_chemical_props.txt
fi

