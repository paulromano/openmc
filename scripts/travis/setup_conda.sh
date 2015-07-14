#!/bin/bash

if [[ "$TRAVIS_PYTHON_VERSION" == "2.7" ]]; then
  if [[ ! -d $HOME/miniconda ]]; then
    wget https://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh -O miniconda.sh
    bash miniconda.sh -b -p $HOME/miniconda
    export PATH="$HOME/miniconda/bin:$PATH"
  fi
else
  if [[ ! -d $HOME/miniconda3 ]]; then
    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
    bash miniconda.sh -b -p $HOME/miniconda3
    export PATH="$HOME/miniconda3/bin:$PATH"
  fi
fi

hash -r
conda config --set always_yes yes --set changeps1 no
conda update -q conda
conda info -a
conda create -q -n test-environment python=$TRAVIS_PYTHON_VERSION numpy scipy h5py
source activate test-environment
