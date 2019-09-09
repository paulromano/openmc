#!/bin/bash
set -ex

# Download HDF5 data
if [[ ! -e $HOME/endfb71_hdf5/cross_sections.xml ]]; then
    wget -q -O - https://anl.box.com/shared/static/9igk353zpy8fn9ttvtrqgzvw1vtejoz6.xz | tar -C $HOME -xJ
fi

# Download ENDF/B-VII.1 distribution
ENDF=$HOME/endf-b-vii.1/
if [[ ! -d $ENDF/neutrons || ! -d $ENDF/photoat || ! -d $ENDF/atomic_relax ]]; then
    wget -q -O - https://anl.box.com/shared/static/4kd2gxnf4gtk4w1c8eua5fsua22kvgjb.xz | tar -C $HOME -xJ
fi
