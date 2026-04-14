#!/usr/bin/env bash

echo "=== Configuring channels ==="
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --add channels defaults
echo "=== Configuring channels Done ==="

echo "=== Installing tools ==="

# PGGB 0.7.2
conda install bioconda::pggb

# fastix 0.1.0
cargo install fastix
cp .cargo/bin/fastix ~/miniconda3/bin/

# mash 2.3
conda install bioconda::mash

# samtools 1.21 + bgzip 1.21
conda install bioconda::htslib
conda install bioconda::samtools

# panacus 0.2.5
conda install bioconda::panacus

# odgi (attempt exact version)
conda install bioconda::odgi

# seqkit 2.6.0
conda install bioconda::seqkit

# Diamond 2.1.10
wget https://github.com/bbuchfink/diamond/releases/download/v2.1.10/diamond-linux64.tar.gz
tar -xzf diamond-linux64.tar.gz
mv diamond ~/miniconda3/bin/

conda install bioconda::repeatmasker

echo "===================================================="
echo " Environment setup completed!"
echo "===================================================="

