#!/usr/bin/env bash
set -euo pipefail

echo "=== Configuring channels ==="
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --add channels defaults
echo "=== Configuring channels Done ==="

echo "=== Installing tools ==="

# PGGB 0.7.2
conda install bioconda::pggb=0.7.2

# fastix 0.1.0
cargo install fastix=0.1.0
cp .cargo/bin/fastix ~/miniconda3/bin/

# mash 2.3
conda install bioconda::mash=2.3

# samtools 1.21 + bgzip 1.21
conda install bioconda::htslib=1.21
conda install bioconda::samtools=1.21

# panacus 0.2.5
conda install bioconda::panacus=0.2.5

# odgi (attempt exact version)
conda install bioconda::odgi=0.9.0

# seqkit 2.6.0
conda install bioconda::seqkit=2.6.0

# Diamond 2.1.10
wget https://github.com/bbuchfink/diamond/releases/download/v2.1.10/diamond-linux64.tar.gz
tar -xzf diamond-linux64.tar.gz
mv diamond ~/miniconda3/bin/

conda install bioconda::repeatmasker=4.1.7

echo "===================================================="
echo " Environment setup completed!"
echo "===================================================="

