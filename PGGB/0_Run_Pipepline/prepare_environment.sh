#!/usr/bin/env bash

ENV_NAME="pggb"

echo "=== Creating conda environment: $ENV_NAME ==="
conda create -y -n $ENV_NAME python=3.10

echo "===================================================="
echo " Environment setup completed!"
echo " Activate it with:"
echo
echo "    conda activate $ENV_NAME"
echo
echo "===================================================="

