#!/bin/bash

# MAKE SURE YOU HAVE bioawk IN YOUR ENVIRONMENT BEFORE RUNNING THIS CSCRIPT
# preface script with `source` to run the script with conda_env


conda activate bioawk

# Get input file
inputFastaFile=$1


# Filtering for lengths greater than 10kb and cleaning up file
bioawk -c fastx '{print $name, length($seq)}' $inputFastaFile | sed 's|chromosome_||g'

# Save output to a text file
