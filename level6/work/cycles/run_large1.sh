#!/bin/bash
# large count>1: cycles + prefix + bridge + tail (solve3f), seeded by current candidate programs
cd "/c/Users/domok/OneDrive/Desktop/ccc 2026/level6/work/cycles"
python publen.py 3-large
OMP_NUM_THREADS=4 ./solve3f.exe ../../in_level-6_3-large.txt pub_lengths_3-large.txt none sols/L_s3_$1.txt 15000 3 $2 $1.idx ../cand/cycles/out_level-6_3-large.txt > sols/L_s3_$1.log 2>&1
