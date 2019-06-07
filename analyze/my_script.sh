#!/bin/bash


set -x

ls -la
pwd
cd Analysis/B2GTTrees/analyze
mkdir output_files
#++ always goes after the .C, before args
echo '.x MuMu/MuMu_reco_background_batch.C++('$2','$3', "EOS_files/2016/combined_back_files_may29.txt");' > cmd.txt
echo ".q" >> cmd.txt
root -l -b < cmd.txt
xrdcp -f output_files/*.root $1
