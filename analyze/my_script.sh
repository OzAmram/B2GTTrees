#!/bin/bash
set -x

cd Analysis/B2GTTrees/analyze
mkdir output_files
echo ".x MuMu/MuMu_reco_mc_batch.C++($2,$3)" > cmd.txt
echo ".q" >> cmd.txt
root -l -b < cmd.txt
xrdcp output_files/*.root $1
