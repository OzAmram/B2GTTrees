#!/bin/bash
set -x

ls -la
pwd
cd Analysis/B2GTTrees/analyze
mkdir output_files
echo ".x FakeRate/ElEl_QCD_MC.C ($2,$3)" > cmd.txt
echo ".q" >> cmd.txt
root -l -b < cmd.txt
xrdcp -f output_files/*.root $1
