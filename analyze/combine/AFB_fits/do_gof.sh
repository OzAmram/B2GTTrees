#!/bin/bash

set -x 
#idx=0
gof_file=GoodnessOfFit/combined_fit.txt
rm $gof_file
for idx in 0 1 2 3 4 5
do

    workspace=workspaces/combined_fit${idx}.root
    combine -M GoodnessOfFit -d $workspace  --algo=KS 
    combine -M GoodnessOfFit -d $workspace  --algo=KS -t 100 -s 123 
    echo "mbin ${idx}" >> $gof_file
    echo ".x GoodnessOfFit/gof_helper.C"  > cmds.txt
    echo ".q"  >> cmds.txt
    root -l -b < cmds.txt  >> $gof_file

done



