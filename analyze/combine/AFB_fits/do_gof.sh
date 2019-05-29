#!/bin/bash

set -x 
idx=0
label=combined
gof_file=GoodnessOfFit/${label}_fit.txt
teststat=AD
rm $gof_file
#for idx in 0 1 2 3 4 5
#do

    workspace=workspaces/${label}_fit${idx}.root
    combine -M GoodnessOfFit -d $workspace  --algo=${teststat} #-S 0
    combine -M GoodnessOfFit -d $workspace  --algo=${teststat} -t 100 -s 123 --expectSignal 0.6 #-S 0 --toysNoSystematics
    echo "mbin ${idx}" >> $gof_file
    echo ".x GoodnessOfFit/gof_helper.C"  > cmds.txt
    echo ".q"  >> cmds.txt
    root -l -b < cmds.txt  >> $gof_file

#done
cat $gof_file



