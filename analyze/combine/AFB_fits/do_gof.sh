#!/bin/bash

set -x 
idx=2
label=ee
label_num=1
gof_file=GoodnessOfFit/${label}_fit.txt
teststat=saturated
rm $gof_file
for idx in 0 1 2 3 4 5
do

    card=cards/${label}_fit_mbin${idx}.txt
    workspace=workspaces/${label}_fit${idx}.root
    text2workspace.py $card --keyword-value M_BIN=${idx} -P Analysis.B2GTTrees.my_model:dy_AFB -o $workspace --channel-masks
    combine -M GoodnessOfFit -d $workspace  --algo=${teststat} --setParametersForEval mask_emu=1,mask_ee_ss=1  #-S 0
    combine -M GoodnessOfFit -d $workspace  --algo=${teststat} -t 100 -s 123 --setParameters Afb=0.6 --setParametersForEval mask_emu=1,mask_ee_ss=1 #-S 0 --toysNoSystematics
    echo "mbin ${idx}" >> $gof_file
    echo ".x GoodnessOfFit/gof_helper.C($label_num, $idx)"  > cmds.txt
    echo ".q"  >> cmds.txt
    root -l -b < cmds.txt  >> $gof_file
    mv higgsCombineTest.GoodnessOfFit.mH120.123.root GoodnessOfFit/${label}_bin${idx}_toys.root
    rm higgsCombineTest.GoodnessOfFit.mH120.root

done
cat $gof_file



