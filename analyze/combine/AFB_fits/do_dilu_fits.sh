#!/bin/bash
set -x 

for idx in 0 1 2 3 4 5
do
    for dilu_bin in 0 1 2 3 4 5 6 7 8 9 10
    do
    #idx=1
    #dilu_bin=10
        card=cards/dilu_fit_mbin${idx}.txt
        workspace=workspaces/dilu_fit${idx}.root
        cp dilus_fit_template.txt $card
        cat cards/mbin${idx}_bins.txt >> $card
        sed -i "s/DLP/${dilu_bin}/" $card
        text2workspace.py $card --keyword-value M_BIN=${idx}  -P Analysis.B2GTTrees.my_model:dy_AFB -o $workspace
        combine -M FitDiagnostics $workspace --setParameters Dilu${dilu_bin}=1.0 --freezeParameters Dilu${dilu_bin}
        echo "fit_s->Print();" > cmd.txt
        echo ".q" >> cmd.txt
        mkdir fit_results/dilu_bin${dilu_bin}
        root -l -b fitDiagnostics.root < cmd.txt > fit_results/dilu_bin${dilu_bin}/mbin${idx}.txt
        rm combine_logger.out
        rm higgsCombineTest.FitDiagnostics.mH120.root 
        rm fitDiagnostics.root
    done
done


