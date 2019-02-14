#!/bin/bash
set -x 

for idx in 0 1 2 3 4 5
do
    card=cards/mumu_fit_mbin${idx}.txt
    workspace=workspaces/mumu_fit${idx}.root
    cp mumu_fit_template.txt $card
    cat cards/mbin${idx}_bins.txt >> $card
    text2workspace.py $card --keyword-value M_BIN=${idx} -P Analysis.B2GTTrees.my_model:dy_AFB -o $workspace
    combine -M MultiDimFit $workspace --saveFitResult --saveWorkspace
    PostFitShapesFromWorkspace -w higgsCombineTest.MultiDimFit.mH120.root -f multidimfit.root:fit_mdf --postfit -o shapes/mumu_fit_shapes_mbin${idx}.root
    echo "fit_mdf->Print();" > cmd.txt
    echo ".q" >> cmd.txt
    root -l -b multidimfit.root < cmd.txt > fit_results/mumu_fit_results_mbin${idx}.txt
    rm combine_logger.out
    rm higgsCombineTest.MultiDimFit.mH120.root 
    rm multidimfit.root
done


