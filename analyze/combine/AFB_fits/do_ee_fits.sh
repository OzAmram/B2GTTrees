#!/bin/bash
set -x 

for idx in 0 1 2 3 4 5
do
    card=cards/ee_fit_mbin${idx}.txt
    workspace=workspaces/ee_fit${idx}.root
    plotdir=postfit_plots/ee_mbin${idx}
    cp ee_fit_template.txt $card
    cat cards/mbin${idx}_bins.txt >> $card
    text2workspace.py $card --keyword-value M_BIN=${idx} -P Analysis.B2GTTrees.my_model:dy_AFB -o $workspace --channel-masks
    combine -M FitDiagnostics $workspace --saveShapes --saveWorkspace --skipBOnlyFit --plot --robustFit 1
    PostFitShapesFromWorkspace -w higgsCombineTest.FitDiagnostics.mH120.root -f fitDiagnostics.root:fit_s --postfit -o shapes/ee_fit_shapes_mbin${idx}.root
    rm -r $plotdir
    mkdir $plotdir
    mv *.png $plotdir
    echo "fit_s->Print();" > cmd.txt
    echo ".q" >> cmd.txt
    root -l -b fitDiagnostics.root < cmd.txt > fit_results/ee_fit_results_mbin${idx}.txt
    rm combine_logger.out
    rm higgsCombineTest.FitDiagnostics.mH120.root 
    rm fitDiagnostics.root
done


