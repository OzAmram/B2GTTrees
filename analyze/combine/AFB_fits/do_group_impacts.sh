#!/bin/bash

#stolen from https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/tree/81x-root606/data/tutorials/groups

set -x 
for idx in 0 1 2 3 4 5
do
    workspace=workspaces/combined_fit${idx}.root
    combine -M MultiDimFit $workspace --robustFit 1 --saveWorkspace -n _nom  &> logs/nomfitinit_${idx}
    combine -M MaxLikelihoodFit  --robustFit 1 -d higgsCombine_nom.MultiDimFit.mH120.root -w w --snapshotName "MultiDimFit" -n _nom &> logs/nomfit_${idx}
    cp impacts/mbin${idx}.json impacts/group_mbin${idx}.json
    for par in emu_fake_shape ee_fake_shape mumu_fake_shape pdfs
        do
        combine -M MaxLikelihoodFit --freezeNuisanceGroups $par  --robustFit 1 -d higgsCombine_nom.MultiDimFit.mH120.root -w w --snapshotName "MultiDimFit" -n _${par} &> logs/${par}fit_${idx}
        ./add_group_impact.py higgsCombine_${par}.MaxLikelihoodFit.mH120.root higgsCombine_nom.MaxLikelihoodFit.mH120.root ${par} impacts/group_mbin${idx}.json
        done

    plotImpacts.py -i impacts/group_mbin${idx}.json -o impacts/impact_plot_group_mbin${idx} -t impacts/rename.json
    rm fitDiagnostics*
    rm higgsCombine*


done

