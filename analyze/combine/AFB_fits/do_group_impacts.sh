#!/bin/bash

#stolen from https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/tree/81x-root606/data/tutorials/groups

set -x 
#idx=0
#for idx in 0 1 2 3 4 5
idx=0
#do
    workspace=workspaces/combined_fit${idx}.root
    combine -M MultiDimFit $workspace --robustFit 1 --saveWorkspace -n _nom  &> logs/nomfit_${idx}
    combine -M MaxLikelihoodFit --robustFit 1 -d higgsCombine_nom.MultiDimFit.mH120.root -w w --snapshotName "MultiDimFit" -v9 -n _${par} &> logs/${par}fit_${idx}
    for par in emu_fake_shape  ee_fake_shape mumu_fake_shape pdfs
        do
        combine -M MaxLikelihoodFit --freezeNuisanceGroups $par  --robustFit 1 -d higgsCombine_nom.MultiDimFit.mH120.root -w w --snapshotName "MultiDimFit" -n _${par} &> logs/${par}fit_${idx}

        done
#done

