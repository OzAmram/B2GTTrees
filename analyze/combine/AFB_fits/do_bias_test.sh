#!/bin/bash

set -x 
#idx=1
nToys=100
chan=ee
label=100_newpdf
freezeParam=pdf,bk_xsec,lumi
for idx in 0 1 2 3 4 5
do

    workspace=workspaces/${chan}_fit${idx}.root

    #combine -M GenerateOnly -d $workspace --toysFrequentist -t $nToys --expectSignal 0.6 --saveToys 
    combine -M GenerateOnly -d $workspace --toysNoSystematics -t $nToys --expectSignal 0.6 --saveToys
    #combine -M FitDiagnostics -d $workspace --toysFile higgsCombineTest.GenerateOnly.mH120.123456.root  -t $nToys --freezeParameters $freezeParam
    combine -M FitDiagnostics -d $workspace --toysFile higgsCombineTest.GenerateOnly.mH120.123456.root  -t $nToys 
    echo 'tree_fit_sb->Draw("(Afb-0.6)/AfbErr>>h(20,-4,4)");' > cmds.txt
    echo 'h->Fit("gaus");' >> cmds.txt
    echo 'c1->Print("bias.png");' >> cmds.txt
    root -l -b fitDiagnostics.root < cmds.txt > bias_tests/${chan}_${label}_mbin${idx}.txt
    mv bias.png bias_tests/plot_${chan}_${label}_mbin${idx}.png
    rm higgsCombineTest.GenerateOnly.mH120.123456.root
    rm higgsCombineTest.FitDiagnostics.mH120.123456.root
    #rm fitDiagnostics.root
    mv fitDiagnostics.root bias_tests/fit_${chan}_${label}_mbin${idx}.root


done



