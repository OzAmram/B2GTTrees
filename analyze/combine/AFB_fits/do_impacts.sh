#!/bin/bash

set -x 
#idx=0
for idx in 0 1 2 3 4 5
do

    workspace=workspaces/combined_fit${idx}.root
    pars="BTAG,alpha,elScale,elSmear,elHLT,elID,elRECO,muRC,muID,muHLT,muTRK,pdf,RENORM,FAC,lumi,dy_xsec,bk_xsec,ee_qcd,mu_qcd,emu_qcd,Pu,Rdy_ee_ss"

    combineTool.py -M Impacts -d $workspace -m 125 --doInitialFit --robustFit 1
    combineTool.py -M Impacts -d $workspace -m 125 --doFits --robustFit 1 --named $pars
    combineTool.py -M Impacts -d $workspace -m 125 -o impacts/mbin${idx}.json --named $pars
    plotImpacts.py -i impacts/mbin${idx}.json -o impacts/impact_plot_mbin${idx}
    rm higgsCombine_initialFit*
    rm higgsCombine_paramFit_Test_*
done



