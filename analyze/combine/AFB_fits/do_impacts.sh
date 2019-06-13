#!/bin/bash

set -x 
#idx=0
for idx in 0 1 2 3 4 5
do

    workspace=workspaces/combined_fit${idx}.root
    pars="Pu,BTAG,elScaleStat,elScaleSyst,elScaleGain,elSmear,elHLT,elID,elRECO,muRC,muID,muHLT,alpha,alphaS,alphaDen,RENORM,FAC,lumi,dy_xsec,bk_xsec,gam_xsec,ee_qcd,mu_qcd,emu_qcd,Rdy_ee_ss,R_ee_os_fakes"

    combineTool.py -M Impacts -d $workspace -m 125 --doInitialFit --robustFit 1
    combineTool.py -M Impacts -d $workspace -m 125 --doFits --robustFit 1 --named $pars --parallel 6
    combineTool.py -M Impacts -d $workspace -m 125 -o impacts/mbin${idx}.json --named $pars
    plotImpacts.py -i impacts/mbin${idx}.json -o impacts/impact_plot_mbin${idx} -t impacts/rename.json
    rm higgsCombine_initialFit*
    rm higgsCombine_paramFit_Test_*
done



