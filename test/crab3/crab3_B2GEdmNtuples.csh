#!/bin/tcsh
# Author: Janos Karancsi
# E-mail: janos.karancsi@cern.ch
# Tool: multicrab-like tool to create, submit, get status of crab3 tasks for B2GEdmNtuples
#       When jobs finish, we can also get report and download all output files locally

echo "Usage:"                                                    >! Usage.txt
echo "  source crab3_B2GEdmNtuples.csh <cmd> <TASKNAME> ..."     >> Usage.txt
echo "  there is a default safety mechanism for each command"    >> Usage.txt
echo "  add --run at the end to excecute them"                   >> Usage.txt
echo ""                                                          >> Usage.txt
echo "Commands:"                                                 >> Usage.txt
echo "1) create <TASKNAME> <PUBLISH_NAME> <DATASETS.txt>"        >> Usage.txt
echo "  <SE_SITE> <SE_USERDIR>"                                  >> Usage.txt
echo "  The txt file should contain a short name and the"        >> Usage.txt
echo "  MINIAOD(SIM) dataset on each line"                       >> Usage.txt
echo "  PUBLISH_NAME should specify the B2GAnaFW version tag"    >> Usage.txt
echo "  SE_SITE and SE_USERDIR should be the location you want"  >> Usage.txt
echo "  to send the output to, eg:"                              >> Usage.txt
echo "  T2_HU_Budapest"                                          >> Usage.txt
echo "  /store/user/jkarancs/SusyAnalysis/B2GEdmNtuple"          >> Usage.txt
echo ""                                                          >> Usage.txt
echo "2) submit <TASKNAME>"                                      >> Usage.txt
echo "  Submit tasks to the grid (check cfg files first)"        >> Usage.txt
echo ""                                                          >> Usage.txt
echo "3) status <TASKNAME>"                                      >> Usage.txt
echo "  Check status of tasks"		                         >> Usage.txt
echo "" 						         >> Usage.txt
echo "4) report <TASKNAME>"                                      >> Usage.txt
echo "  Show commands to get report of tasks"		         >> Usage.txt
echo "" 						         >> Usage.txt
echo "5) get_datasets <TASKNAME>"                                >> Usage.txt
echo "  Get a list of produced datasets"		         >> Usage.txt
echo "" 						         >> Usage.txt
echo "6) make_ttrees <TASKNAME>"                                 >> Usage.txt
echo "  Generate a table in general cern twiki format that"      >> Usage.txt
echo "  contains useful infos for the produced datasets. Try it:">> Usage.txt
echo "  https://twiki.cern.ch/twiki/bin/view/Main/<username>"    >> Usage.txt
echo "" 						         >> Usage.txt
echo "7) download <TASKNAME> <DLDIR>"                            >> Usage.txt
echo "  Download all root output files from the Storage"         >> Usage.txt
echo "  element to a local directory DLDIR"                      >> Usage.txt
echo "  Check that the SRM site name and url is defined for"     >> Usage.txt
echo "  SE_SITE in the se_util.csh script (look for SITE_INFO)"  >> Usage.txt
echo "" 						         >> Usage.txt
echo "8) make_ttrees <TASKNAME> <DIR> <NPARALLEL>"               >> Usage.txt
echo "  Produce locally TTreeNtuples from EdmNtuples"            >> Usage.txt
echo "  Files need to be downloaded with the download command"   >> Usage.txt
echo "  Multiple background cmsRuns can be run in parallel"      >> Usage.txt
echo "  Useful if you have small samples (Few GB)"               >> Usage.txt

if ( $1 == "-h" || $1 == "--help" || $#argv < 2 ) then
    cat Usage.txt; rm Usage.txt; exit
endif

# options
set cmd=$1
set TASKNAME=$2
# Working directory name (created in the current directory)
set TASKDIR="B2G_edm_"$TASKNAME

# Aliases
if ( ! (-e $PWD/source_parallel.csh) || ! (-e $PWD/se_util.csh) ) then
    echo "Please run this script from the same directory where\nsource_parallel.csh and se_util.csh is or copy them here" 
    exit
endif
# script that takes another script as argument and runs n lines in parallel
# used by se_util.csh
alias par_source 'source source_parallel.csh \!*'
# Storage element utility (ls, cp, dl, mkdir etc commands)
alias se         'source se_util.csh \!*'

set rest_args=""
foreach n ( `seq 3 $#argv` )
    set rest_args="$rest_args "$argv[$n]
end
if ( `echo "$rest_args" | grep "\-\-run" | wc -l` ) then
    alias eval_or_echo   'echo \!* >! temp.csh; source temp.csh; rm temp.csh'
    alias source_or_peek 'echo "> Add --run option or Run this script:\nsource "\!*"\n\n> head -5 "\!*":"; head -5 \!*; if ( `cat \!* | wc -l` > 10 ) echo "...\n> tail -5 "\!*":"; if ( `cat \!* | wc -l` > 10 ) tail -5 \!*'
    set dry="0"
else
    alias eval_or_echo   'echo \!*'
    set dry="1"
endif


if ( `echo $cmd | grep "create" | wc -l` ) then
    if ( $#argv < 6 ) then
	cat Usage.txt; rm Usage.txt; exit
    endif
    set PUBNAME=$3
    set TXT_FILE=$4
    set SE_SITE=$5
    set SE_USERDIR=$6
    mkdir $TASKDIR
    echo "SE_SITE "$SE_SITE >! $TASKDIR/config.txt
    echo "SE_USERDIR "$SE_USERDIR >> $TASKDIR/config.txt
    if ( !(-f $TXT_FILE) ) then
	echo "$TXT_FILE doesn't exist"; exit
    endif
    grep "/MINIAOD" $TXT_FILE >! $TASKDIR/input_datasets.txt
    sed "s;TASKDIR;$TASKDIR;;s;SE_SITE;$SE_SITE;;s;SE_USERDIR;$SE_USERDIR;" crab_template_edmntuple_Data_py.txt > $TASKDIR/crab_template_edmntuple_MC_py.txt
    sed "s;TASKDIR;$TASKDIR;;s;SE_SITE;$SE_SITE;;s;SE_USERDIR;$SE_USERDIR;" crab_template_edmntuple_MC_py.txt > $TASKDIR/crab_template_edmntuple_Data_py.txt
    set N=`cat $TASKDIR/input_datasets.txt | wc -l`
    foreach i ( `seq 1 $N` )
	set line=`sed -n "$i"p $TASKDIR/input_datasets.txt`
	set SHORT=`echo $line | awk '{ print $1 }'`
	set DATASET=`echo $line | awk '{ print $2 }'`
	set isData=`echo $DATASET | grep 'MINIAOD$' | wc -l`
	set PUBNAME2=`echo $DATASET | sed "s;/; ;g" | awk '{ print $2 }'`
	if ( `echo $PUBNAME2 | grep "Run2015B-17Jul2015" | wc -l ` ) then
	    set GLOBALTAG="74X_dataRun2_Prompt_v0"
	    set DATAPROC="ReReco17Jul50ns"
	    set JSON="/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_246908-251883_13TeV_PromptReco_Collisions15_JSON_v2.txt"
	    set RUNS="251162-251562"
	else if ( `echo $PUBNAME2 | grep "Run2015B-PromptReco" | wc -l ` ) then
	    set GLOBALTAG="74X_dataRun2_Prompt_v0"
	    set DATAPROC="PromptReco50ns"
	    set JSON="/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_246908-251883_13TeV_PromptReco_Collisions15_JSON_v2.txt"
	    set RUNS="251585-252501"
	else if ( `echo $PUBNAME2 | grep "Asympt50ns" | wc -l ` ) then
	    set GLOBALTAG="MCRUN2_74_V9A"
	    set DATAPROC="MC50ns"
	else	    
	    set GLOBALTAG="MCRUN2_74_V9"
	    set DATAPROC="MC25ns"
	endif
	if ( $isData ) then
	    sed "s;REQNAME;$SHORT;;s;PUBNAME;$PUBNAME"_"$PUBNAME2;;s;DATASET;$DATASET;;s;GLOBALTAG;$GLOBALTAG;;s;DATAPROC;$DATAPROC;;s;JSON;$JSON;" $TASKDIR/crab_template_edmntuple_Data_py.txt > $TASKDIR/crab_$SHORT.py
	else
	    sed "s;REQNAME;$SHORT;;s;PUBNAME;$PUBNAME"_"$PUBNAME2;;s;DATASET;$DATASET;;s;GLOBALTAG;$GLOBALTAG;;s;DATAPROC;$DATAPROC;" $TASKDIR/crab_template_edmntuple_MC_py.txt > $TASKDIR/crab_$SHORT.py
	endif
    end
    rm $TASKDIR/crab_template_edmntuple_Data_py.txt $TASKDIR/crab_template_edmntuple_MC_py.txt
    echo "Config files ready in directory: "$TASKDIR
    ls -ltr $TASKDIR

else if ( `echo $cmd | grep "submit" | wc -l` ) then
    if ( $dry == "1" ) echo "Add --run after command to excecute following lines:\n"
    foreach cfg_file ( `ls -ltr $TASKDIR/*.py | awk '{ print $NF }'`)
        eval_or_echo "crab submit -c $cfg_file --wait"
    end

else if ( `echo $cmd | grep "status" | wc -l` ) then
    if ( $dry == "1" ) echo "Add --run after command to excecute following lines:\n"
    foreach dir ( `ls -ltrd $TASKDIR/* | grep "^d" | awk '{ print $NF }'`)
        eval_or_echo "crab status -d $dir"
    end

else if ( `echo $cmd | grep "finished" | wc -l` ) then
    if ( $dry == "1" ) echo "Add --run after command to excecute following lines:\n"
    foreach dir ( `ls -ltrd $TASKDIR/* | grep "^d" | awk '{ print $NF }'`)
        eval_or_echo "crab status -d $dir | grep finished"
    end

else if ( `echo $cmd | grep "report" | wc -l` ) then
    if ( $dry == "1" ) echo "Add --run after command to excecute following lines:\n"
    foreach dir ( `ls -ltrd $TASKDIR/* | grep "^d" | awk '{ print $NF }'`)
        eval_or_echo "crab report -d $dir"
    end

else if ( `echo $cmd | grep "get_datasets" | wc -l` ) then
    if ( (-f $TASKDIR/output_datasets.txt) ) then
	if ( `cat $TASKDIR/output_datasets.txt | wc -l` != `cat $TASKDIR/input_datasets.txt | wc -l` ) then
	    rm $TASKDIR/output_datasets.txt
	    if ( $dry == "1" ) echo "Add --run after command to excecute following lines:\n"
		foreach dir ( `ls -ltrd $TASKDIR/* | grep "^d" | awk '{ print $NF }'`)
		eval_or_echo "crab status -d $dir | grep 'Output dataset:' | awk '{ print "'$NF'" }' >>! $TASKDIR/output_datasets.txt"
	    end
	endif
    else
	if ( $dry == "1" ) echo "Add --run after command to excecute following lines:\n"
	foreach dir ( `ls -ltrd $TASKDIR/* | grep "^d" | awk '{ print $NF }'`)
	    eval_or_echo "crab status -d $dir | grep 'Output dataset:' | awk '{ print "'$NF'" }' >>! $TASKDIR/output_datasets.txt"
	end
    endif
    eval_or_echo "cat $TASKDIR/output_datasets.txt"
    set Nchar_max="0"
    foreach short ( `cat $TASKDIR/input_datasets.txt | awk '{ print $1 }'` )
	set Nchar=`echo $short | wc -m`
	if ( $Nchar > $Nchar_max ) set Nchar_max=$Nchar
    end
    set N=`cat $TASKDIR/input_datasets.txt | wc -l`
    if ( -f EdmNtuple_"$TASKNAME"_input.txt ) rm EdmNtuple_"$TASKNAME"_input.txt
    foreach i ( `seq 1 $N` )
	set SHORT=`sed -n "$i"p $TASKDIR/input_datasets.txt | awk '{ print $1 }'`
	set OUT_DATASET=`sed -n "$i"p $TASKDIR/output_datasets.txt`
	printf "%-"$Nchar_max"s %s\n" $SHORT $OUT_DATASET >>! EdmNtuple_"$TASKNAME"_input.txt
    end

else if ( `echo $cmd | grep "getoutput" | wc -l` ) then
    if ( $#argv < 3 ) then
	cat Usage.txt; rm Usage.txt; exit
    endif
    if ( $dry == "1" ) echo "Add --run after command to excecute following lines:\n"
    set DLDIR=`echo $3"/"$TASKNAME | sed "s;//;/;"`
    foreach dir ( `ls -ltrd $TASKDIR/* | grep "^d" | awk '{ print $NF }'`)
	eval_or_echo "mkdir -p $DLDIR/"`echo $dir | sed "s;$TASKDIR/crab_;;"`
        eval_or_echo "crab getoutput -d $dir --outputpath=$DLDIR/"`echo $dir | sed "s;$TASKDIR/crab_;;"`"/"
    end

else if ( `echo $cmd | grep "download" | wc -l` ) then
    if ( $#argv < 3 ) then
	cat Usage.txt; rm Usage.txt; exit
    endif
    set DLDIR=`echo $3"/"$TASKNAME | sed "s;//;/;"`
    set SE_SITE=`grep SE_SITE $TASKDIR/config.txt | awk '{ print $2 }'`
    set SE_USERDIR=`grep SE_USERDIR $TASKDIR/config.txt | awk '{ print $2 }'`
    if ( $dry == "1" ) echo "source dl_$TASKNAME.csh $DLDIR\nOR add --run after command to excecute following lines:\n"
    echo 'if ( $#argv < 1 ) echo "Please specify directory where you want to download files"' >! dl_$TASKNAME.csh
    echo 'alias par_source "source $CMSSW_BASE/src/Analysis/B2GTTrees/test/crab3/source_parallel.csh \\!*"' >> dl_$TASKNAME.csh
    echo 'alias se         "source $CMSSW_BASE/src/Analysis/B2GTTrees/test/crab3/se_util.csh \\!*"\n' >> dl_$TASKNAME.csh
    foreach taskdir ( `ls -ltrd $TASKDIR/*/ | sed "s;/; ;g" | awk '{ print $NF }'`)
        set primary_dataset=`grep inputDataset $TASKDIR/$taskdir.py | sed "s;/; ;g" | awk '{ print $4 }'`
        set pubname=`grep publishDataName $TASKDIR/$taskdir.py | sed "s;'; ;g" | awk '{ print $3 }'`
        set SAMPLEDIR=`echo $taskdir | sed "s;crab_;;"`
        set time=`se ls "$SE_SITE":"$SE_USERDIR/$primary_dataset/$pubname"`
        eval_or_echo "mkdir -p $DLDIR/$SAMPLEDIR"
        eval_or_echo "se dl_mis $SE_SITE":"$SE_USERDIR/$primary_dataset/$pubname/$time/0000 $DLDIR/$SAMPLEDIR --par 4 --run"
        echo "mkdir -p "'$1'"/$SAMPLEDIR" >> dl_$TASKNAME.csh
        echo "se dl_mis $SE_SITE":"$SE_USERDIR/$primary_dataset/$pubname/$time/0000 "'$1'"/$SAMPLEDIR --par 4 --run" >> dl_$TASKNAME.csh
    end
    if ( ( `grep "EDM_NTUPLE" $TASKDIR/config.txt | wc -l` == 1 ) && ( $dry == "0" ) ) echo "EDM_NTUPLE $DLDIR" >> $TASKDIR/config.txt

else if ( `echo $cmd | grep "make_ttrees" | wc -l` ) then
    if ( $#argv < 4 ) then
	cat Usage.txt; rm Usage.txt; exit
    endif
    if ( `grep "EDM_NTUPLE" $TASKDIR/config.txt | wc -l` == 0 ) then
	echo "EdmNtuple not yet downloaded, issue command:\ndownload <TASKNAME> <DLDIR>"
    else
	if ( $dry == "1" ) echo "Add --run after command to excecute following lines:\n"
        set EDM_NTUPLE=`grep EDM_NTUPLE $TASKDIR/config.txt | awk '{ print $2 }'`
	set TTREEDIR=$3
        set Nparallel=$4
        foreach dir ( `ls -ltrd $EDM_NTUPLE/* | grep "^d" | sed "s;/; ;g" | awk '{ print $NF }'` )
            eval_or_echo "mkdir -p $TTREEDIR/$dir"
            if ( -f $TASKDIR/make_ttrees_"$TASKNAME"_$dir.csh ) rm $TASKDIR/make_ttrees_"$TASKNAME"_$dir.csh
            foreach file ( `ls $EDM_NTUPLE/$dir | grep ".root"` )
                set outfile=`echo "$file" | sed "s;B2GEDMNtuple;B2GTTreeNtupleExtra;"`
                echo "nice cmsRun $CMSSW_BASE/src/Analysis/B2GTTrees/test/B2GEdmToTTreeNtupleExtra_cfg.py sample=file:$EDM_NTUPLE/$dir/$file outputLabel=$TTREEDIR/$dir/$outfile" >>! $TASKDIR/make_ttrees_"$TASKNAME"_$dir.csh
            end
            eval_or_echo "source source_parallel.csh $TASKDIR/make_ttrees_"$TASKNAME"_$dir.csh $Nparallel"
        end
    endif

else if ( `echo $cmd | grep "make_twiki" | wc -l` ) then
    set infile=$TASKDIR/input_datasets.txt
    set outfile=$TASKDIR/output_datasets.txt
    echo "|  *Dataset*  |  *B2GEdmNtuple*  |  *Nevents*  |  *Nfile*  |  *LO cross section (pb)*  |"
    foreach in_dataset (`awk '{ print $2 }' $infile`)
	set primary_dataset=`echo $in_dataset | sed "s;/; ;g" | awk '{ print $1 }'`
	set out_dataset=`grep "/$primary_dataset/" $outfile`
	set nevents=`das_client.py --query="dataset=$in_dataset | grep dataset.nevents" | tail -1`
	set nfiles=`das_client.py --query="dataset=$in_dataset | grep dataset.nfiles" | tail -1`
	set parent=`das_client.py --query="parent dataset=$in_dataset" | tail -1`
	set ntry=0
	while ( (`echo $parent | grep '/GEN-SIM$' | wc -l` == 0) && ( $ntry < 10 ) )
	    set ntry=`expr $ntry + 1`
	    set parent2=`das_client.py --query="parent dataset=$parent" | tail -1`
	    if ( $parent2 != "None" ) then
		set parent=$parent2
	    else 
		set ntry=10
	    endif

	end
	set prep_id=`das_client.py --query="dataset=$parent | grep dataset.prep_id" | tail -1`
	#set LO_XSec=`das_client.py --query="mcm prepid=$prep_id | grep mcm.generator_parameters.cross_section" | tail -1`
	set LO_XSec=`das_client.py --format=json --query="mcm prepid=$prep_id | grep mcm.generator_parameters.cross_section" | tr "," "\n" | grep '"cross_section"' | tail -1 | sed "s;};;g;s;];;g" | awk '{ print $NF }'`
	echo -n "| [[https://cmsweb.cern.ch/das/request?input=$in_dataset&instance=prod%2Fglobal]["`echo $in_dataset | cut -d '/' -f-3`"]] "
	echo -n "|  [[https://cmsweb.cern.ch/das/request?input=$out_dataset&instance=prod%2Fphys03][DAS link]]  "
	echo -n "|  $nevents "
	echo -n "|  $nfiles "
	#echo -n "|  $LO_XSec [[https://cmsweb.cern.ch/das/request?view=list&limit=10&instance=prod%2Fglobal&input=mcm+prepid%3D$prep_id+|+grep+mcm.generator_parameters.cross_section][(Link)]] |\n"
	echo -n "|  $LO_XSec [[https://cms-pdmv.cern.ch/mcm/requests?dataset_name=$primary_dataset&page=0&shown=262163][(Link)]] |\n"
    end

endif
rm Usage.txt
