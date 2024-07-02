#!/bin/bash

export USER="$(id -u -n)"
export LOGNAME=${USER}
export HOME=/sphenix/u/${LOGNAME}/sphnx_software/analysis/JS-Jet/JetValidation/macro/
export MYINSTALL=/sphenix/u/jlnliu/install/

source /opt/sphenix/core/bin/sphenix_setup.sh -n
source /opt/sphenix/core/bin/setup_local.sh $MYINSTALL
 
 # print the environment - needed for debugging
printenv

#echo $dataFileList
#simFileList=/sphenix/user/bkimelman/analysis/CaloEmbedAnalysis/macro/fileLists/simLists/simDST_$(printf "%04d" $1).list
caloList=/sphenix/u/jlnliu/sphnx_software/analysis/JS-Jet/JetValidation/macro/dst_calo_cluster.list
globalList=/sphenix/u/jlnliu/sphnx_software/analysis/JS-Jet/JetValidation/macro/dst_global.list
truthList=/sphenix/u/jlnliu/sphnx_software/analysis/JS-Jet/JetValidation/macro/dst_truth_jet.list
#dataFileList=/sphenix/user/bkimelman/Run23745_ana399_DSTs/Run23745_ana399_DST_$(printf "%04d" $1).list

#echo simFileList: $simFileList
echo caloList: $caloList
echo globalList: $globalList
echo truthList: $truthList
#echo dataFileList: $dataFileList

 # use this for sim
root.exe -q -b Fun4All_JetVal.C\(\"$truthList\",\"$caloList\",\"$globalList\",\"/sphenix/u/jlnliu/sphnx_software/analysis/JS-Jet/JetValidation/macro/output.root\"\)

 
 # use this for data

#root.exe -q -b Fun4All_JetVal.C\(\"\",\"$dataFileList\",\"\",\"/sphenix/u/jlnliu/analysis/JS-Jet/JetValidation/macro/Run23745_ana399_DST_$(printf "%04d" $1).root\"\)
 
echo all done
