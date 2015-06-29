#!/bin/bash

#moves to directory where txt are located
cd /afs/cern.ch/work/j/jcoulter/WORK/CMSSW_7_4_0/src/rootfiles/

#assures that you don't forget to get a proxy
voms-proxy-init --voms cms

#echo > pbpb_Data.txt

#runs through all the lines in the txt
while read -r line 
do
    #chops off the file name from the end of the file address, then saves to a var
    STR=`echo "$line" | awk -F"/" '{print $NF}'`
    #copies the file to the indicated directory by providing the address and the root file name
    lcg-cp -v srm://se01.cmsaf.mit.edu:8443///"$line" "$STR"
    #sends the file to CERNBox
    xrdcp "$STR" root://eosuser.cern.ch://eos/user/j/jcoulter/Data/pp_MC/"$STR"
    #Removes created root file from current directory
    rm "$STR"

done < "PP_MC_ntuples_MIT.txt" #finishes the while loop and indicates which file is to be read
