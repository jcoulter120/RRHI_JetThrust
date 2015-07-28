
#!/bin/bash
cd /afs/cern.ch/work/j/jcoulter/WORK/CMSSW_5_3_20/src/rootfiles/Batch_Scripts
#cmsenv
eval `scramv1 runtime -sh`

cd /afs/cern.ch/work/j/jcoulter/WORK/CMSSW_5_3_20/src/rootfiles/Batch_Scripts
voms-proxy-init --voms cms

which root

echo $ROOTSYS

nJobs=1
i=0
while [ $i -lt $nJobs ];
do 
   let "start=i*1"
   let "end=(i+1)*1"
  export FIRST=$start  
  export LAST=$end
  export JOBNUM=$i
  echo "First = $FIRST and last file = $LAST"   
  bsub -R "pool>5000" -M 3000000 -q 1nd -J merge_job_${i} < /afs/cern.ch/work/j/jcoulter/WORK/CMSSW_5_3_20/src/rootfiles/Batch_Scripts/submit.sh
  let "i++"
done

echo "submit all jobs!"
