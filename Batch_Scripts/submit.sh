#!/bin/bash

cd /afs/cern.ch/work/j/jcoulter/WORK/CMSSW_5_3_20/src/rootfiles/
#cmsenv
eval `scramv1 runtime -sh`

#Added by Ian
export X509_USER_PROXY=~/x509_user_proxy/proxy
voms-proxy-init --noregen
#</Ian>

cd /afs/cern.ch/work/j/jcoulter/WORK/CMSSW_5_3_20/src/rootfiles/

echo "root -l -b -q thrust_HiForest.C++"
echo "First = $FIRST and last file = $LAST"   

root -b > thrust_data${FIRST}-${LAST}.log <<EOF
.x thrust_HiForest.C(${FIRST},${LAST},${JOBNUM});
.q
EOF

echo "Done all jobs!"
