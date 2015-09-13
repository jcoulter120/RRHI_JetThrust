'''
Ian Hunt-Isaak
ihuntisa@oberlin.edu -Primary email un may 2017
ianhuntisaak@gmail.com - Permanent email


ASSUMES THAT FILES ARE NAMED as Pythia_[minFile-maxFile]_[numEventString].log   -  change input section below to adjust log names to your log files.
This file extracts the sum of cross section weights from the below section of the Pythia output log file. SUM = sum of all cross sections. 


Searches for "PYMAXI" in a line and then gets the cross section from that which is SUMString
Also searches for "PYSTAT" in a line then gets the cross from that which is SUM_STATString


PYMAXI:

 ******** PYMAXI: summary of differential cross-section maximum search ********

           ==========================================================
           I                                      I                 I
           I  ISUB  Subprocess name               I  Maximum value  I
           I                                      I                 I
           ==========================================================
           I                                      I                 I
           I   11   f + f' -> f + f' (QCD)        I    5.1634D-02   I
           I   12   f + fbar -> f' + fbar'        I    5.1324D-04   I
           I   13   f + fbar -> g + g             I    4.5178D-04   I
           I   28   f + g -> f + g                I    4.1964D-01   I
           I   53   g + g -> f + fbar             I    1.0936D-02   I
           I   68   g + g -> g + g                I    4.2013D-01   I
           I   96   Semihard QCD 2 -> 2           I    5.5054D+03   I
           I                                      I                 I
           ==========================================================

PYSTAT:
********* PYSTAT:  Statistics on Number of Events and Cross-sections *********

 ==============================================================================
 I                                  I                            I            I
 I            Subprocess            I      Number of points      I    Sigma   I
 I                                  I                            I            I
 I----------------------------------I----------------------------I    (mb)    I
 I                                  I                            I            I
 I N:o Type                         I    Generated         Tried I            I
 I                                  I                            I            I
 ==============================================================================
 I                                  I                            I            I
 I   0 All included subprocesses    I        50000        237080 I  1.925D-01 I
 I  11 f + f' -> f + f' (QCD)       I         2661         13643 I  9.931D-03 I
 I  12 f + fbar -> f' + fbar'       I           43           141 I  1.523D-04 I
 I  13 f + fbar -> g + g            I           46           109 I  1.451D-04 I
 I  28 f + g -> f + g               I        20155        110478 I  7.705D-02 I
 I  53 g + g -> f + fbar            I          990          2967 I  3.689D-03 I
 I  68 g + g -> g + g               I        26105        109742 I  1.015D-01 I
 I                                  I                            I            I
 ==============================================================================

'''
import sys
#======================INPUTS===================
minFile=0
maxFile=12
numEventString="12"
#=======================end inputs
skipSix=False
SUMString=""
SUM_STATString=""

for fileNum in range(minFile,maxFile):
    SUM=0
    SUM_STAT=0

    print "\n\n\nFILE   :"+str(fileNum)+""
    File = open('cd ../GeneratorInterface/Pythia_'+str(fileNum)+'_numEvents'+numEventString+'/Pythia_'+str(fileNum)+'_numEvents'+numEventString+'.log')
    #f = open('out'+str(fileNum)+'.txt','w')
    print 'Pythia_Sum_'+str(fileNum)+'_'+numEventString+'.log'
    for i in range(0,4000000):
        line= File.readline()
        if "PYMAXI" in line:
            for z in range(0,7):
                File.readline()
            for x in range(0,7):
                line=File.readline().replace('D','E')
      #          sys.stdout.write(line)
                #print(float(line[55:-2]))
                SUM += float(line[55:-2])
                #print SUM
            SUMString+=str(SUM)+","
        if "PYSTAT" in line:
            print "found PYSTAT"
            #trashes first few useless lines
            for c in range(0,11):
                File.readline()
                #print line
            for z in range(11,18):
                #print (File.readline())
                #File.readline()
                line=File.readline().replace('D','E')
                print line[68:-2]
            #if line[68:-2] in ['\n','\r\n',"***","==="]:
            #if line is whitespace, continue
            #if line[68:-2] and not line[68:-2].isspace():
                if not line[68:-2]:
                    continue
                if '*' or '=' or ' ' or '0*' in line[68:-2]:
                    continue
            #  sys.stdout.write(line)
                print line[68:-2]
                print float(line[68:-2])
                print str(float(line[68:-2]))
                SUM_STAT+=float(line[68:-2])
            SUM_STATString+=str(float(line[68:-2]))+","

            
print "SUM cross sections:"
print SUMString[:-1]#slice removes final comma

print "SUM STAT cross sections:"
print SUM_STATString[:-1]#slice removes final comma
