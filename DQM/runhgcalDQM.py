#!/usr/bin/env python

from __future__ import print_function, division

import os, re, time, subprocess, glob

cmsswSource = "/home/tmudholk/HGCal_DQM/CMSSW_8_0_17/src/HGCal"
# dataFolder = "/home/tmudholk/mount/pchgcal01/hgcaldata/data" # Important to NOT have a trailing slash
dataFolder = "/home/tmudholk/HGCal_DQM/testInput" # Important to NOT have a trailing slash
# outputFolder = "/home/tmudholk/mount/pchgcal01/hgcaldata/PromptFeedback/processingOutput" # Important to NOT have a trailing slash
outputFolder = "/home/tmudholk/HGCal_DQM/testOutput" # Important to NOT have a trailing slash
# dqmPlotsRecoFolder = "/hgcaldata/PromptFeedback/DQMPlots_Reco"
# dqmPlotsDigiFolder = "/hgcaldata/PromptFeedback/DQMPlots_Digi"
# dqmPlotsFolder = "/var/www/html/dqm"
# dqmPlotsFolder = "/home/tmudholk/mount/pchgcal01/var/www/html/dqm"
dqmPlotsFolder = "/home/tmudholk/HGCal_DQM/testDQMPlots"
# commonPrefix = "Ped_Output"
chainSequence = 3
nSpills = 1
pedestalsHighGain = "%s/CondObjects/data/Ped_HighGain_OneLayer_H2CERN.txt"%(cmsswSource)
pedestalsLowGain = "%s/CondObjects/data/Ped_LowGain_OneLayer_H2CERN.txt"%(cmsswSource)
pathToProcessingStatusLogger = "%s/hgcalDQMProcessingStatusLog"%(outputFolder)
pathToScriptLogger = "%s/runhgcalDQMLog"%(outputFolder)
listOfRunsAlreadyProcessed = []
latestListOfRunsInDataFolder = []
latestTypesOfRunsInDataFolder = {} # empty dictionary
listOfRunsToProcess = []
sleepTime = 6

def initiateListOfRunsAlreadyProcessed():
    global listOfRunsAlreadyProcessed
    if (os.path.isfile(pathToProcessingStatusLogger)):
        processingStatusLoggerFile = open(pathToProcessingStatusLogger,'r')
        for runNumberProcessedStr in processingStatusLoggerFile:
            runNumberProcessed = int(runNumberProcessedStr)
            listOfRunsAlreadyProcessed += [runNumberProcessed]
        processingStatusLoggerFile.close()
    else:
        os.system("touch %s"%(pathToProcessingStatusLogger))

def updateListOfRunsInDataFolder():
    global latestListOfRunsInDataFolder, latestTypesOfRunsInDataFolder
    # listOfAllFilesInDataFolder = os.listdir(dataFolder)
    listOfdaqdoneFiles = glob.glob(dataFolder+"/.daqdone.*")
    runNumberMatchPattern = r"^\.daqdone\.([0-9]*)$"
    compiledRunNumberMatchPattern = re.compile(runNumberMatchPattern)
    for fullFilePath in listOfdaqdoneFiles:
        filename = fullFilePath[1+len(dataFolder):]
        runNumberMatchAttempt = compiledRunNumberMatchPattern.match(filename)
        if runNumberMatchAttempt:
            runNumberFoundStr = runNumberMatchAttempt.group(1)
            runNumberFound = int(runNumberFoundStr)
            if not(runNumberFound in latestListOfRunsInDataFolder):
                latestListOfRunsInDataFolder += [runNumberFound]
                if (os.path.isfile("%s/HGCRun_Output_%06d.txt"%(dataFolder, runNumberFound))):
                    latestTypesOfRunsInDataFolder[runNumberFound] = "HGCRun"
                elif (os.path.isfile("%s/PED_Output_%06d.txt"%(dataFolder, runNumberFound))):
                    latestTypesOfRunsInDataFolder[runNumberFound] = "PED"
                else:
                    latestTypesOfRunsInDataFolder[runNumberFound] = "Unknown"

def getListOfRunsToProcess():
    global listOfRunsToProcess
    listOfRunsToProcess = []
    updateListOfRunsInDataFolder()
    for finished_run in latestListOfRunsInDataFolder:
        if not(finished_run in listOfRunsAlreadyProcessed):
            listOfRunsToProcess += [int(finished_run)]

if __name__ == "__main__":
    initiateListOfRunsAlreadyProcessed()
    print ("listOfRunsAlreadyProcessed:%s"%(listOfRunsAlreadyProcessed))
    while True:
        getListOfRunsToProcess()
        print ("List of runs in DQM queue:%s"%(listOfRunsToProcess))
        processingStatusLoggerFile = open(pathToProcessingStatusLogger,'a')
        for runToProcess in listOfRunsToProcess:
            runTypeToProcess = latestTypesOfRunsInDataFolder[runToProcess]
            if (runTypeToProcess == "PED"):
                chainSequence = 1
                nSpills = 1
                pedestalsHighGain = "%s/Ped_HighGain_OneLayer_%06d.txt"%(outputFolder, runToProcess)
                pedestalsLowGain = "%s/Ped_LowGain_OneLayer_%06d.txt"%(outputFolder, runToProcess)
            elif (runTypeToProcess == "HGCRun"):
                chainSequence = 5
                nSpills = 1
                pedestalsHighGain = "%s/CondObjects/data/Ped_HighGain_L8.txt"%(cmsswSource)
                pedestalsLowGain = "%s/CondObjects/data/Ped_LowGain_L8.txt"%(cmsswSource)

            if (runTypeToProcess == "PED" or runTypeToProcess == "HGCRun"):
                cmsswArguments = "print dataFolder=%s outputFolder=%s runNumber=%d runType=%s chainSequence=%d nSpills=%d pedestalsHighGain=%s pedestalsLowGain=%s"%(dataFolder, outputFolder, runToProcess, runTypeToProcess, chainSequence, nSpills, pedestalsHighGain, pedestalsLowGain)
                print ("For run number %d, runTypeToProcess = %s and cmsswArguments = %s"%(runToProcess, runTypeToProcess, cmsswArguments))
                print ("About to execute: cmsRun test_cfg.py %s"%(cmsswArguments))
                subprocess.call("cd %s && eval `scram runtime -sh` && cmsRun test_cfg.py %s && cd -"%(cmsswSource, cmsswArguments), shell=True)
                dqmOutputFolder = dqmPlotsFolder + "/%d"%(runToProcess)
                if (not(os.path.isdir(dqmOutputFolder))):
                    os.system("mkdir -p "+dqmOutputFolder)
                if (not(os.path.isdir(dqmOutputFolder+"/Detailed"))):
                    os.system("mkdir -p "+dqmOutputFolder+"/Detailed")
                if (runTypeToProcess == "HGCRun"):
                    subprocess.call("cd %s && eval `scram runtime -sh` && cd - && root -b -q \"DumpPlotsReco.C++(\\\"%s/HGCRun_Output_%06d_Reco.root\\\", \\\"%s\\\", %d, %d)\""%(cmsswSource, outputFolder, runToProcess, dqmOutputFolder, runToProcess, nSpills), shell=True)
                elif (runTypeToProcess == "PED"):
                    subprocess.call("cd %s && eval `scram runtime -sh` && cd - && root -b -q \"DumpPlotsDigi.C++(\\\"%s/PED_Output_%06d_Digi.root\\\", \\\"%s\\\", %d, %d)\""%(cmsswSource, outputFolder, runToProcess, dqmOutputFolder, runToProcess, nSpills), shell=True)
            else:
                print ("Only runtypes PED and HGCRun supported for now")
            listOfRunsAlreadyProcessed += [runToProcess]
            processingStatusLoggerFile.write("%d\n"%(runToProcess))
        processingStatusLoggerFile.close()
            
        time.sleep(sleepTime)
