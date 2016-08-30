#!/usr/bin/env python

from __future__ import print_function, division

import os, re, time, subprocess, glob

cmsswSource = "/hgcaldata/PromptFeedback/CMSSW_8_0_17/src/HGCal"
dataFolder = "/hgcaldata/data" # Important to NOT have a trailing slash
outputFolder = "/hgcaldata/PromptFeedback/testOutput" # Important to NOT have a trailing slash
# dqmPlotsRecoFolder = "/hgcaldata/PromptFeedback/DQMPlots_Reco"
# dqmPlotsDigiFolder = "/hgcaldata/PromptFeedback/DQMPlots_Digi"
dqmPlotsFolder = "/var/www/html/dqm"
# commonPrefix = "Ped_Output"
chainSequence = 3
nSpills = 6
pedestalsHighGain = "%s/CondObjects/data/Ped_HighGain_OneLayer_H2CERN.txt"%(cmsswSource)
pedestalsLowGain = "%s/CondObjects/data/Ped_LowGain_OneLayer_H2CERN.txt"%(cmsswSource)
pathToProcessingStatusLogger = "%s/hgcalDQMProcessingStatusLog"%(outputFolder)
pathToScriptLogger = "%s/runhgcalDQMLog"%(outputFolder)
listOfRunsAlreadyProcessed = []
latestListOfRunsInDataFolder = []
latestTypesOfRunsInDataFolder = {} # empty dictionary
listOfRunsToProcess = []
sleepTime = 60

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
        # DQMHelperOptions = "%s %s %s %s"%(cmsswSource, dataFolder, outputFolder, commonPrefix)
        for runToProcess in listOfRunsToProcess:
            runTypeToProcess = latestTypesOfRunsInDataFolder[runToProcess]
            # subprocess.call("bash runhgcalDQMHelper.sh %s %d %s"%(DQMHelperOptions, runToProcess, pathToScriptLogger), shell=True)
            if (runTypeToProcess == "PED"):
                chainSequence = 1
                nSpills = 6
                pedestalsHighGain = "%s/Ped_HighGain_OneLayer_%06d.txt"%(outputFolder, runToProcess)
                pedestalsLowGain = "%s/Ped_LowGain_OneLayer_%06d.txt"%(outputFolder, runToProcess)
            elif (runTypeToProcess == "HGCRun"):
                chainSequence = 3
                nSpills = 34
                pedestalsHighGain = "%s/CondObjects/data/Ped_HighGain_OneLayer_H2CERN.txt"%(cmsswSource)
                pedestalsLowGain = "%s/CondObjects/data/Ped_LowGain_OneLayer_H2CERN.txt"%(cmsswSource)

            if (runTypeToProcess == "PED" or runTypeToProcess == "HGCRun"):
                cmsswArguments = "print dataFolder=%s outputFolder=%s runNumber=%d runType=%s chainSequence=%d nSpills=%d pedestalsHighGain=%s pedestalsLowGain=%s"%(dataFolder, outputFolder, runToProcess, runTypeToProcess, chainSequence, nSpills, pedestalsHighGain, pedestalsLowGain)
                print ("For run number %d, runTypeToProcess = %s and cmsswArguments = %s"%(runToProcess, runTypeToProcess, cmsswArguments))
                # subprocess.call("bash runhgcalDQMHelper.sh %s %d %s %d %s %s"%(DQMHelperOptions, runToProcess, runTypeToProcess, chainSequence, pedestalsHighGain, pedestalsLowGain), shell=True)
                # cmsRunCommand = "bash runhgcalDQMHelper.sh %s %s"%(cmsswSource, cmsswArguments)
                # print ("About to execute %s"%(executeCommand))
                subprocess.call("cd %s && eval `scram runtime -sh` && cmsRun test_cfg.py %s && cd -"%(cmsswSource, cmsswArguments), shell=True)
                if (runTypeToProcess == "HGCRun"):
                    subprocess.call("cd %s && eval `scram runtime -sh` && cd - && root -b -q \"DumpPlotsReco.C++(\\\"%s/HGCRun_Output_%06d_Reco.root\\\", \\\"%s\\\", %d)\""%(cmsswSource, outputFolder, runToProcess, dqmPlotsFolder, runToProcess), shell=True)
                elif (runTypeToProcess == "PED"):
                    subprocess.call("cd %s && eval `scram runtime -sh` && cd - && root -b -q \"DumpPlotsDigi.C++(\\\"%s/PED_Output_%06d_Digi.root\\\", \\\"%s\\\", %d)\""%(cmsswSource, outputFolder, runToProcess, dqmPlotsFolder, runToProcess), shell=True)
            else:
                print ("Only runtypes PED and HGCRun supported for now")
            listOfRunsAlreadyProcessed += [runToProcess]
            processingStatusLoggerFile.write("%d\n"%(runToProcess))
        processingStatusLoggerFile.close()
            
        time.sleep(sleepTime)
