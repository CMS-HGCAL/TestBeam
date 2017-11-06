#!/usr/bin/env python

from __future__ import print_function, division

from math import pow

dEdX_for = {'W': 2.21, 'Cu': 1.257, 'Air': 0, 'Si': 0.3876, 'AcSi': 0.3876, 'Pb': 1.274}
#dEdX_for = {'W': 2.21, 'Cu': 1.257, 'Air': 0, 'Si': 0., 'AcSi': 0.3876, 'Pb': 1.274}
dEdX_for['WCu'] = 0.75*dEdX_for['W'] + 0.25*dEdX_for['Cu']

X0_for = {'W': 3.504, 'Cu': 14.36, 'Air': 303900, 'Si': 93.70, 'AcSi': 93.70, 'Pb': 5.612}
X0_for['WCu'] = pow(0.75*pow(X0_for['W'],-1) + 0.25*pow(X0_for['Cu'],-1),-1)

#path_to_setup_data = "./TM_setup_data_FNAL.txt"
#path_to_setup_data = "./SJ_setup_data_CERN.txt"
#path_to_setup_data = "./setup_data_CERN_setup8Layers_centerShower.txt"
#path_to_setup_data = "./setup_data_CERN_setup28Layers.txt"

##28 layers
#path_to_setup_data = "./setup_data_test.txt"
##8layers central shower
#path_to_setup_data = "./GeometrySetup/setup_data_CERN_setup8Layers_centerShower.txt"

##8layers tail shower
#path_to_setup_data = "./GeometrySetup/setup_data_CERN_setup8Layers_tailShower.txt"
 
#path_to_setup_data = "./GeometrySetup/setup_data_FNAL.txt"

## 28 layers FNAL
path_to_setup_data = "./GeometrySetup/setup_data_FNAL_28layers.txt"

dEdX_times_deltaX_before_first_layer = 0
radiation_lengths_before_first_layer = 0 # this means 0 radiation lengths before the beam encounters the first sublayer

def getAveragedWeights(weights_raw):
    averagedWeights = []
    for layerCounter in range(0,len(weights_raw)-1): # first to (n-1)th layers, where n is the total number of layers
        averagedWeights += [0.5*(weights_raw[layerCounter]+weights_raw[layerCounter+1])]
    averagedWeights += [weights_raw[-1]] # last layer
    return averagedWeights

def getNormalized(weights):
    normalizationFactor = len(weights)/sum(weights)
    normalizedWeights = []
    for weight in weights:
        normalizedWeights.append(weight*normalizationFactor)
    return normalizedWeights

def printGapWithText(text):
    print ("\n"+"-"*100)
    print ("%s"%(text))

def getPrettierArray(weights):
    prettierArray = []
    for weight in weights:
        prettierArray += [float("%.3f"%(weight))]
    return prettierArray

printGapWithText("dEdXs for materials: %s"%(dEdX_for))
printGapWithText("X0s for materials: %s"%(X0_for))

dEdXs = []
X0s = []

setup_data = open(path_to_setup_data)
lastLayer = ""
dEdX_layer = 0
X0_layer = 0

# Initialize with some number of radiation lengths/dEdX times deltaX:
dEdX_layer = dEdX_times_deltaX_before_first_layer
X0_layer = radiation_lengths_before_first_layer 

for sublayer_data in setup_data:
    sublayer_data_array = (sublayer_data.strip()).split(',')
    sublayer = [float(sublayer_data_array[0]), sublayer_data_array[1]]
    print ("Sublayer: %s"%(sublayer))
    if (sublayer[1] == "AcSi"):
        #print (" AAAA sublayer[0]*dEdX_for[sublayer[1]] = %s"%(sublayer[0]))
        dEdXs.append(dEdX_layer)
        X0s.append(X0_layer)
        dEdX_layer = 0
        X0_layer = 0
    if (sublayer[1] != "AcSi"):
        dEdX_layer += sublayer[0]*dEdX_for[sublayer[1]]
        #print (" #### sublayer[0]*dEdX_for[sublayer[1]] = %s"%(sublayer[0]*dEdX_for[sublayer[1]]))
        X0_layer += sublayer[0]/X0_for[sublayer[1]]
        lastLayer = sublayer[1]

dEdXs.append(dEdX_layer)
X0s.append(X0_layer)

printGapWithText("Checking... sum_over_sublayers(dE/dX(sublayer)*thickness(sublayer)) = %s"%(getPrettierArray(dEdXs)))
printGapWithText("Checking... sum_over_sublayers(thickness(sublayer)/X0(sublayer)) = %s"%(getPrettierArray(X0s)))

printGapWithText("Checking... total number of radiation lengths: %f"%(sum(X0s)))

dEdX_normalization_factor = len(dEdXs)/sum(dEdXs)
X0_normalization_factor = len(X0s)/sum(X0s)

dEdXs_normalized = getNormalized(dEdXs)
X0s_normalized = getNormalized(X0s)

dEdXs_averaged_nn = getAveragedWeights(dEdXs)
X0s_averaged_nn = getAveragedWeights(X0s)

dEdXs_averaged = getNormalized(getAveragedWeights(dEdXs))
X0s_averaged = getNormalized(getAveragedWeights(X0s))

printGapWithText("Checking... dEdX weights, unnormalized: %s"%(getPrettierArray(dEdXs)))
printGapWithText("Checking... X0 weights, unnormalized: %s"%(getPrettierArray(X0s)))

printGapWithText("Checking... sum of dEdX weights after normalization: %f"%(sum(dEdXs_normalized)))
printGapWithText("Checking... sum of X0 weights after normalization: %f"%(sum(X0s_normalized)))
printGapWithText("Checking... sum of dEdX averaged weights after normalization: %f"%(sum(dEdXs_averaged)))
printGapWithText("Checking... sum of X0 averaged weights after normalization: %f"%(sum(X0s_averaged)))

printGapWithText("Normalized dEdX weights: %s"%(getPrettierArray(dEdXs_normalized)))
printGapWithText("Normalized X0 weights: %s"%(getPrettierArray(X0s_normalized)))

printGapWithText(" >>>> Non Normalized averaged dEdX weights: %s"%(getPrettierArray(dEdXs_averaged_nn)))
printGapWithText("Non Normalized averaged X0 weights: %s"%(getPrettierArray(X0s_averaged_nn)))

printGapWithText("Normalized averaged dEdX weights: %s"%(getPrettierArray(dEdXs_averaged)))
printGapWithText("Normalized averaged X0 weights: %s"%(getPrettierArray(X0s_averaged)))
printGapWithText("")
