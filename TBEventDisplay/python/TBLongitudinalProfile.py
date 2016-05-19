#!/usr/bin/env python
#-----------------------------------------------------------------------------
# File:        TBLongitudinalProfile.py
# Description: TB 2016 - display of longitudinal shower profile
# Created:     19-May-2016 Harrison B. Prosper
#-----------------------------------------------------------------------------
import sys, os, re
from string import atof, lower, replace, strip, split, joinfields, find
from HGCal.TBEventDisplay.TBUtil import *
from HGCal.TBStandaloneSimulator.TBGeometryUtil import *
from math import *
from ROOT import *
#------------------------------------------------------------------------------
class LongitudinalProfile:

    def __init__(self, parent, page):
        self.cellmap  = parent.cellmap
        self.geometry = parent.geometry
        self.sensitive= parent.sensitive
        self.canvas   = page.canvas
        self.nlayers  = len(self.sensitive)

        nbins = self.nlayers
        xmin  = 1
        xmax  = self.nlayers+1
        name  = 'lprofile'
        h = TH1F(name, "", nbins, xmin, xmax)
        h.SetFillStyle(3001)
        h.SetFillColor(kRed)
        h.SetLineWidth(1)
        h.GetXaxis().CenterTitle()
        h.GetXaxis().SetTitle("silicon layer number")
        h.GetYaxis().CenterTitle()
        h.GetYaxis().SetTitle("ADC counts or energy (keV)")
        self.hist = h

    def __del__(self):
        pass

    def Draw(self, parent):
        if parent.hits == None: return

        # check if we are in accumulate mode
        if not parent.accumulate:
            self.hist.Reset()

        for ii in xrange(self.nlayers):
            layer = ii + 1
            cells = parent.cells[layer]
            total = 0.0
            for ii in xrange(cells.size()):
                total += cells[ii].count
            self.hist.Fill(layer, total)

        self.canvas.cd()
        self.hist.Draw('hist')
        self.canvas.Update()

        if parent.shutterOpen:
            filename = "longprofile%5.5d.png" % parent.eventNumber
            self.canvas.SaveAs(filename)
