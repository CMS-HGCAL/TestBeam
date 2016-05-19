#!/usr/bin/env python
#-----------------------------------------------------------------------------
# File:        TBADCCounts.py
# Description: TB 2016 - display of ADC distribution
# Created:     10-Apr-2016 Jeremy Thomas, Harrison B. Prosper
#-----------------------------------------------------------------------------
import sys, os, re
from string import atof, lower, replace, strip, split, joinfields, find
from HGCal.TBEventDisplay.TBUtil import *
from HGCal.TBStandaloneSimulator.TBGeometryUtil import *
from math import *
from ROOT import *
#------------------------------------------------------------------------------
class ADCCounts:

    def __init__(self, parent, page):
        self.cellmap = parent.cellmap
        self.geometry  = parent.geometry
        self.sensitive = parent.sensitive
        self.canvas  = page.canvas

        # try to figure out an arrangement of plots on the
        # canvas
        self.nlayers = len(self.sensitive)
        self.nplots  = divideCanvas(self.nlayers, self.canvas)

        # create a histogram for each sensor
        nbins = 128
        xmin  = 0
        xmax  = 128

        self.hist = []
        for layer in xrange(self.nplots):
            name = 'ADClayer%3d' % (layer+1)
            h = TH1F(name, "", nbins, xmin, xmax)
            h.SetFillStyle(3001)
            h.SetFillColor(kRed)
            h.SetLineWidth(1)
            h.GetXaxis().CenterTitle()
            h.GetXaxis().SetTitle("channel")
            h.GetYaxis().CenterTitle()
            h.GetYaxis().SetTitle("count")
            self.hist.append(h)

    def __del__(self):
        pass

    def Draw(self, parent):
        if parent.hits == None: return

        # check if we are in accumulate mode
        if not parent.accumulate:
            for h in self.hist:
                h.Reset()

        gStyle.SetOptStat("")
        ymax = 0.0
        for ii, h in enumerate(self.hist):
            layer = ii + 1
            cells = parent.cells[layer]
            for ii in xrange(cells.size()):
                cell   = cells[ii]
                skiroc = cell.skiroc
                channel= cell.channel
                ski = (skiroc+1) % 2
                jj = channel + 64 * ski + 1
                h.SetBinContent(jj, cell.count)
                if cell.count > ymax:
                    ymax = cell.count
            self.canvas.cd(layer)

        ymax *= 1.1
        for ii, h in enumerate(self.hist):
            layer = ii + 1
            self.canvas.cd(layer)
            h.SetMaximum(ymax)
            h.Draw()
        self.canvas.Update()

        if parent.shutterOpen:
            filename = "channels%5.5d.png" % parent.eventNumber
            self.canvas.SaveAs(filename)
