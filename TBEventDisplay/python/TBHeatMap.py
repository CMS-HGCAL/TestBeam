#!/usr/bin/env python
#-----------------------------------------------------------------------------
# File:        TBHeatMap.py
# Description: TB 2016 simple heat map display of HGCal
# Created:     10-Apr-2016 Jeremy Thomas, Harrison B. Prosper
#-----------------------------------------------------------------------------
import sys, os, re
from string import atof, lower, replace, strip, split, joinfields, find
from HGCal.TBEventDisplay.TBUtil import *
from HGCal.TBStandaloneSimulator.TBGeometryUtil import *
from math import *
from ROOT import *
#------------------------------------------------------------------------------
class HeatMap:

    def __init__(self, parent, page):
        self.parent  = parent
        self.page    = page
        self.canvas  = page.canvas
        self.cellmap = parent.cellmap
        self.hist    = parent.hist

        self.geometry  = parent.geometry
        self.sensitive = parent.sensitive

        # try to figure out an arrangement of plots on the
        # canvas
        self.nlayers = len(self.sensitive)
        self.nplots  = divideCanvas(self.nlayers, self.canvas)

        # construct a hexagon centered at the origin to
        # represent sensor
        layer   = 1
        element = self.geometry[self.sensitive[layer]]
        side    = element['side']
        self.wafer = TH2Poly()
        self.wafer.SetName('wafer')
        self.wafer.SetTitle('wafer')
        self.wafer.GetXaxis().CenterTitle()
        self.wafer.GetXaxis().SetTitle("#font[12]{x} axis")
        self.wafer.GetYaxis().CenterTitle()
        self.wafer.GetYaxis().SetTitle("#font[12]{y} axis")
        xv, yv  = computeHexVertices(side)
        self.wafer.AddBin(len(xv), xv, yv)

        # draw a wafer outline for each wafer
        for l in xrange(self.nplots):
            layer = l + 1
            self.canvas.cd(layer)
            self.wafer.Draw()
        self.canvas.Update()

    def __del__(self):
        pass

    def Draw(self, parent):
        if parent.hits == None: return

        gStyle.SetOptStat("")
        self.text = TText()
        self.text.SetTextSize(0.02)
        self.text.SetTextAlign(22)  # centered

        for l, h in enumerate(self.hist):
            layer = l + 1
            cells = parent.cells[layer]

            self.canvas.cd(layer)
            h.Draw("colz")
            self.wafer.Draw("same")
            if len(self.hist) > 2: continue

            for ii in xrange(cells.size()):
                cell = cells[ii]
                if cell.count < parent.ADCmin: continue
                self.text.DrawText(cell.x, cell.y, "%d" % cell.count)

        self.canvas.Update()

        if parent.shutterOpen:
            filename = "heatmap%5.5d.png" % parent.eventNumber
            self.canvas.SaveAs(filename)

