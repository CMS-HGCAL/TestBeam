#!/usr/bin/env python
#-----------------------------------------------------------------------------
# File:        TBLego.py
# Description: TB 2016 simple lego display of HGCal
# Created:     10-Apr-2016 Jeremy Thomas, Harrison B. Prosper
#-----------------------------------------------------------------------------
import sys, os, re
from string import atof, lower, replace, strip, split, joinfields, find
from HGCal.TBEventDisplay.TBUtil import *
from HGCal.TBStandaloneSimulator.TBGeometryUtil import *
from math import *
from ROOT import *
#------------------------------------------------------------------------------
class Lego:

    def __init__(self, parent, page):
        self.parent  = parent
        self.page    = page
        self.canvas  = page.canvas
        self.cellmap = parent.cellmap        
        self.geometry  = parent.geometry
        self.sensitive = parent.sensitive

        # try to figure out an arrangement of plots on the
        # canvas
        self.nlayers = len(self.sensitive)
        self.nplots  = divideCanvas(self.nlayers, self.canvas)

        # create honeycomb histograms
        self.hist = []
        for l in xrange(self.nplots):
            layer =  l + 1
            element = self.geometry[self.sensitive[layer]]
            if not element.has_key('cellsize'):
                sys.exit('** TBLego - cellsize not found')
            if not element.has_key('side'):
                sys.exit('** TBLego - side not founds')

            cellside= element['cellsize']
            side    = element['side']
            poly = TH2Poly()
            poly.SetName('lego %3d' % layer)
            poly.SetTitle('lego %3d' % layer)
            poly.GetXaxis().CenterTitle()
            poly.GetXaxis().SetTitle("#font[12]{x} axis")
            poly.GetYaxis().CenterTitle()
            poly.GetYaxis().SetTitle("#font[12]{y} axis")

            # populate histogram with cells
            cells = parent.cells[layer]
            for ii, cell in enumerate(cells):
                xv, yv = computeBinVertices(cellside, cell)
                poly.AddBin(len(xv), xv, yv)
            self.hist.append(poly)

    def __del__(self):
        pass

    def Draw(self, parent):
        if parent.hits == None: return

        for l, h in enumerate(self.hist):
            layer = l + 1
            cells = parent.cells[layer]
            for ii, cell in enumerate(cells):
                h.SetBinContent(ii+1, cell.count)
            h.SetMaximum(parent.maxCount)
            self.canvas.cd(layer)
            h.Draw("legogl")
        self.canvas.Update()

        if parent.shutterOpen:
            filename = "lego%5.5d.png" % parent.eventNumber
            self.canvas.SaveAs(filename)

