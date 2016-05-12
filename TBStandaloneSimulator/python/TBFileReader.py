#!/usr/bin/env python
#-----------------------------------------------------------------------------
# File: TBFileReader.py
# Description: File reader for HGC Test Beam 2016
# Created: 09-Apr-2016 Harrison B. Prosper
#-----------------------------------------------------------------------------
import os, sys, re
from string import *
from ROOT import *
from DataFormats.FWLite import Events, Handle
#-----------------------------------------------------------------------------
gSystem.Load("libFWCoreFWLite.so")
FWLiteEnabler.enable()

simplify = re.compile(",edm::Strict.*[>] [>]")
extract  = re.compile("(?<=[<]).*(?=[>])")
dequote  = re.compile("\"")
#-----------------------------------------------------------------------------
class TBFileReader:
    def __init__(self, filename, **optional):
        
        if not os.path.exists(filename):
            sys.exit("** file not found: %s" % filename)

        # cache input variables 
        self.filename = filename

        # create options by scanning file
        print "reading file %s ..." % filename
        os.system("edmDumpEventContent %s >& .dumpevent" % filename)
        records = open(".dumpevent").readlines()
        for ii, record in enumerate(records):
            if record[:4] != "----":   continue
            records = records[ii+1:]
            break

        objects = {}
        for ii, record in enumerate(records):
            record = simplify.sub(">", strip(record))
            t = split(record)
            s = extract.findall(record)
            if len(s) > 0:
                name = s[0]
            else:
                name = t[0]
            objects[name] = (t[0], dequote.sub("",t[1]), dequote.sub("",t[2]))

        # create event buffer, get event iterator,
        # and create handles
        self.events = Events(filename)
        self.iter   = self.events.__iter__()
        self.event  = self.iter.next()
        self.numberOfEvents = self.event.size()
        print "\n  Number of events: %d" % self.numberOfEvents
        self.handles= {}

        print "  %-20s %s" % ("Key", "Type")
        print "-"*78
        keys = objects.keys()
        keys.sort()
        for key in keys:
            edmtype, label1, label2 = objects[key]
            h = Handle(edmtype)
            self.handles[key] = h
            print "  %-20s %s" % (key, edmtype)

        self.entry   = 0
        self.buffer  = {}
        self.objects = objects
        self.keys    = keys

    def __del__(self):
        pass

    # read event at specified index in TFile and 
    # cache requested objects
    def read(self, index):
        if index < 0: return False
        if index > self.numberOfEvents-1: return False
  
        self.event.to(index)

        for key in self.keys:
            handle = self.handles[key]
            if handle == None: continue
            edmtype, module, label = self.objects[key]
            self.event.getByLabel(module, label, handle)
            self.buffer[key] = handle.product()
        self.entry += 1
        return True

    # get next event and cache requested objects
    def next(self):
        try:
            self.event = self.iter.next()
        except StopIteration:
            return False

        for key in self.keys:
            handle = self.handles[key]
            if handle == None: continue
            edmtype, module, label = self.objects[key]
            self.event.getByLabel(module, label, handle)
            self.buffer[key] = handle.product()
        self.entry += 1
        return True

    # return object by key
    def __call__(self, key):
        try:
            return self.buffer[key]
        except:
            return None

    def entries(self):
        return self.numberOfEvents

    def names(self):
        return self.keys()

    def __len__(self):
        return self.entries()

#-----------------------------------------------------------------------------
def main():
    print "\n\t<== TBFileReader ==>\n"

    filename = "HGCal_digi_32GeV_electrons.root"
    cellmap  = HGCCellMap()
    reader   = TBFileReader( filename)
    index    = 0

    while reader.read(index):
        skiroc = reader("SKIROC2DataFrame")
        if skiroc == None: sys.exit("\n** cannot find skiroc collection\n")
        print "%d\tnumber of skiroc hits: %d" % (index, skiroc.size())
        for ii in xrange(skiroc.size()):
            digi   = SKIROC2DataFrame(skiroc2[ii])
            nsamples = digi.samples()
            detid = digi.detid()
            iu = detid.iu()
            iv = detid.iv()
            xy = cellmap.uv2xy(iu, iv)
            x  = xy.first
            y  = xy.second
            for jj in xrange(nsamples):
                adc = digi[jj].adcHigh()
                print "\t%5d: cell(%6.2f,%6.2f):        ADC = %d" % \
                    (jj, x, y, adc)
            print
        print
        index += 1
        break
#-----------------------------------------------------------------------------
if __name__ == "__main__": main()

