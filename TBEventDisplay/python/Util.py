#-----------------------------------------------------------------------------
# File:        Util.py
# Package:     TB
# Description: Some gui utilities
# Created:     09-May-2014 Harrison B. Prosper & Sam Bein
#-----------------------------------------------------------------------------
import sys, os, re, platform
from ROOT import *
from string import lower, replace, strip, split, joinfields, find
from array import array
#-----------------------------------------------------------------------------
ICONDIR = "%s/src/HGCal/TBEventDisplay/icons" % os.environ["CMSSW_BASE"]
#-----------------------------------------------------------------------------
status = gSystem.Load('libFWCoreFWLite')
if status != 0:
     sys.exit('**unable to load libFWCoreFWLite')
else:
     print 'enabling auto loader'
     FWLiteEnabler.enable()

# GUI widget Layouts
RIGHT        = TGLayoutHints(kLHintsRight)
RIGHT_X      = TGLayoutHints(kLHintsRight | kLHintsExpandX)
RIGHT_X_Y    = TGLayoutHints(kLHintsRight | kLHintsExpandX | kLHintsExpandY)

LEFT         = TGLayoutHints(kLHintsLeft)
LEFT_X       = TGLayoutHints(kLHintsLeft | kLHintsExpandX)
LEFT_X_Y     = TGLayoutHints(kLHintsLeft | kLHintsExpandX | kLHintsExpandY)


TOP          = TGLayoutHints(kLHintsTop)
TOP_X        = TGLayoutHints(kLHintsTop  | kLHintsExpandX)
TOP_X_SUNKEN = TGLayoutHints(kLHintsTop  | kLHintsExpandX | kSunkenFrame)
TOP_LEFT_X   = TGLayoutHints(kLHintsTop  | kLHintsLeft | kLHintsExpandX)
TOP_LEFT_X_Y = TGLayoutHints(kLHintsTop  | kLHintsLeft | kLHintsExpandX |
							 kLHintsExpandY)
TOP_LEFT     = TGLayoutHints(kLHintsTop  | kLHintsLeft)
TOP_LEFT_PADDED     = TGLayoutHints(kLHintsTop  | kLHintsLeft,  5, 5, 2, 2)
TOP_RIGHT_PADDED    = TGLayoutHints(kLHintsTop  | kLHintsRight, 5, 5, 2, 2)
TOP_RIGHT_X_PADDED  = TGLayoutHints(kLHintsTop  | kLHintsRight | \
					    kLHintsExpandX, 5, 5, 2, 2)


TOP_RIGHT    = TGLayoutHints(kLHintsTop  | kLHintsRight)
TOP_X_Y      = TGLayoutHints(kLHintsTop  | kLHintsExpandX | kLHintsExpandY)
TOP_X_Y_RAISED = TGLayoutHints(kLHintsTop | kLHintsExpandX | kLHintsExpandY | \
				       kRaisedFrame)

K_PROG_MAX  = 20.0
C_STANDARD  = 0
C_OPENGL    = 1

BLACK       = root.Color("black")
WHITE       = root.Color("white")
RED         = root.Color("red")
ORANGE      = root.Color("orange")
YELLOW      = root.Color("yellow")
GREEN       = root.Color("green")
BLUE        = root.Color("blue")
DARKRED     = root.Color("darkred")
LIGHTYELLOW = root.Color("lightyellow")
LIGHTGREEN  = root.Color("lightgreen")
SKYBLUE     = (0.62, 0.57, 0.98)
#-----------------------------------------------------------------------------
GLCANVAS = '''
#include "TRootEmbeddedCanvas.h"
#include "TStyle.h"
TRootEmbeddedCanvas* TRootEmbeddedGLCanvas(const char* name, 
                                           TGWindow*    p, 
                                           unsigned int w, 
                                           unsigned int h)
{
  gStyle->SetCanvasPreferGL(true);
  TRootEmbeddedCanvas* embedded = new TRootEmbeddedCanvas(name, p, w, h);
  return embedded;
}
'''
gROOT.ProcessLine(GLCANVAS)

class Element:
	pass

class MenuBar(TGMenuBar):
	
	def __init__(self, obj, frame, layout=TOP_LEFT_PADDED):
		TGMenuBar.__init__(self, frame)
		frame.AddFrame(self, TOP_X)

		self.object   = obj
		self.frame    = frame
		self.layout   = layout
		self.number   = 0
		self.callbacks= {}
		self.elements = []

		
	def __del__(self):
		pass

	def Add(self, name, items):
		menu = TGPopupMenu(root.GetRoot())  	
		connection = Connection(menu, "Activated(Int_t)",
					self, "menu")
		self.AddPopup(name, menu, self.layout)
		self.elements.append((menu, connection))

		# add menu items
		for item in items:
			if type(item) == type(0):
				menu.AddSeparator()
			else:
				self.number += 1
				namen, method = item
				menu.AddEntry(namen, self.number)
				self.callbacks[self.number] = \
				    'self.object.%s()' % method

	# Responds to: Activated(Int_t)
	def menu(self, number):
		if self.callbacks.has_key(number):
			exec(self.callbacks[number])
		else:
			print "Unrecognized menu id = %d" % number


#-----------------------------------------------------------------------------
class ProgressBar(TGHProgressBar):
	
     def __init__(self, obj, toolBar, seconds=0.2):
          TGHProgressBar.__init__(self, toolBar,
                                  TGProgressBar.kFancy, 100)
		
          toolBar.AddFrame(self,
                           TGLayoutHints(100, 600, 12, 80))

          self.SetBarColor("green")
          self.SetRange(0, K_PROG_MAX)

          # Set up a timer for progress bars
          #self.timer = TTimer()
          #self.connection = Connection(self.timer, "Timeout()", object, method)

     def __del__(self):
          pass
#-----------------------------------------------------------------------------
buttonNumber=-1
class TextButton(TGTextButton):
     def __init__(self, obj, toolBar,
                  label,
                  method,
                  text='',
                  layout=kLHintsRight):
          global buttonNumber
          buttonNumber += 1
          number = buttonNumber
		
          TGTextButton.__init__(self, toolBar, 
                                label, number)
          toolBar.AddFrame(self,
                           TGLayoutHints(layout, 
                                         2, 2, 2, 2))
          self.SetToolTipText(text)
          self.connection = Connection(self, "Clicked()",
                                       obj, method)
     def __del__(self):
          pass
#-----------------------------------------------------------------------------
class PictureButton(TGPictureButton):
	
	def __init__(self, obj, toolBar,
		     picture,
		     method,
		     text='',
		     layout=kLHintsRight):
		global buttonNumber
		buttonNumber += 1
		number = buttonNumber
		
		self.pool = TGPicturePool(root.Client(), ICONDIR)
		if self.pool:
			self.picture = self.pool.GetPicture(picture)
			if self.picture:
				TGPictureButton.__init__(self, toolBar, 
							 self.picture, number)
				toolBar.AddFrame(self,
						 TGLayoutHints(layout, 
							       2, 2, 2, 2))
				self.SetToolTipText(text)
				self.connection = Connection(self, "Clicked()",
							     obj, method)
	def __del__(self):
		pass
#-----------------------------------------------------------------------------

class CheckButton(TGCheckButton):
	
	def __init__(self, obj, toolBar,
		     hotstring,
		     method,
		     text='',
		     layout=kLHintsRight):
		global buttonNumber
		buttonNumber += 1
		number = buttonNumber
		
		self.hotstring = hotstring
		if self.hotstring:
            		TGCheckButton.__init__(self, toolBar, 
					       self.hotstring, number)
            		toolBar.AddFrame(self,
                        	         TGLayoutHints(layout, 2, 2, 2, 2))
            	self.SetToolTipText(text)
            	self.connection = Connection(self, "Clicked()",
					     obj, method)
                
	def __del__(self):
		pass
#-----------------------------------------------------------------------------

class RadioButton(TGRadioButton):
	
	def __init__(self, obj, toolBar,
		     hotstring,
		     method,
		     text='',
		     layout=kLHintsRight):
		global buttonNumber
		buttonNumber += 1
		number = buttonNumber
		
		self.hotstring = hotstring
		if self.hotstring:
            		TGRadioButton.__init__(self, toolBar, 
					       self.hotstring, number)
            		toolBar.AddFrame(self,
                        	         TGLayoutHints(layout, 2, 2, 2, 2))
            	self.SetToolTipText(text)
            	self.connection = Connection(self, "Clicked()",
                                         obj, method)
	def __del__(self):
		pass

class VSlider(TGVSlider):
     def __init__(self, obj):
          pass
#-----------------------------------------------------------------------------
class NoteBook(TGTab):
	
	def __init__(self, parent, frame, method, width=800, height=600):
		TGTab.__init__(self, frame, 1, 1)
		
		frame.AddFrame(self, TOP_X_Y)
		self.connection = Connection(self, "Selected(Int_t)",
					     parent, method)
		self.parent= parent
		self.number=-1
		self.pages = {}
		self.names = {}
		self.width = width
		self.height= height
		self.pageNumber = 0
		self.page = None
		
	def __del__(self):
		pass

	def Add(self, name, sidebar=None):
		self.number += 1
		self.names[name] = self.number
		self.pages[self.number] = Element()
		element = self.pages[self.number]
		self.page = element
		
		element.name   = name
		element.redraw = True
		element.tab    = self.AddTab(name)
                # structure
                # +---+------------------+
                # |+--++----------------+|
                # ||  ||                ||
                # ||  ||                ||
                # ||  ||                ||
                # ||  ||                ||
                # |+--++----------------+|
                # +----------------------+
                element.hframe = TGHorizontalFrame(element.tab, 1, 1)
                element.tab.AddFrame(element.hframe, TOP_X_Y)

		# check for menu items
		if sidebar:
			element.sidebar = TGVerticalFrame(element.hframe, 
                                                          1, 1)
			element.hframe.AddFrame(element.sidebar, TOP_LEFT)
                        for code in sidebar: 
                             exec(code)

                element.display = TGHorizontalFrame(element.hframe, 1, 1)
                element.hframe.AddFrame(element.display, TOP_X_Y)

		from string import upper
                
		standardCanvas = find(upper(name), '3D') < 0
		if standardCanvas:
                     # set up a regular ROOT  canvas
                     if find(upper(name), 'LEGO') > -1:
                          element.ecanvas \
                              = TRootEmbeddedGLCanvas("c%s" % name,
                                                      element.display,
                                                      self.width,
                                                      self.height)
                     else:
                          element.ecanvas = TRootEmbeddedCanvas("c%s" % name,
                                                                element.display,
                                                                self.width,
                                                                self.height)
                     element.canvas  = element.ecanvas.GetCanvas()
                     element.display.AddFrame(element.ecanvas, TOP_X_Y)
		else:
                     # set up a canvas that can handle 3D displays
			
                     element.viewer  = TGLEmbeddedViewer(element.display)
                     element.display.AddFrame(element.viewer.GetFrame(),
                                              TOP_X_Y)
                     # set sky blue background color
                     bkg = element.viewer.ColorSet().Background()
                     bkg.SetColor(SKYBLUE[0], SKYBLUE[1], SKYBLUE[2])
                     # draw axes
                     root.DrawAxes(element.viewer)
                     
                     # we want our own simplified gui
                     #TEveManager.Create(kFALSE, 'l')
                     TEveManager.Create(kFALSE)

                     element.canvas = TEveViewer("Viewer")
                     element.canvas.SetGLViewer(element.viewer,
                                                element.viewer.GetFrame())
                     element.canvas.AddScene(gEve.GetEventScene())
                     gEve.GetViewers().AddElement(element.canvas)
			
                     element.shapes = []
                     element.fixedelements = TEveElementList("fixed")
                     element.elements = TEveElementList("transients")
                     gEve.AddElement(element.fixedelements)
                     gEve.AddElement(element.elements)

	def SetPage(self, idd):

		# Before changing current tab's color,
		# re-color previous tab to the
		# TGMainframe's default color

		tab = self.GetTabTab(self.pageNumber)
		tab.ChangeBackground(self.GetDefaultFrameBackground())

		# Now change tab and re-color tab

		if self.pages.has_key(idd):
			self.pageNumber = idd
			
		elif self.names.has_key(idd):
			self.pageNumber = self.names[idd]

		self.page = self.pages[self.pageNumber]
		self.GetTabTab(self.pageNumber).ChangeBackground(YELLOW)
		self.SetTab(self.pageNumber)
