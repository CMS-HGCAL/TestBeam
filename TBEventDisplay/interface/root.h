#ifndef ROOT_H
#define ROOT_H
// -*- C++ -*-
//
// Package:    PhysicsTools/PyGui
/**
 Description: A class of Root utilities. These functions are placed in a class
              so that Reflex can handle overloading automatically. This is
	      just a collection of simple boilerplate code to lessen
	      clutter that I've written over the years.

 Implementation:
     As simple as possible
*/
//
// Original Author:  Harrison B. Prosper
//         Created:  Fri Apr 04 2008
// $Id: root.h,v 1.1.1.1 2011/05/04 13:04:28 prosper Exp $
//
//
//-----------------------------------------------------------------------------
#include <iostream>
#include <vector>
#include <string>
//-----------------------------------------------------------------------------
#include "TGClient.h"
#include "TGWindow.h"
#include "TRootHelpDialog.h"
#include "TGFont.h"
#include "TGListBox.h"
#include "TGResourcePool.h"
#include "TGLUtil.h"
#include "TGLViewer.h"
//-----------------------------------------------------------------------------

struct root {
	virtual ~root() {}
	///
	static
	const TGWindow* GetRoot();

	///
	static
	const TGClient* Client();

	///
	static
	Pixel_t Color(std::string name);

	///
	static
	TGLBEntry* GLBEntry(TGListBox* listbox, std::string str, int id,
	                    std::string font = "helvetica-medium-r", int fontsize = 14);

	///
	static
	void DrawAxes(TGLViewer* viewer);

	///
	static
	int SetSpectrumPalette();
};

#endif
