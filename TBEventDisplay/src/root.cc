//-----------------------------------------------------------------------------
/** PhysicsTools/PyGui/src/root.cc

 Description: A collection of simple Root utilities.

 Implementation:
     As simple as possible
*/
// Created: Summer-2008 Harrison B. Prosper
//-----------------------------------------------------------------------------
//$Revision: 1.1.1.1 $

#include <cassert>
#include <algorithm>

//#include "PhysicsTools/PyGui/interface/root.h"
#include "HGCal/TBEventDisplay/interface/root.h"
#include "TRootHelpDialog.h"
#include "TGFont.h"
#include "TGListBox.h"
#include "TGResourcePool.h"
#include "TGLUtil.h"
#include "TColor.h"

using namespace std;


const TGClient*
root::Client()
{
	return gClient;
}

const TGWindow*
root::GetRoot()
{
	return gClient->GetRoot();
}

Pixel_t
root::Color(std::string name)
{
	Pixel_t pixel;
	gClient->GetColorByName(name.c_str(), pixel);
	return pixel;
}

TGLBEntry*
root::GLBEntry(TGListBox* listbox,
               std::string str, int id,
               std::string font, int fontsize)
{
	// Create font object

	char fontstr[256];
	sprintf(fontstr, "-adobe-%s-*-*-%d-*-*-*-*-*-iso8859-1",
	        font.c_str(), fontsize);

	const TGFont* ufont = gClient->GetFont(fontstr);
	if (!ufont)
		ufont = gClient->GetResourcePool()->GetDefaultFont();

	// Create graphics context object

	GCValues_t val;
	val.fMask = kGCFont;
	val.fFont = ufont->GetFontHandle();
	TGGC* uGC = gClient->GetGC(&val, kTRUE);

	TGTextLBEntry* entry = new TGTextLBEntry(listbox->GetContainer(),
	        new TGString(str.c_str()),
	        id,
	        uGC->GetGC(),
	        ufont->GetFontStruct());
	return (TGLBEntry*)entry;
}

void
root::DrawAxes(TGLViewer* viewer)
{
	viewer->SetGuideState(TGLUtil::kAxesOrigin,
	                      kTRUE,
	                      kFALSE,
	                      0);
}

int
root::SetSpectrumPalette()
{
	TColor::SetPalette(1, 0);
	return TColor::GetNumberOfColors();
}
