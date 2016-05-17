//-----------------------------------------------------------------------------
// Created: Summer-2008 Harrison B. Prosper
//-----------------------------------------------------------------------------
#include <algorithm>
#include "TRootHelpDialog.h"
#include "TGFileDialog.h"
#include "TGInputDialog.h"
#include <TGFont.h>
#include <TGListBox.h>
#include <TGResourcePool.h>

//#include "PhysicsTools/PyGui/interface/Dialog.h"
#include "HGCal/TBEventDisplay/interface/Dialog.h"

using namespace std;


Dialog::Dialog(const TGWindow* window,
               const TGWindow* main)
	: window_(window),
	  main_(main)
{}

string
Dialog::SelectFile(EFileDialogMode dlg_type,
                   string IniDir_,
                   string IniFilename_)
{
	TGFileInfo file_info;

	char* inidir = new char[IniDir_.size() + 1];
	copy(IniDir_.begin(), IniDir_.end(), inidir);
	inidir[IniDir_.size()] = 0;
	file_info.fIniDir = inidir;

	if ( IniFilename_ != "" ) {
		char* inifilename = new char[IniFilename_.size() + 1];
		copy(IniFilename_.begin(), IniFilename_.end(), inifilename);
		inifilename[IniFilename_.size()] = 0;
		file_info.fFilename = inifilename;
	}

	vector<string> ftypes;
	if ( dlg_type == kFDOpen ) {
		ftypes.push_back("Root files");
		ftypes.push_back("*.root");
		ftypes.push_back("All files");
		ftypes.push_back("*");
	} else {
		ftypes.push_back("All files");
		ftypes.push_back("*");
		ftypes.push_back("Root files");
		ftypes.push_back("*.root");
	}
	const char* filetypes[] = {ftypes[0].c_str(), ftypes[1].c_str(),
	                           ftypes[2].c_str(), ftypes[3].c_str(),
	                           0,            0
	                          };

	file_info.fFileTypes = filetypes;

	new TGFileDialog(window_, main_, dlg_type, &file_info);
	_filename = string(file_info.fFilename);
	_inidir   = string(file_info.fIniDir);
	return _filename;
}

Dialog::~Dialog() {}

string
Dialog::IniDir()
{
	return _inidir;
}


void Dialog::SetText(string title,
                     string text,
                     UInt_t w,
                     UInt_t h)
{
	TRootHelpDialog* d = new TRootHelpDialog(main_,
	        title.c_str(),
	        w, h);
	d->SetText(text.c_str());
	d->Popup();
}

string
Dialog::GetInput(string prompt, string defstr)
{
	char retstr[120];
	new TGInputDialog(window_, main_, prompt.c_str(), defstr.c_str(), retstr);
	return string(retstr);
}


FileDialog::FileDialog() {}

FileDialog::FileDialog(const TGWindow* window,
                       const TGWindow* main,
                       EFileDialogMode dlg_type,
                       string IniDir,
                       string IniFilename)
{
	TGFileInfo file_info;

	char* inidir = new char[IniDir.size() + 1];
	copy(IniDir.begin(), IniDir.end(), inidir);
	inidir[IniDir.size()] = 0;
	file_info.fIniDir = inidir;

	if ( IniFilename != "" ) {
		char* inifilename = new char[IniFilename.size() + 1];
		copy(IniFilename.begin(), IniFilename.end(), inifilename);
		inifilename[IniFilename.size()] = 0;
		file_info.fFilename = inifilename;
	}

	vector<string> ftypes;
	if ( dlg_type == kFDOpen ) {
		ftypes.push_back("Root files");
		ftypes.push_back("*.root");
		ftypes.push_back("All files");
		ftypes.push_back("*");
	} else {
		ftypes.push_back("All files");
		ftypes.push_back("*");
		ftypes.push_back("Root files");
		ftypes.push_back("*.root");
	}
	const char* filetypes[] = {ftypes[0].c_str(), ftypes[1].c_str(),
	                           ftypes[2].c_str(), ftypes[3].c_str(),
	                           0,            0
	                          };

	file_info.fFileTypes = filetypes;

	new TGFileDialog(window, main, dlg_type, &file_info);
	_filename = string(file_info.fFilename);
	_inidir   = string(file_info.fIniDir);
}

FileDialog::~FileDialog() {}

string
FileDialog::Filename()
{
	return _filename;
}

string
FileDialog::IniDir()
{
	return _inidir;
}

HelpDialog::HelpDialog() {}

HelpDialog::HelpDialog(const TGWindow* main,
                       string title,
                       string text,
                       UInt_t w,
                       UInt_t h)
{
	TRootHelpDialog* d = new TRootHelpDialog(main,
	        title.c_str(),
	        w, h);
	d->SetText(text.c_str());
	d->Popup();
}

HelpDialog::~HelpDialog() {}
