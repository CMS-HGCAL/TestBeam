
//-------------------------------------------------------------------------------------
// Authors: Shin-Shan Eiko Yu <syu@cern.ch>, Vieri Candelise <vieri.candelise@cern.ch>
//-------------------------------------------------------------------------------------

#include <fstream>
#include <iostream>
#include <TH2.h>
#include <TH1.h>
#include <TProfile.h>
#include <TF1.h>
#include <string>
#include <TFitResult.h>
#include <TMatrixD.h>
#include <TFile.h>
#include <TSystem.h>
#include <TProfile.h>

using namespace std;



void cutoff(double* locutoff, double* hicutoff, const TFitResult& low, const TFitResult& hi)
{
  double a=low.GetParams()[0];
  double b=low.GetParams()[1];

  TMatrixD cov_lo= low.GetCovarianceMatrix();
  double aerr2 = cov_lo[0][0];
  double aberr = cov_lo[0][1];
  double berr2 = cov_lo[1][1];

  double c=hi.GetParams()[0];
  double d=hi.GetParams()[1];

  TMatrixD cov_hi= hi.GetCovarianceMatrix();
  double cerr2 = cov_hi[0][0];
  double cderr = cov_hi[0][1];
  double derr2 = cov_hi[1][1];

  //  cout << aerr2 << " " << aberr << " " << berr2 << endl;
  //  cout << cerr2 << " " << cderr << " " << derr2 << endl;
  

  if(fabs(a-c)<1e-6)
    {
      locutoff[0]=-999;
      locutoff[1]=-999;
      hicutoff[0]=-999;
      hicutoff[1]=-999;
      return;
    }
  else{
    locutoff[0]= (d-b)/(a-c);
    
    double deriv_c =  (d-b)/(a-c)/(a-c);
    double deriv_a = -deriv_c;
    double deriv_d =  1./(a-c);
    double deriv_b = -deriv_d;

    double err2 = pow(deriv_a,2)*aerr2 + pow(deriv_b,2)*berr2+ 2*deriv_a*deriv_b*aberr +
      pow(deriv_c,2)*cerr2 + pow(deriv_d,2)*derr2 + 2*deriv_c*deriv_d*cderr;

    locutoff[1]= sqrt(err2);

    hicutoff[0]= (a*d-b*c)/(a-c);
    deriv_a = -c*(d-b)/(a-c)/(a-c);
    deriv_b = -c/(a-c);
    deriv_c = a*(d-b)/(a-c)/(a-c);
    deriv_d = a/(a-c);

    err2 = pow(deriv_a,2)*aerr2 + pow(deriv_b,2)*berr2+ 2*deriv_a*deriv_b*aberr +
      pow(deriv_c,2)*cerr2 + pow(deriv_d,2)*derr2 + 2*deriv_c*deriv_d*cderr;


    hicutoff[1]= sqrt(err2);


    return;
  }
}


void computeSaturation(string inputFile, bool fitProfile=false,string histoPrefix="HighGain_LowGain_2D_lct")
{
  bool hasType=false;
  if(histoPrefix.find("lct")!=std::string::npos)hasType=true;

  const int NTYPES= hasType? 6:64;
  const int NLAYERS=8;
  const int NCHIPS=2;

  TH1* h2D[NLAYERS][NCHIPS][NTYPES];

  TFile *f1 = TFile::Open(inputFile.data());

  TF1* flo = new TF1("flo","[0]*x+[1]");
  TF1* fhi = new TF1("fhi","[0]*x+[1]");
  TString prefix=gSystem->GetFromPipe(Form("file=%s; test=${file##*/}; test2=${test%%_HGC*}; echo \"${test2}\"",inputFile.data()));

  if(fitProfile)
    prefix = "profile_" + prefix;

  TString runNumber=gSystem->GetFromPipe(Form("file=%s; test2=${file%%_Reco.root*}; test=${test2##*_}; echo \"${test}\"",inputFile.data()));

  string runNumber_string = runNumber.Data();

  // profile can only be done for large amount of data

  ofstream fout;
  fout.open(Form("%s_%s.dat",prefix.Data(),histoPrefix.data()),ios::out | ios::app);

  ofstream fout_slope;
  fout_slope.open(Form("slope_%s_%s.dat",prefix.Data(),histoPrefix.data()),ios::out | ios::app);

  for(int il=0; il< NLAYERS; il++){
    for(int ic=0; ic< NCHIPS; ic++){
      for(int it=0; it< NTYPES; it++){

	cout << "runNumber = " << runNumber << "   layers = " << il+1 << ", skiroc chip " << ic+1 << " it = " << it << endl;

      h2D[il][ic][it] = (TH1*)(f1->FindObjectAny(Form("%s%d%02i%02i",histoPrefix.data(),il+1,ic+1,it)));
      h2D[il][ic][it] ->SetName(Form("h2D%d%02i%02i",il,ic,it));
      double par[4];
      double parerr[4];
      TFitResultPtr fitptr_lo;
      TFitResultPtr fitptr_hi;
      bool ispf = h2D[il][ic][it]->InheritsFrom(TProfile::Class());
      TProfile* pfx;
      if(ispf) // if the input is already TProfile
	{
	  pfx = (TProfile*)h2D[il][ic][it];

	  fitptr_lo= pfx->Fit("flo","S","",0,150); 
	  if(!(TF1*)pfx->GetFunction("flo"))continue;

	  fitptr_hi= pfx->Fit("fhi","S","",250,400); 
	  if(!(TF1*)pfx->GetFunction("fhi"))continue;
	}
      else // if the input is a TH2F
	{
	  if(fitProfile)
	    {
	      pfx = ((TH2F*)h2D[il][ic][it])->ProfileX();	  
	      fitptr_lo= pfx->Fit("flo","S","",0,150); 
	      if(!(TF1*)pfx->GetFunction("flo"))continue;
	      fitptr_hi= pfx->Fit("fhi","S","",250,400); 
	      if(!(TF1*)pfx->GetFunction("fhi"))continue;
	    }
	  else
	    {
	      fitptr_lo= ((TH2F*)h2D[il][ic][it])->Fit("flo","S","",0,150); 
	      if(!(TF1*)(((TH2F*)h2D[il][ic][it])->GetFunction("flo")))continue;
	      fitptr_hi= ((TH2F*)h2D[il][ic][it])->Fit("fhi","S","",250,400); 
	      if(!(TF1*)(((TH2F*)h2D[il][ic][it])->GetFunction("fhi")))continue;
	    }
	}

      TFitResult fitresult_lo = (*fitptr_lo);
      TFitResult fitresult_hi = (*fitptr_hi);
      if(fitresult_lo.IsEmpty())continue;
      if(fitresult_hi.IsEmpty())continue;

      double LG_cutoff[2], HG_cutoff[2];
      cutoff(LG_cutoff, HG_cutoff,fitresult_lo, fitresult_hi);
      cout << "LG cut off = " << LG_cutoff[0] << " +- " << LG_cutoff[1] << endl;
      cout << "HG cut off = " << HG_cutoff[0] << " +- " << HG_cutoff[1] << endl;

      fout << runNumber << " " << il+1 << " " << ic+1 << " " << it << " " << LG_cutoff[0] << " " << LG_cutoff[1] << " " << HG_cutoff[0] << " " << HG_cutoff[1] << endl;

      double slope = flo->GetParameter(0);
      double slopeerr = flo->GetParError(0);

      fout_slope << runNumber << " " << il+1 << " " << ic+1 << " " << it << " " << slope << " " << slopeerr << endl;
      
      } // loop over number of cell types
    } // end of loop over skiroc chips
  } // end of loop over layers
  fout.close();
  fout_slope.close();
  f1->Close();
}
