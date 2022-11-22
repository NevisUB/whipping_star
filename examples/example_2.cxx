#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <vector>
#include <unistd.h>
#include <getopt.h>
#include <cstring>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TString.h"
#include "TNtuple.h"
#include "TChain.h"
#include "TMath.h"
#include "TSystem.h"
#include "TMatrixT.h"
#include "TRandom.h"
#include "TStyle.h"
#include "TError.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TGraph.h"

#include "SBNconfig.h"
#include "SBNchi.h"
#include "SBNspec.h"
#include "SBNosc.h"
#include "SBNfit.h"
#include "SBNfit3pN.h"
#include "SBNgenerate.h"
#include "prob.h"
#include "loghelper.h"
log_level_t GLOBAL_LEVEL = LOG_DEBUG;

#define no_argument 0
#define required_argument 1
#define optional_argument 2



using namespace sbn;

int main(int argc, char* argv[])
{
	std::string xml = "example.xml";
	int iarg = 0;
	opterr=1;
	int index;
	bool gen = false;
	bool numudis = false;
	bool combined = false;
	int mass_start = -1;

	const struct option longopts[] =
	{
		{"xml", 		required_argument, 	0, 'x'},
		{"gen",	no_argument, 0, 'g'},
		{"dis",	no_argument,0,'d'},
		{"comb", no_argument,0,'c'},
		{"part", required_argument,0,'p'},
		{0,			no_argument, 		0,  0},
	};

	while(iarg != -1)
	{
		iarg = getopt_long(argc,argv, "x:bscp:g", longopts, &index);

		switch(iarg)
		{
			case 'x':
				xml = optarg;
				break;
			case 'g':
				gen = true;
				break;
			case 'd':
				numudis = true;
				break;
			case 'c':
				combined = true;
				break;
			case 'p':
				mass_start = atoi(optarg);
				break;
			case '?':
			case 'h':
				std::cout<<"Allowed arguments:"<<std::endl;
				std::cout<<"\t-x\t--xml\t\tInput .xml file for SBNconfig"<<std::endl;
				return 0;
		}
	}

	std::string tag = "EXAMPLE2";

	//Load up the central value spectra we computed in example 1, to act as a signal spectra
    	SBNspec sig_spectra("EXAMPLE1.SBNspec.root",xml);
    	
	//Repeat but this time scale the LEE signal subchannel to 0, to act as our background spectra
	SBNspec bkg_spectra("EXAMPLE1.SBNspec.root",xml);
	bkg_spectra.Scale("nu_uBooNE_nue_leesignal", 0.0); // Scale() will scale all that match so bkg_spectra.Scale("leesignal",0.0); would work too 

	//Load up our covariance matricies we calculated in example1 (we could also load up single variation ones)
	TFile * fsys = new TFile("EXAMPLE1.SBNcovar.root","read");
	TMatrixD * cov = (TMatrixD*)fsys->Get("frac_covariance");
	
	//Create two SBNchi objects, if you pass a TMatrixD it will set that as systematics, otherwise it will be a statonly.
	SBNchi *chi_statonly = new SBNchi(bkg_spectra);
	SBNchi *chi = new SBNchi(bkg_spectra,*cov);

	std::vector<double> ans_chi;
	std::vector<double> ans_chi_statonly;
	std::vector<double> scaling;

	for(double k=-2; k<=1; k+=0.025){	
		//Create a tempoary copy of our signal, scaled up or down by some amount
		SBNspec tmp = sig_spectra;
		tmp.Scale("nu_uBooNE_nue_leesignal",pow(10,k));	

		//calculate the chi^2 between this temp signal and the bkg_spectra with and without spectra
		double chi2 = chi->CalcChi(&tmp);
		double chi2_statonly = chi_statonly->CalcChi(&tmp);

		std::cout<<"Signal Scaling: "<<pow(10,k)<<" chi^2: "<<chi2<<" chi^2(statonly): "<<chi2_statonly<<std::endl;

		//And just fill a few things for plotting below	
		ans_chi.push_back(chi2);			
		ans_chi_statonly.push_back(chi2_statonly);			
		scaling.push_back(pow(10,k));			
	}



	//and just plot it quickly
	TCanvas *c = new TCanvas();
	c->cd();

	TGraph *g_chi = new TGraph(scaling.size(), &scaling[0],&ans_chi[0]);
	TGraph *g_chi_statonly = new TGraph(scaling.size(), &scaling[0],&ans_chi_statonly[0]);

	g_chi->SetLineColor(kRed-7);	
	g_chi_statonly->SetLineColor(kBlue-7);	

	g_chi->Draw("acl");
	g_chi->GetXaxis()->SetTitle("Signal Scaling");
	g_chi->GetYaxis()->SetTitle(" #chi^{2}");

	g_chi_statonly->Draw("cl same");


	TLegend *l = new TLegend(0.49,0.11,0.89,0.49);
	l->AddEntry(g_chi,"Sys+Stats","l");
	l->AddEntry(g_chi_statonly,"Stat Only","l");	
	l->SetLineColor(kWhite);
	l->SetFillStyle(0);
	l->Draw();

	c->SaveAs((tag+"_simple_scaling.pdf").c_str(),"pdf");

	return 0;
}
