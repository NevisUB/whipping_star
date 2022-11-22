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
#include "TError.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TGraph.h"

#include "params.h"
#include "SBNconfig.h"
#include "SBNchi.h"
#include "SBNspec.h"
#include "SBNosc.h"
#include "SBNfit.h"
#include "SBNfit3pN.h"
#include "SBNcovariance.h"
#include "loghelper.h"
log_level_t GLOBAL_LEVEL = LOG_DEBUG;


#define no_argument 0
#define required_argument 1
#define optional_argument 2

using namespace sbn;

/*************************************************************
 *************************************************************
 *		BEGIN unit1.cxx
 ************************************************************
 ************************************************************/
int main(int argc, char* argv[])
{


	/*************************************************************
	 *************************************************************
	 *		Command Line Argument Reading
	 ************************************************************
	 ************************************************************/
	const struct option longopts[] =
	{
		{0,			no_argument, 		0,  0},
	};

	int iarg = 0;
	opterr=1;
	int index;
    int mode = 1;

    //a tag to identify outputs and this specific run. defaults to EXAMPLE1
    std::string tag = "EXAMPLE1";

	while(iarg != -1)
	{
		iarg = getopt_long(argc,argv, "", longopts, &index);

		switch(iarg)
		{
            case '?':
			case 'h':
				std::cout<<"No Allowed arguments for this unit test."<<std::endl;
				return 0;
		}
	}

    /*************************************************************
	 *************************************************************
	 *			Main Program Flow
	 ************************************************************
	 ************************************************************/
	time_t start_time = time(0);
	
	std::cout<<"Begining unit test 1. "<<std::endl;


    SBNspec signal_A("unit1a.root","unit1a.xml");
    SBNspec bkg_A("unit1a.root","unit1a.xml");
    bkg_A.Scale("leesignal",0.0);
 
    SBNspec signal_B("unit1b.root","unit1b.xml");
    SBNspec bkg_B("unit1b.root","unit1b.xml");
    bkg_B.Scale("leesignal",0.0);
 
    SBNspec signal_C("unit1c.root","unit1c.xml");
    SBNspec bkg_C("unit1c.root","unit1c.xml");
    bkg_C.Scale("leesignal",0.0);

    signal_A.WriteOut("unit1_signal_A");
    signal_B.WriteOut("unit1_signal_B");
    signal_C.WriteOut("unit1_signal_C");


   	SBNchi *sbnchi_A_statonly = new SBNchi(bkg_A);
   	SBNchi *sbnchi_B_statonly = new SBNchi(bkg_B);
   	SBNchi *sbnchi_C_statonly = new SBNchi(bkg_C);

    double chi_A_statonly = sbnchi_A_statonly->CalcChi(&signal_A);
    double chi_B_statonly = sbnchi_B_statonly->CalcChi(&signal_B);
    double chi_C_statonly = sbnchi_C_statonly->CalcChi(&signal_C);


    std::cout<<"Unit test 1: Stat only, Analytical answer should be: 2.336902 for all THREE below."<<std::endl;
    std::cout<<"A: "<<chi_A_statonly<<" B: "<<chi_B_statonly<<" C: "<<chi_C_statonly<<std::endl;
    
    
    TFile * fsysa = new TFile("unit1a_matrix.root","read");
	TMatrixD * cova = (TMatrixD*)fsysa->Get("TMatrixT<double>;1");
	 
    TFile * fsysb = new TFile("unit1b_matrix.root","read");
	TMatrixD * covb = (TMatrixD*)fsysb->Get("TMatrixT<double>;1");
	

    SBNchi *sbnchi_A_statplussys = new SBNchi(bkg_A,*cova);
	SBNchi *sbnchi_B_statplussys = new SBNchi(bkg_B,*covb);
    SBNchi *sbnchi_C_statplussys = new SBNchi(bkg_C,*covb);

    double chi_A_statplussys = sbnchi_A_statplussys->CalcChi(&signal_A);
    double chi_B_statplussys = sbnchi_B_statplussys->CalcChi(&signal_B);
    double chi_C_statplussys = sbnchi_C_statplussys->CalcChi(&signal_C);

    std::cout<<"Unit test 1: Stat+ Systematics, Analytical answer should be: 0.3222 for all THREE below."<<std::endl;
    std::cout<<"A: "<<chi_A_statplussys<<" B: "<<chi_B_statplussys<<" C: "<<chi_C_statplussys<<std::endl;
    




	std::cout << "Total wall time: " << difftime(time(0), start_time)/60.0 << " Minutes.\n";
	return 0;

}
