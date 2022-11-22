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
#include "SBNcls.h"
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
 *		BEGIN unit2.cxx
 *		Some simple frequentist studies.
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
	
	std::cout<<"Begining unit test 2. "<<std::endl;

    int num_MC_events = 500000;


    SBNspec signal("unit1a.root","unit1a.xml");
    signal.Scale("leesignal",2.0);
    SBNspec bkg("unit1a.root","unit1a.xml");
    bkg.Scale("leesignal",0.0);

    TFile * fsys = new TFile("unit1a_matrix.root","read");
	TMatrixD * cov = (TMatrixD*)fsys->Get("TMatrixT<double>;1");
	TMatrixD * stat = (TMatrixD*)fsys->Get("TMatrixT<double>;1");
    
    stat->Zero();
    signal.CalcFullVector();
    for(int i=0; i<signal.num_bins_total; i++){
        (*stat)(i,i) = signal.full_vector[i];
    }


    SBNcls cls_factory_pois(&bkg, &signal,* cov);
    SBNcls cls_factory_cov(&bkg, &signal,* cov);
   	//cls_factory_pois.SetSampleCovariance();
    cls_factory_cov.SetSampleCovariance();


    cls_factory_pois.maxchival = 55;
    cls_factory_cov.maxchival  = 55;


    cls_factory_pois.CalcCLS(num_MC_events, "unit_test2_pois");
    cls_factory_cov.CalcCLS(num_MC_events, "unit_test2_cov");
 


	std::cout << "Total wall time: " << difftime(time(0), start_time)/60.0 << " Minutes.\n";
	return 0;

}
