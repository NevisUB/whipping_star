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

#define no_argument 0
#define required_argument 1
#define optional_argument 2

using namespace sbn;

/*************************************************************
 *************************************************************
 *		BEGIN sbnfit_make_covariance.cxx
 ************************************************************
 ************************************************************/
int main(int argc, char* argv[])
{

	std::string xml = "example.xml";
	bool print_mode = false;

	/*************************************************************
	 *************************************************************
	 *		Command Line Argument Reading
	 ************************************************************
	 ************************************************************/
	const struct option longopts[] =
	{
		{"xml", 		required_argument, 	0, 'x'},
		{"printall", 		no_argument, 		0, 'p'},
		{"tag", 		required_argument,	0, 't'},
		{"help", 		no_argument,	0, 'h'},
		{0,			    no_argument, 		0,  0},
	};

	int iarg = 0;
	opterr=1;
	int index;

    //a tag to identify outputs and this specific run. defaults to EXAMPLE1
    std::string tag = "TEST";

	while(iarg != -1)
	{
		iarg = getopt_long(argc,argv, "x:t:ph", longopts, &index);

		switch(iarg)
		{
			case 'x':
				xml = optarg;
				break;
			case 'p':
				print_mode=true;
				break;
		    	case 't':
				tag = optarg;
				break;
		    	case '?':
			case 'h':
				std::cout<<"---------------------------------------------------"<<std::endl;
				std::cout<<"sbnfit_make_covariance allows for the building of covariance matricies from input root files containing reconstructed variables and the EventWeight class std::map<std::string, std::vector<double>>."<<std::endl;
				std::cout<<"---------------------------------------------------"<<std::endl;
				std::cout<<"--- Required arguments: ---"<<std::endl;
				std::cout<<"\t-x\t--xml\t\tInput configuration .xml file for SBNconfig"<<std::endl;
				std::cout<<"\t-t\t--tag\t\tA unique tag to identify the outputs [Default to TEST]"<<std::endl;
				std::cout<<"--- Optional arguments: ---"<<std::endl;
				std::cout<<"\t-p\t--printall\tRuns in BONUS print mode, making individual spectra plots for ALLVariations. (warning can take a while!) "<<std::endl;
                std::cout<<"\t-h\t--help\t\tThis help menu."<<std::endl;
				std::cout<<"---------------------------------------------------"<<std::endl;
	
				return 0;
		}
	}

	//std::string dict_location = "../libio/libEventWeight.so";
	//std::cout<<"Trying to load dictionary: "<<dict_location<<std::endl;
	//gSystem->Load(  (dict_location).c_str());

	/*************************************************************
	 *************************************************************
	 *			Main Program Flow
	 ************************************************************
	 ************************************************************/
	time_t start_time = time(0);
	
	std::cout<<"Begining Covariance Calculation for tag: "<<tag<<std::endl;

	//Create a SBNcovariance object initilizing with the inputted xml
	//This will load all the files and weights as laid out
	SBNcovariance example_covar(xml);

	//Form the covariance matrix from loaded weights and MC events
	example_covar.FormCovarianceMatrix(tag);

    //and make some plots of the resulting things
	//Will be outputted in the form: SBNfit_covariance_plots_TAG.root
	example_covar.PrintMatricies(tag);

	if(print_mode){
		//This takes a good bit longer, and prints every variation to file. 
		example_covar.PrintVariations(tag);
	}

	std::cout << "Total wall time: " << difftime(time(0), start_time)/60.0 << " Minutes.\n";
	return 0;

}
