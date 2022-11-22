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
 *		BEGIN sbnfit_scale_spec.cxx
 ************************************************************
 ************************************************************/
int main(int argc, char* argv[])
{

	std::string xml = "example.xml";
    std::string input;
    std::string output_tag;
    std::string scale_string;
    double scale_value;


	/*************************************************************
	 *************************************************************
	 *		Command Line Argument Reading
	 ************************************************************
	 ************************************************************/
	const struct option longopts[] =
	{
		{"xml", 		required_argument, 	0, 'x'},
		{"input", 		required_argument,	0, 'i'},
		{"scalestring", required_argument,	0, 's'},
		{"value", 		required_argument,	0, 'v'},
		{"tag", 		required_argument,	0, 't'},
		{"help", 		no_argument,	0, 'h'},
		{0,			    no_argument, 		0,  0},
	};

	int iarg = 0;
	opterr=1;
	int index;

	while(iarg != -1)
	{
		iarg = getopt_long(argc,argv, "x:i:t:s:v:h", longopts, &index);

		switch(iarg)
		{
			case 'x':
				xml = optarg;
				break;
			case 'i':
				input = optarg;
				break;
            case 't':
				output_tag = optarg;
				break;
            case 's':
				scale_string= optarg;
				break;
            case 'v':
                scale_value = strtod(optarg,NULL);
                break;
            case '?':
			case 'h':
		        std::cout<<"---------------------------------------------------"<<std::endl;
				std::cout<<"sbnfit_scale_spec allows for the simple scaling of one or more subchannels in a existing SBNspec.root file to produce another."<<std::endl;
				std::cout<<"---------------------------------------------------"<<std::endl;
				std::cout<<"--- Required arguments: ---"<<std::endl;
				std::cout<<"\t-x\t--xml\t\t\tInput configuration .xml file for SBNconfig"<<std::endl;
				std::cout<<"\t-t\t--tag\t\t\tA unique tag to identify the outputs, will be saved as TAG.SBNspec.root "<<std::endl;
                std::cout<<"\t-i\t--input\t\t\tInput SBNspec.root file"<<std::endl;
                std::cout<<"\t-s\t--scalestring\t\tAny subchannel that contains this string will be scaled by `value` "<<std::endl;
	            std::cout<<"\t-v\t--value\t\t\tWhat value do you want to scale by?"<<std::endl;
                std::cout<<"\t-h\t--help\t\t\tThis help menu."<<std::endl;
				std::cout<<"---------------------------------------------------"<<std::endl;
        return 0;
		}
	}
	/*************************************************************
	 *************************************************************
	 *			Main Program Flow
	 ************************************************************
	 ************************************************************/
	time_t start_time = time(0);

    SBNspec input_spectrum(input,xml);

    if(scale_string=="all") {
            input_spectrum.ScaleAll(scale_value);
}else{
            input_spectrum.Scale(scale_string,scale_value);
        }
    input_spectrum.WriteOut(output_tag);


	std::cout << "Total wall time: " << difftime(time(0), start_time)/60.0 << " Minutes.\n";
	return 0;

}
