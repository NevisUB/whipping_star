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
 *		BEGIN sbnfit_analyze_covariance.cxx
 ************************************************************
 ************************************************************/
int main(int argc, char* argv[])
{

    std::string xml = "example.xml";

    /*************************************************************
     *************************************************************
     *		Command Line Argument Reading
     ************************************************************
     ************************************************************/
    const struct option longopts[] =
    {
        {"xml", 		required_argument, 	0, 'x'},
        {"tag", 		required_argument,	0, 't'},
        {"signal", 		required_argument,	0, 's'},
        {"help", 		no_argument,	0, 'h'},
        {"covar",		required_argument,    0, 'c'},
        {"flat", required_argument,0,'f'},
        {"zero",no_argument,0,'z'},
        {"cmin",required_argument,0,'k'},
        {"cmax",required_argument,0,'p'},
        {0,			    no_argument, 		0,  0},
    };

    int iarg = 0;
    opterr=1;
    int index;
    std::vector<std::string> covar_files;
    std::vector<std::string> covar_names;
    bool stats_only = true;
    std::string signal_file;

    bool bool_flat_det_sys = false;
    double flat_det_sys_percent = 0.0;

    bool remove_correlations = false;

    double cmin = 0;
    double cmax = -9;

    //a tag to identify outputs and this specific run. defaults to EXAMPLE1
    std::string tag = "TEST";

    while(iarg != -1)
    {
        iarg = getopt_long(argc,argv, "f:n:x:s:t:c:k:p:zh", longopts, &index);

        switch(iarg)
        {
            case 'f':
                bool_flat_det_sys = true;
                flat_det_sys_percent = (double)strtod(optarg,NULL);
                break;
            case 'p':
                cmax = (double)strtod(optarg,NULL);
                break;
            case 'k':
                cmin = (double)strtod(optarg,NULL);
                break;
            case 'x':
                xml = optarg;
                break;
            case 'z':
                remove_correlations = true; 
                break;
            case 'c':
                covar_files.push_back(optarg);
                for (int i = optind; i < argc; i++) {
                    covar_files.push_back(argv[i]);
                }
                break;
            case 'n':
                covar_names.push_back(optarg);
                for (int i = optind; i < argc; i++) {
                    covar_names.push_back(argv[i]);
                }
                break;

            case 't':
                tag = optarg;
                break;
            case 's':
                signal_file = optarg;
                break;

            case '?':
            case 'h':
                std::cout<<"---------------------------------------------------"<<std::endl;
                std::cout<<"sbnfit_anayze_covariance allows for the analysis of covariance matricies from input root files containing reconstructed variables and covariance matricies. "<<std::endl;
                std::cout<<"---------------------------------------------------"<<std::endl;
                std::cout<<"--- Required arguments: ---"<<std::endl;
                std::cout<<"\t-x\t--xml\t\tInput configuration .xml file for SBNconfig"<<std::endl;
                std::cout<<"\t-t\t--tag\t\tA unique tag to identify the outputs [Default to TEST]"<<std::endl;
                std::cout<<"\t-s\t--signal\t\tInput signal SBNspec.root file"<<std::endl;
                std::cout<<"\t-c\t--covar\t\tList of covariances to plot (a.root b.root..z.root)"<<std::endl;
                std::cout<<"\t-n\t--names\t\tList of names to plot (a b z)"<<std::endl;
                std::cout<<"--- Optional arguments: ---"<<std::endl;
                std::cout<<"\t-h\t--help\t\tThis help menu."<<std::endl;
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

    int p = sbn::analyzeCovariance(xml,signal_file, tag, covar_files, covar_names);


    std::cout << "Total wall time: " << difftime(time(0), start_time)/60.0 << " Minutes.\n";
    return 0;

}
