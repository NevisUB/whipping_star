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
 *		BEGIN sbnfit_plot_covariance.cxx
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
        {"mcstat",no_argument,0,'m'},
        {"cmin",required_argument,0,'k'},
        {"cmax",required_argument,0,'p'},
        {0,			    no_argument, 		0,  0},
    };

    int iarg = 0;
    opterr=1;
    int index;
    std::string covar_file = "Stats_Only";
    bool stats_only = true;
   std::string signal_file;

    bool bool_incl_mcstat = false;

    bool bool_flat_det_sys = false;
    double flat_det_sys_percent = 0.0;

    bool remove_correlations = false;

    double cmin = 0;
    double cmax = -9;

    //a tag to identify outputs and this specific run. defaults to EXAMPLE1
    std::string tag = "TEST";

    while(iarg != -1)
    {
        iarg = getopt_long(argc,argv, "f:x:s:t:c:k:p:mzh", longopts, &index);

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
            case 'm':
                bool_incl_mcstat = true;
                break;

            case 't':
                tag = optarg;
                break;
     	    case 'c':
    	        covar_file=optarg;
                stats_only = false;
	          break;
            case 's':
                signal_file = optarg;
                break;

            case '?':
            case 'h':
                std::cout<<"---------------------------------------------------"<<std::endl;
                std::cout<<"sbnfit_plot_covariance allows for the plotting of covariance matricies from input root files containing reconstructed variables and covariance matricies. "<<std::endl;
                std::cout<<"---------------------------------------------------"<<std::endl;
                std::cout<<"--- Required arguments: ---"<<std::endl;
                std::cout<<"\t-x\t--xml\t\tInput configuration .xml file for SBNconfig"<<std::endl;
                std::cout<<"\t-t\t--tag\t\tA unique tag to identify the outputs [Default to TEST]"<<std::endl;
                std::cout<<"\t-s\t--signal\t\tInput signal SBNspec.root file"<<std::endl;
                std::cout<<"--- Optional arguments: ---"<<std::endl;
                std::cout<<"\t-f\t--flat\t\tAdd a flat percent systematic to fractional covariance matrix (all channels) (default false, pass in percent, i.e 5.0 for 5\% experimental)"<<std::endl;
                std::cout<<"\t-z\t--zero\t\tZero out all off diagonal elements of the systematics covariance matrix (default false, experimental!)"<<std::endl;
                std::cout<<"\t--cmax\t max for fractional covariance plot" << std::endl;
                std::cout<<"\t--cmin\t min for fractional covariance plot" << std::endl;
                std::cout<<"\t-m\t--mcstat\tInclude MC stat error?(Default No)"<<std::endl;
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

    std::cout<<"Begining Covariance Plotting for tag: "<<tag<<std::endl;
    std::cout<<"Loading SBNspec file : "<<signal_file<<" with xml "<<xml<<std::endl;
    SBNspec sig(signal_file,xml);
    if(!bool_incl_mcstat){
        std::cout<<"Not including MCStat error!"<<std::endl;
        sig.RemoveMCError();
    }else{
        std::cout<<"We are including MCStat error."<<std::endl;
    }
    sig.CalcFullVector();

    std::cout<<"Loading fractional covariance matrix from "<<covar_file<<std::endl;

    TFile * fsys;
    TMatrixD * cov;
    
    if(!stats_only){
        fsys = new TFile(covar_file.c_str(),"read");
        cov = (TMatrixD*)fsys->Get("frac_covariance");
   
        TMatrixD frac_flat_matrix(sig.num_bins_total, sig.num_bins_total);
        if(bool_flat_det_sys){
            std::cout << "RUNNING with flat systematics: " << flat_det_sys_percent << "%!" << std::endl;
            frac_flat_matrix.ResizeTo(sig.num_bins_total,sig.num_bins_total);
            frac_flat_matrix.Zero();
            for(int i=0 ; i< sig.num_bins_total; i++){
                   frac_flat_matrix(i,i)=flat_det_sys_percent*flat_det_sys_percent/10000.;
            }
            (*cov) = (*cov)+(frac_flat_matrix);
        }


        if(remove_correlations){
            std::cout<<"WARNING! We are running in   `Remove All Off Diagional Covariances/Correlations Mode` make sure this is what you want. "<<std::endl;
            for(int i=0; i<sig.num_bins_total;i++){ 
                for(int j=0; j<sig.num_bins_total;j++){ 
                    if(i==j)continue;
                    (*cov)(i,j) =0.0;
                }
            }
        }

       int num_nans = 0;
       for(int i=0; i<sig.num_bins_total;i++){ 
                for(int j=0; j<sig.num_bins_total;j++){ 
                    double val = (*cov)(i,j);
                    if( isinf(val) || isnan(val) || val!=val){
                            (*cov)(i,j) = 0.0;
                            num_nans++;
                    }
                }
            }
        std::cout<<"We have removed "<<num_nans<<" / "<<pow(sig.num_bins_total,2)<<" nans in the fractional covariance matrix"<<std::endl;
    }

    SBNchi SigChi(sig, *cov);
    SigChi.SetFracPlotBounds(cmin,cmax);
    SigChi.PrintMatricies(tag);
    sig.WriteOut(tag);

    
    std::cout << "Total wall time: " << difftime(time(0), start_time)/60.0 << " Minutes.\n";
    return 0;

}
