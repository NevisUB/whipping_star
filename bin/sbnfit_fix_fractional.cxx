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
 *		BEGIN sbnfit_fix_fractional.cxx
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
        {"help", 		no_argument,	0, 'h'},
        {"zero", 		no_argument,	0, 'z'},
    	{"covar",		required_argument,    0, 'c'},
        {0,			    no_argument, 		0,  0},
    };

    int iarg = 0;
    opterr=1;
    int index;
    std::string covar_file = "Stats_Only";
    bool stats_only = true;
   std::string signal_file;

    bool bool_flat_det_sys = false;
    double flat_det_sys_percent = 0.0;
    bool run_in_silly = false;

    bool remove_correlations = false;

    //a tag to identify outputs and this specific run. defaults to EXAMPLE1
    std::string tag = "TEST";

    while(iarg != -1)
    {
        iarg = getopt_long(argc,argv, "f:x:s:t:c:zh", longopts, &index);

        switch(iarg)
        {
            case 'x':
                xml = optarg;
                break;
            case 'z':
                run_in_silly = true;
                    break;

            case 't':
                tag = optarg;
                break;
     	    case 'c':
    	        covar_file=optarg;
	          break;
            case '?':
            case 'h':
                std::cout<<"---------------------------------------------------"<<std::endl;
                std::cout<<"sbnfit_fix_fractional allows for the fixing and modification of a fractional covariance matrix . "<<std::endl;
                std::cout<<"---------------------------------------------------"<<std::endl;
                std::cout<<"--- Required arguments: ---"<<std::endl;
                std::cout<<"\t-x\t--xml\t\tInput configuration .xml file for SBNconfig"<<std::endl;
                std::cout<<"\t-t\t--tag\t\tA unique tag to identify the outputs [Default to TEST]"<<std::endl;
                std::cout<<"\t-c\t--covarance\t use this root file for a covariance matrix" << std::endl;
                std::cout<<"--- Optional arguments: ---"<<std::endl;
                std::cout<<"\t-h\t--help\t\tThis help menu."<<std::endl;
                std::cout<<"\t-z\t--zero\t\t insead of doign anything clever, zero instead"<<std::endl;
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

    SBNspec core(xml);

    std::cout<<"Loading fractional covariance matrix from "<<covar_file<<std::endl;
    TFile * fsys = new TFile(covar_file.c_str(),"read");
    TMatrixD cov =  *(TMatrixD*)fsys->Get("frac_covariance");
    TMatrixD beforeCov = cov; 
    
    std::cout<<"Matrix is of size: "<<cov.GetNrows()<<" and core spec: "<<core.num_bins_total<<std::endl;

    for(int i=0 ; i< core.num_bins_total; i++){
                    double val = (cov)(i,i);

                    if( isinf(val) || isnan(val) || val!=val){
                        int wh = core.GetHistNumber(i);
                        std::cout<<val<<" @ "<<i<<" in hist "<<wh<<" and subchannel "<<core.hist[wh].GetName()<<std::endl;
                        std::vector<int> bins = core.map_tag_to_covariance_index[core.hist[wh].GetName()];
                           
                        //Calculate the max number in the bins. 
                        double max_err = -999;
                        for(int j=bins[0]; j<=bins[1];j++){
                               double dval = cov(j,j);
                               if(isinf(dval) || isnan(dval) || dval!=dval) continue;
                               if(dval>max_err)max_err = dval; 
                        }

                        if(max_err <0){
                             //Add a check for BNBext or subchannels we do not want to set.
                            if(core.vec_is_data[wh]==1){
                                std::cout<<"Although 0 in subchannel, this is assigned as data so no scaling." <<std::endl;
                                max_err = 0.0;
                            }else{
                                std::cout<<"There are no non-zero in subchannel, setting to 100% "<<std::endl;
                                max_err = 1.0;
                            }
                        }else{
                            std::cout<<"The max error in other bins is "<<max_err<<std::endl;
                        }


                        std::cout<<"Setting them to be uncorrelated"<<std::endl;
                        for(int j=0 ; j< core.num_bins_total; j++){
                            if(i==j)continue;
                            cov(i,j) = 0.0;
                            cov(j,i) = 0.0;
                        }
                        
                        if(!run_in_silly){
                            cov(i,i) =  max_err;
                        }else{
                            cov(i,i) = 0.0;
                        }
                    }
    }
    
    TFile *fout = new TFile((tag+"_fracfixed.SBNcovar.root").c_str(),"recreate");
    fout->cd();
    cov.Write("frac_covariance", TObject::kWriteDelete);
    beforeCov.Write("original_frac_covariance", TObject::kWriteDelete);
    fout->Close();

    std::cout << "Total wall time: " << difftime(time(0), start_time)/60.0 << " Minutes.\n";
    return 0;

}
