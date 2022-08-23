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
 *		BEGIN sbnfit_merge_fractional.cxx
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
        {"tag", 		required_argument,	0, 't'},
        {"force", 		no_argument,	0, 'f'},
        {"zero", 		no_argument,	0, 'z'},
        {"covars", 		required_argument,	0, 'c'},
        {"help", 		no_argument,	0, 'h'},
        {0,			    no_argument, 		0,  0},
    };

    int iarg = 0;
    opterr=1;
    int index;
    std::vector<std::string> covar_files;
    bool stats_only = true;
    std::string signal_file;
    bool force = false;

    bool bool_flat_det_sys = false;
    double flat_det_sys_percent = 0.0;

    bool remove_correlations = false;

    //a tag to identify outputs and this specific run. defaults to EXAMPLE1
    std::string tag = "TEST";
    int oi = 0;
    while(iarg != -1)
    {
        iarg = getopt_long(argc,argv, "fs:c:t:zh", longopts, &index);
        switch(iarg)
        {
            case 't':
                tag = optarg;
                break;
            case 'c':
                covar_files.push_back(optarg);
                for (int i = optind; i < argc; i++) {
                    covar_files.push_back(argv[i]);
                }
                break;
            case 'f':
                force = true;
                break;
            case 'z':
                remove_correlations=true;
                break;
            case '?':
            case 'h':
                std::cout<<"---------------------------------------------------"<<std::endl;
                std::cout<<"sbnfit_merge_fractional allows for the fixing and modification of a fractional covariance matrix . "<<std::endl;
                std::cout<<"---------------------------------------------------"<<std::endl;
                std::cout<<"--- Required arguments: ---"<<std::endl;
                std::cout<<"\t-t\t--tag\t\tA unique tag to identify the outputs [Default to TEST]"<<std::endl;
                std::cout<<"\t-c\t--covars\t\tWhat SBNcovars to merge"<<std::endl;
                std::cout<<"--- Optional arguments: ---"<<std::endl;
                std::cout<<"\t-f\t--force\t\t Forces even if nans present"<<std::endl;
                std::cout<<"\t-f\t--zero\t\t Zeros out off diag"<<std::endl;
                std::cout<<"\t-h\t--help\t\tThis help menu."<<std::endl;
                std::cout<<"---------------------------------------------------"<<std::endl;

                return 0;
        }
        oi++;
    }

    /*************************************************************
     *************************************************************
     *			Main Program Flow
     ************************************************************
     ************************************************************/
    time_t start_time = time(0);

    std::cout<<"The Covar File string is length: "<<covar_files.size()<<std::endl;
    for(auto &f: covar_files) std::cout<<" "<<f<<std::endl;


    std::cout<<"Writing now"<<std::endl;
    TFile *f  = new TFile((tag+".SBNcovar.root").c_str(),"RECREATE");
    f->cd();


    std::vector<TFile*> files;
    std::vector<TMatrixD> mats;
    std::vector<TMatrixD> fmats;

    for(auto s: covar_files){
        std::cout<<"Loading "<<s.c_str()<<std::endl;
        files.push_back(new TFile(s.c_str(),"read") );
        TMatrixD fm = *(TMatrixD*)files.back()->Get("frac_covariance");
        std::cout<<s<<" has Dimensions: "<<fm.GetNrows()<<" "<<fm.GetNcols()<<std::endl;
        fmats.push_back(fm);

        files.back()->Close();
        int mnans = 0;
        for(int i=0; i< fmats.back().GetNrows(); i++){
            for(int j=0; j< fmats.back().GetNrows(); j++){
                if(isnan(fmats.back()(i,j)) || fmats.back()(i,j)!=fmats.back()(i,j) || isinf(fmats.back()(i,j))){

                    if(!force){
                        std::cout<<"ERROR ERROR We have a NAN or INF , at "<<i<<" "<<j<<" "<<fmats.back()(i,j)<<std::endl;
                        std::cout<<"In Matrix "<<s<<std::endl;
                        std::cout<<"Shouldnt be the case, fix before merging or it will break"<<std::endl;
                        exit(EXIT_FAILURE);
                    }else{
                        fmats.back()(i,j)=0.0;
                    }
                }
            }
        }  
    }

    f->cd();
    TMatrixD m_fsum = fmats[0];
    std::cout<<fmats.size()<<std::endl;

    for(int i=1; i< fmats.size(); i++){
        m_fsum = m_fsum + fmats[i];
        std::cout<<"Now adding : "<<i<<std::endl;
    }

    if(remove_correlations){
        std::cout<<"Removing Correlation, i.e off diagonals"<<std::endl;
        for(int i=0; i< fmats.back().GetNrows(); i++){
            for(int j=0; j< fmats.back().GetNrows(); j++){
                if(i!=j) m_fsum(i,j)=0.0;
            }
        }
    }

    f->cd();
    m_fsum.Write("frac_covariance",TObject::kWriteDelete);
    //m_fsum.Write("frac_covariance");

    f->Close();
    std::cout << "Total wall time: " << difftime(time(0), start_time)/60.0 << " Minutes.\n";
    return 0;

}
