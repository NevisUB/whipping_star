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
#include "SBNfeld.h"

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

    std::string xml = "oscillate_example.xml";

    /*************************************************************
     *************************************************************
     *		Command Line Argument Reading
     ************************************************************
     ************************************************************/
    const struct option longopts[] =
    {
        {"xml", 		required_argument, 	0, 'x'},
        {"printall", 		no_argument, 		0, 'p'},
        {"stat", 		no_argument, 		0, 's'},
        {"number", 		required_argument,	0,'n'},
        {"tag", 		required_argument,	0, 't'},
        {"mode",        required_argument, 0 ,'m'},
        {"randomseed",        required_argument, 0 ,'r'},
        {"help", 		no_argument,	0, 'h'},
        {0,			    no_argument, 		0,  0},
    };

    int iarg = 0;
    opterr=1;
    int index;

    //a tag to identify outputs and this specific run. defaults to EXAMPLE1
    std::string tag = "TEST";
    std::string mode_option;
    bool bool_stat_only = false;
    int number = 2500;
    int grid_pt = 0;
    double random_number_seed = -1;

    while(iarg != -1)
    {
        iarg = getopt_long(argc,argv, "x:t:m:n:r:p:sh", longopts, &index);

        switch(iarg)
        {
            case 'x':
                xml = optarg;
                break;
            case 'n':
                number = (int)strtod(optarg,NULL);
                break;
            case 'p':
                grid_pt = (int)strtod(optarg,NULL);
                break;
            case 't':
                tag = optarg;
                break;
            case 'm':
                mode_option = optarg;
                break;
            case 'r':
                random_number_seed = (double)strtod(optarg,NULL);
                std::cout<<"Reading in random seed argument: "<<random_number_seed<<std::endl;
                break;
            case 's':
                bool_stat_only = true;
                break;
            case '?':
            case 'h':
                std::cout<<"---------------------------------------------------"<<std::endl;
                std::cout<<"sbnfit_feldman_cousins is a work in progress."<<std::endl;
                std::cout<<"---------------------------------------------------"<<std::endl;
                std::cout<<"--- Required arguments: ---"<<std::endl;
                std::cout<<"\t-x\t--xml\t\tInput configuration .xml file for SBNconfig"<<std::endl;
                std::cout<<"\t-t\t--tag\t\tA unique tag to identify the outputs [Default to TEST]"<<std::endl;
                std::cout<<"\t-m\t--mode\t\tWhat mode you want to run in. Arguments are:"<<std::endl;
                std::cout<<"\t\t\t--\t gen : Generate the preoscillated spectra for all mass-splittings"<<std::endl;  
                std::cout<<"\t\t\t--\t genbkg : Generate a background only spectra for all mass-splittings"<<std::endl;  
                std::cout<<"\t\t\t--\t feldman : Perform a full feldman cousins analysis"<<std::endl;  
                std::cout<<"\t-p\t--point\t\tWhat Grid Point to run over. -1 for a Full run (WARNING takes forever) [Defaults to 0]"<<std::endl;
                std::cout<<"\t\t\t--\t test : Just a testbed. Can be ignored"<<std::endl;
                std::cout<<"--- Optional arguments: ---"<<std::endl;
                std::cout<<"\t-s\t--stat\t\tStat only runs"<<std::endl;
                std::cout<<"\t-n\t--number\t\tNumber of pseudo-experiments to simulate (default 2500)"<<std::endl; 
                std::cout<<"\t-r\t--randomseed\t\tRandomNumber Seed (default from machine)"<<std::endl; 
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

    std::cout<<"Begining FeldmanCousins for tag: "<<tag<<std::endl;

    NGrid mygrid;
    mygrid.AddDimension("scale", 0.01, 20, 0.2);//0.1

    //Print the grid interesting bits
    mygrid.Print();
    SBNfeld myfeld(mygrid,tag,xml);

    if(mode_option == "feldman"){
    
        std::cout<<"Begininning a full Feldman-Cousins analysis for tag : "<<tag<<std::endl;
      
        NeutrinoModel nullModel(0,0,0);
        SBNgenerate  gen(xml,nullModel);
        gen.WriteCVSpec(tag);

        return 0;

        TMatrixT<double> *Msys = new TMatrixT<double>(8,8);
        Msys->Zero();

        //myfeld.SetFractionalCovarianceMatrix(tag+".SBNcovar.root","frac_covariance");
        myfeld.SetFractionalCovarianceMatrix(Msys);
        myfeld.SetStatOnly();
        myfeld.m_subchannel_to_scale = "nu_uBooNE_trk_signal";

        myfeld.SetCoreSpectrum(tag+"_CV.SBNspec.root");
        myfeld.GenerateScaledSpectra();
        myfeld.LoadBackgroundSpectrum(tag+"_CV.SBNspec.root");


        std::cout<<"Setting random seed "<<random_number_seed<<std::endl;
        myfeld.SetRandomSeed(random_number_seed);
        myfeld.SetNumUniverses(number);

        std::cout<<"Calculating the necessary SBNchi objects"<<std::endl;
        myfeld.CalcSBNchis();
        std::cout <<"DONE calculating the necessary SBNchi objects at : " << difftime(time(0), start_time)/60.0 << " Minutes.\n";

        std::cout<<"Beginning to peform FullFeldmanCousins analysis"<<std::endl;
        myfeld.FullFeldmanCousins();


    }else if(mode_option == "test"){

        myfeld.SetCoreSpectrum(tag+"_BKG_ONLY.SBNspec.root");
        if(bool_stat_only){
            myfeld.SetEmptyFractionalCovarianceMatrix();
            myfeld.SetStatOnly();
            std::cout<<"RUNNING Stat Only!"<<std::endl;
        }else{
            myfeld.SetFractionalCovarianceMatrix(tag+".SBNcovar.root","frac_covariance");
        }

        std::cout<<"Setting random seed "<<random_number_seed<<std::endl;
        myfeld.SetRandomSeed(random_number_seed);
        std::cout<<"Loading precomputed spectra"<<std::endl;
        myfeld.LoadPreOscillatedSpectra();
        myfeld.LoadBackgroundSpectrum();

        std::cout<<"Calculating the necessary SBNchi objects"<<std::endl;
        myfeld.CalcSBNchis();

        std::cout<<"Beginning to peform a globalScan analysis"<<std::endl;
        myfeld.GlobalScan("out_scaling.txt");

    }

    std::cout << "Total wall time: " << difftime(time(0), start_time)/60.0 << " Minutes.\n";
    return 0;

}

