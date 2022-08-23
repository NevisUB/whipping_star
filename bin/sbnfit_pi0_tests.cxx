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
        {"data", 		required_argument,	0, 'd'},
        {"mode",        required_argument, 0 ,'m'},
        {"flat",        required_argument, 0 ,'f'},
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
    int number = -1;
    int grid_pt = 0;
    double random_number_seed = -1;

    bool input_data = false;
    std::string data_filename;

    bool bool_flat_det_sys = false;
    double flat_det_sys_percent = 0.0;

    while(iarg != -1)
    {
        iarg = getopt_long(argc,argv, "d:x:t:m:n:r:p:f:sh", longopts, &index);

        switch(iarg)
        {
            case 'x':
                xml = optarg;
                break;
            case 'd':
                input_data = true;
                data_filename = optarg;
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
            case 'f':
                bool_flat_det_sys = true;
                flat_det_sys_percent = (double)strtod(optarg,NULL);
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
                std::cout<<"\t-d\t--data\t\tInput observed data for a global scan"<<std::endl;
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

    mygrid.AddDimension("NCPi0Coh", 0.0, 5, 0.05);
    mygrid.AddDimension("NCPi0NotCoh",0.5, 1.0, 0.025); 

    //Print the grid interesting bits
    mygrid.Print();
    SBNfeld myfeld(mygrid,tag,xml);

    if(mode_option == "genbkg"){
        myfeld.SetCoreSpectrum(tag+"_CV.SBNspec.root");
        myfeld.GenerateBackgroundScaledSpectrum();
    
    }else if(mode_option == "fit"){

        std::cout<<"Setting Core Spectrum"<<std::endl;
        myfeld.SetCoreSpectrum(tag+"_CV.SBNspec.root");
        if(bool_stat_only){
            myfeld.SetEmptyFractionalCovarianceMatrix();
            myfeld.SetStatOnly();
            std::cout<<"RUNNING Stat Only!"<<std::endl;
        }else{

          std::cout<<"Loading Covariance Matrix"<<std::endl;
          if(number >= 0) {
            myfeld.SetFractionalCovarianceMatrix(tag+".SBNcovar.root","frac_covariance_"+std::to_string(number));
          }else{
            myfeld.SetFractionalCovarianceMatrix(tag+".SBNcovar.root","frac_covariance");
          }
        }

        if(bool_flat_det_sys){
            myfeld.AddFlatDetSystematic(flat_det_sys_percent);
        }

        std::cout<<"Setting random seed "<<random_number_seed<<std::endl;
        myfeld.SetRandomSeed(random_number_seed);

        
        std::cout<<"Loading Scaled Spectra"<<std::endl;
        myfeld.GenerateScaledSpectra();
        myfeld.SetBackgroundSpectrum(tag+"_CV.SBNspec.root","",1.0);

        std::cout<<"Calculating the necessary SBNchi objects"<<std::endl;
        myfeld.CalcSBNchis();


        //everything up to here settng things up.
        std::cout<<"Beginning to peform a globalScan analysis"<<std::endl;
        
        if(input_data){
            std::cout<<"Loading Data "<<data_filename<<std::endl;
            SBNspec * observed_spectrum = new SBNspec(data_filename, xml, false);
            myfeld.GlobalScan(observed_spectrum);
        }else{

            if(grid_pt==0){
                myfeld.GlobalScan();//1503
            }else{
                myfeld.GlobalScan(grid_pt);
            }
            
        }

    }else if(mode_option=="genie"){

        SBNspec * observed_spectrum = new SBNspec(data_filename, xml, false);
        SBNspec * cv = new SBNspec(tag+"_CV.SBNspec.root",xml,false);
        SBNchi * chi = new SBNchi(*cv, true);

        TFile *f = new TFile("ALL.SBNcovar.root");
        std::vector<std::string> All_Genie= {"AGKYpT1pi_UBGenie_frac_covariance","AGKYxF1pi_UBGenie_frac_covariance","AhtBY_UBGenie_frac_covariance","All_UBGenie_frac_covariance","AxFFCCQEshape_UBGenie_frac_covariance","BhtBY_UBGenie_frac_covariance","CV1uBY_UBGenie_frac_covariance","CV2uBY_UBGenie_frac_covariance","CoulombCCQE_UBGenie_frac_covariance","DecayAngMEC_UBGenie_frac_covariance","EtaNCEL_UBGenie_frac_covariance","FrAbs_N_UBGenie_frac_covariance","FrAbs_pi_UBGenie_frac_covariance","FrCEx_N_UBGenie_frac_covariance","FrCEx_pi_UBGenie_frac_covariance","FrInel_pi_UBGenie_frac_covariance","FrPiProd_N_UBGenie_frac_covariance","FrPiProd_pi_UBGenie_frac_covariance","FracDelta_CCMEC_UBGenie_frac_covariance","FracPN_CCMEC_UBGenie_frac_covariance","MFP_N_UBGenie_frac_covariance","MFP_pi_UBGenie_frac_covariance","MaCCQE_UBGenie_frac_covariance","MaCCRES_UBGenie_frac_covariance","MaNCEL_UBGenie_frac_covariance","MaNCRES_UBGenie_frac_covariance","MvCCRES_UBGenie_frac_covariance","MvNCRES_UBGenie_frac_covariance","NonRESBGvbarnCC1pi_UBGenie_frac_covariance","NonRESBGvbarnCC2pi_UBGenie_frac_covariance","NonRESBGvbarnNC1pi_UBGenie_frac_covariance","NonRESBGvbarnNC2pi_UBGenie_frac_covariance","NonRESBGvbarpCC1pi_UBGenie_frac_covariance","NonRESBGvbarpCC2pi_UBGenie_frac_covariance","NonRESBGvbarpNC1pi_UBGenie_frac_covariance","NonRESBGvbarpNC2pi_UBGenie_frac_covariance","NonRESBGvnCC1pi_UBGenie_frac_covariance","NonRESBGvnCC2pi_UBGenie_frac_covariance","NonRESBGvnNC1pi_UBGenie_frac_covariance","NonRESBGvnNC2pi_UBGenie_frac_covariance","NonRESBGvpCC1pi_UBGenie_frac_covariance","NonRESBGvpCC2pi_UBGenie_frac_covariance","NonRESBGvpNC1pi_UBGenie_frac_covariance","NonRESBGvpNC2pi_UBGenie_frac_covariance","NormCCCOH_UBGenie_frac_covariance","NormCCMEC_UBGenie_frac_covariance","NormNCCOH_UBGenie_frac_covariance","NormNCMEC_UBGenie_frac_covariance","RPA_CCQE_Reduced_UBGenie_frac_covariance","RPA_CCQE_UBGenie_frac_covariance","RootinoFix_UBGenie_frac_covariance","ThetaDelta2NRad_UBGenie_frac_covariance","Theta_Delta2Npi_UBGenie_frac_covariance","TunedCentralValue_UBGenie_frac_covariance","VecFFCCQEshape_UBGenie_frac_covariance","XSecShape_CCMEC_UBGenie_frac_covariance"};

        
        for(int i=0; i< All_Genie.size(); i++){
            std::cout<<"Starting on "<<All_Genie[i]<<std::endl;
            TMatrixT<double>* M = (TMatrixT<double>*)f->Get(("individualDir/"+All_Genie[i]).c_str());
             std::cout<<"Loaded M "<<M->IsValid()<<std::endl;
            cv->CollapseVector();
            std::cout<<"AAGH1"<<std::endl;
            TMatrixT<double> background_full_covariance_matrix = chi->CalcCovarianceMatrix(M, cv->full_vector);
            std::cout<<"AAGH2"<<std::endl;
            TMatrixT<double> background_collapsed_covariance_matrix(cv->num_bins_total_compressed, cv->num_bins_total_compressed);
            std::cout<<"AAGH3"<<std::endl;
            chi->CollapseModes(background_full_covariance_matrix, background_collapsed_covariance_matrix);    

            std::cout<<"And WriteOut"<<std::endl;
            cv->CompareSBNspecs(background_collapsed_covariance_matrix,observed_spectrum,All_Genie[i]);   

        }



    }else if(mode_option == "plot"){
        if(number<0)number =1503; 
        myfeld.SetCoreSpectrum(tag+"_BKG_ONLY.SBNspec.root");
        myfeld.SetEmptyFractionalCovarianceMatrix();
        myfeld.SetStatOnly();
        myfeld.SetRandomSeed(random_number_seed);
        myfeld.LoadPreOscillatedSpectrum(number);

    }

    std::cout << "Actually done with everything! " << std::endl;
    std::cout << "Total wall time: " << difftime(time(0), start_time)/60.0 << " Minutes.\n";
    return 0;

}

