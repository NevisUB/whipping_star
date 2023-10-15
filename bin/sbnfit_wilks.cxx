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
#include "loghelper.h"
log_level_t GLOBAL_LEVEL = LOG_ERROR;


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
    std::string cov = "covariance_example.root";

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
        {"um4",        required_argument, 0 ,'u'},
        {"ue4",        required_argument, 0 ,'a'},
        {"begin",        required_argument, 0 ,'b'},
        {"end",        required_argument, 0 ,'e'},
        {"dlist",        required_argument, 0 ,'l'},
        {"cov", 		required_argument, 	0, 'c'},
        {"help", 		no_argument,	0, 'h'},
        {0,			    no_argument, 		0,  0},
    };

    int iarg = 0;
    opterr=1;
    int index;

    //a tag to identify outputs and this specific run. defaults to EXAMPLE1
    std::string tag = "TEST";
    std::string mode_option;
    std::string detector_list = "";
    bool bool_stat_only = false;
    int number = 10;
    int grid_pt = 0;
    double random_number_seed = -1;
    double begin = -1.0;
    double end = -1.0;
    double um4 = -1.0;
    double ue4 = -10.0;

    bool input_data = false;
    std::string data_filename;

    bool bool_flat_det_sys = false;
    double flat_det_sys_percent = 0.0;

    while(iarg != -1)
    {
        iarg = getopt_long(argc,argv, "d:x:t:c:m:n:r:p:f:s:b:e:u:a:l:h", longopts, &index);

        switch(iarg)
        {
            case 'x':
                xml = optarg;
                break;
            case 'd':
                input_data = true;
                data_filename = optarg;
                break;
            case 'c':
                cov = optarg;
                std::cout<<"If number < 10, using cov file as "<<cov<<std::endl;
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
            case 'b':
                begin = (double)strtod(optarg,NULL);
                std::cout<<"Setting m41 begin to "<<begin<<std::endl;
                break;
            case 'e':
                end = (double)strtod(optarg,NULL);
                std::cout<<"Setting m41 end to "<<end<<std::endl;
                break;
            case 'u':
                um4 = (double)strtod(optarg,NULL);
                std::cout<<"Setting Um4 is fixed to "<<um4<<std::endl;
                break;
            case 'a':
                ue4 = (double)strtod(optarg,NULL);
                std::cout<<"Setting Ue4 IS FIXED to "<<ue4<<std::endl;
                break;
            case 's':
                bool_stat_only = true;
                break;
            case 'l':
                detector_list = optarg;
                break;
            case '?':
            case 'h':
                std::cout<<"---------------------------------------------------"<<std::endl;
                std::cout<<"sbnfit_wilks is a work in progress."<<std::endl;
                std::cout<<"---------------------------------------------------"<<std::endl;
                std::cout<<"--- Required arguments: ---"<<std::endl;
                std::cout<<"\t-x\t--xml\t\tInput configuration .xml file for SBNconfig"<<std::endl;
                std::cout<<"\t-c\t--cov\t\tInput covariance file name"<<std::endl;
                std::cout<<"\t-d\t--data\t\tInput observed data for a global scan"<<std::endl;
                std::cout<<"\t-t\t--tag\t\tA unique tag to identify the outputs [Default to TEST]"<<std::endl;
                std::cout<<"\t-m\t--mode\t\tWhat mode you want to run in. Arguments are:"<<std::endl;
                std::cout<<"\t\t\t--\t gen : Generate the preoscillated spectra for all mass-splittings"<<std::endl;  
                std::cout<<"\t\t\t--\t genbkg : Generate a background only spectra for all mass-splittings"<<std::endl;  
                std::cout<<"\t-p\t--point\t\tWhat Grid Point to run over. -1 for a Full run (WARNING takes forever) [Defaults to 0]"<<std::endl;
                std::cout<<"\t\t\t--\t test : Just a testbed. Can be ignored"<<std::endl;
                std::cout<<"--- Optional arguments: ---"<<std::endl;
                std::cout<<"\t-s\t--stat\t\tStat only runs"<<std::endl;
                std::cout<<"\t-n\t--number\t\tNumber of pseudo-experiments to simulate (default 2500)"<<std::endl; 
                std::cout<<"\t-r\t--randomseed\t\tRandomNumber Seed (default from machine)"<<std::endl; 
                std::cout<<"\t-e\t--end\t\tEnding value for m4 scan"<<std::endl; 
                std::cout<<"\t-b\t--begin\t\tBegining value for m4 scan"<<std::endl; 
                std::cout<<"\t-u\t--um4\t\tValue for Um4 fixed (default -1.0)"<<std::endl; 
                std::cout<<"\t-a\t--ue4\t\tValue for Ue4 fixed (default scans)"<<std::endl; 
                std::cout<<"\t-l\t--dlist\t\tDetector list (MIS)"<<std::endl; 
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

    std::cout<<"Begining Wilks for tag: "<<tag<<std::endl;

    NGrid mygrid;
    //mygrid.AddDimension("m4", -1.0, 2, 0.05);   //jessie
    mygrid.AddDimension("m4", begin, end, 0.05);   //jessie
    //mygrid.AddFixedDimension("m4", begin);   //jessie
    if(ue4 == -10){
        //mygrid.AddDimension("ue4", -5/4., 0.7, 0.1);
        mygrid.AddDimension("ue4", -2.3, -0.3, 0.1);
        std::cout<<"Grid m4 start "<<begin << " and end " << end << ". ue4 start -5/4 and end 0.7. um4 = " << um4 <<std::endl;
    } 
    else{
        mygrid.AddFixedDimension("ue4", ue4);
        std::cout<<"Grid m4 start "<<begin << " and end " << end << ". ue4 = " << ue4 << ". um4 = " << um4 <<std::endl;
    }
    //mygrid.AddDimension("um4", -2, 0.05, 0.2); //jessie
    mygrid.AddFixedDimension("um4", um4); //jessie

    //Print the grid interesting bits
    mygrid.Print();
    cout << xml << endl;
    SBNfeld myfeld(mygrid,tag,xml);

    string spec_file;
    string cov_file;
    if(mode_option == "gen"){
        //myfeld.SetCoreSpectrum(tag+"_CV.SBNspec.root")
        spec_file = "/cluster/home/jmical01/uboone/analysis/zz_DataRelease_EventList_BNBNuMI/sbnfit_merge_completed_split_error_cosmics_3det.root";
        /*
        if(number == 11 || detector_list == "I" || detector_list == "MI"){
            spec_file = "/uboone/app/users/micallef/zz_DataRelease_EventList_BNBNuMI/sbnfit_merge_wICARUS.root";
        }
        else if(number == 12 || detector_list == "S" || detector_list == "MIS" || detector_list == "MS"){
            spec_file = "/uboone/app/users/micallef/zz_DataRelease_EventList_BNBNuMI/sbnfit_merge_completed_split_error_cosmics_3det.root";
        else{
            //spec_file = "/uboone/app/users/micallef/zz_DataRelease_EventList_BNBNuMI/sbnfit_merge.root";
            spec_file = "/uboone/app/users/micallef/zz_DataRelease_EventList_BNBNuMI/sbnfit_merge_completed_split_error_cosmics.root";
        }
        }*/
        myfeld.SetCoreSpectrum(spec_file);
        std::cout<< "Generate osc & background spec using "<< spec_file << std::endl;
        myfeld.GenerateOscillatedSpectra();
        myfeld.GenerateBackgroundSpectrum();


    }else if(mode_option == "genbkg"){
        spec_file = tag+"_CV.SBNspec.root";
        myfeld.SetCoreSpectrum(spec_file);
        myfeld.GenerateBackgroundSpectrum();
        std::cout<< "Set generate background spec using "<< spec_file << std::endl;

    }else if(mode_option == "test"){
        spec_file = tag+"_BKG_ONLY.SBNspec.root";
        myfeld.SetCoreSpectrum(spec_file);

        if(bool_stat_only){
            myfeld.SetEmptyFractionalCovarianceMatrix();
            myfeld.SetStatOnly();
            std::cout<<"RUNNING Stat Only!"<<std::endl;
        }
        else{
          if(number >= 0 && number < 10) {
            //myfeld.SetFractionalCovarianceMatrix(tag+".SBNcovar.root","frac_covariance_"+std::to_string(number));
            myfeld.SetFractionalCovarianceMatrix("/cluster/home/jmical01/uboone/analysis/zz_DataRelease_EventList_BNBNuMI/"+cov,"fracCOV_sum");
          }
          else if(detector_list == "M" || detector_list == "I" || detector_list == "S"){
              //cov_file = "/uboone/app/users/micallef/zz_DataRelease_EventList_BNBNuMI/sbnfit_file_systematics.root";
              //cov_file = "/uboone/app/users/micallef/zz_DataRelease_EventList_BNBNuMI/sbnfit_file_systematics_offdiagonal_TEST.root";
              cov_file = "/cluster/home/jmical01/uboone/analysis/zz_DataRelease_EventList_BNBNuMI/sbnfit_1detectors_file_systematics_offdiagonal_split_int_over_nodirt.root";
            myfeld.SetFractionalCovarianceMatrix(cov_file,"fracCOV_sum");
            std::cout<<"Using 1 TEST cov matrix" << std::endl;
          }
          else if(detector_list == "IS" || detector_list == "MI" || detector_list == "MS"){
            cov_file = "/cluster/home/jmical01/uboone/analysis/zz_DataRelease_EventList_BNBNuMI/sbnfit_2detectors_file_systematics_offdiagonal_split_int_over_nodirt.root";
            myfeld.SetFractionalCovarianceMatrix(cov_file,"fracCOV_sum");
            std::cout<<"Using 2 cov matrix" << std::endl;
          }
          else if(detector_list == "MIS"){
            cov_file = "/cluster/home/jmical01/uboone/analysis/zz_DataRelease_EventList_BNBNuMI/sbnfit_3detectors_file_systematics_offdiagonal_split_int_over_nodirt.root";
            myfeld.SetFractionalCovarianceMatrix(cov_file,"fracCOV_sum");
            std::cout<<"Using 2 cov matrix" << std::endl;
          }
          else{
              cov_file = tag+".SBNcovar.root";
            myfeld.SetFractionalCovarianceMatrix(cov_file,"frac_covariance");
          }
            std::cout<< "Testing using covariance matrix in "<< cov_file << std::endl;
        }


        if(bool_flat_det_sys){
            myfeld.AddFlatDetSystematic(flat_det_sys_percent);
        }

        std::cout<<"Setting random seed "<<random_number_seed<<std::endl;
        myfeld.SetRandomSeed(random_number_seed);
        std::cout<<"Loading precomputed spectra"<<std::endl;
        myfeld.LoadPreOscillatedSpectra();
        std::cout<<"Loading Background spectra"<<std::endl;
        myfeld.LoadBackgroundSpectrum();
	
        std::cout<<"Calculating the necessary SBNchi objects"<<std::endl;
        myfeld.CalcSBNchis();


        //everything up to here settng things up.



        std::cout<<"Beginning to peform a globalScan analysis"<<std::endl;
        //myfeld.GlobalScan(884);
        
        if(input_data){
            SBNspec * observed_spectrum = new SBNspec(data_filename, xml, false);
            myfeld.GlobalScan(observed_spectrum);
        }else{

            if(grid_pt==0){
		std::cout<<"globalScan analysis on grid pt == 0 "<<std::endl;
                myfeld.GlobalScan();//1503
            }else{
                myfeld.GlobalScan(grid_pt);
            }
            
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

