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
        {"m4min",   optional_argument, 0, 'a'}, 
        {"m4max",   optional_argument, 0, 'b'}, 
        {"m4step",  optional_argument, 0, 'c'}, 
        {"m4index",   required_argument, 0, 'i'}, 
        {"u4index",  required_argument, 0, 'u'}, 
        {"chiname", required_argument, 0, 'o'},
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
    double m4min = -1.0;
    double m4max = 1.05;
    double m4step= 0.05;
    int m4index = -1;
    int u4index = -1; 
    // int number = 2500;
    int grid_pt = 0;
    double random_number_seed = -1;

    bool input_data = false;
    std::string data_filename;
    std::string chi_name;
    bool bool_flat_det_sys = false;
    double flat_det_sys_percent = 0.0;

    while(iarg != -1)
    {
        iarg = getopt_long(argc,argv, "d:x:t:m:n:r:p:f:a:b:c:i:u:o:sh", longopts, &index);

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
            case 'a':
                m4min =  (double)strtod(optarg,NULL);
                break;
            case 'b':
                m4max =  (double)strtod(optarg,NULL);
                break;
            case 'c':
                m4step = (double)strtod(optarg,NULL);
                break;
            case 'i':
                m4index = (int)strtod(optarg,NULL);
                break;
            case 'u':
                u4index = (int)strtod(optarg, NULL);
                break;
            case 'o':
                chi_name = optarg;
                std::cout<<chi_name<<std::endl;
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
 
    //grid for numu disappearance   
    //mygrid.AddDimension("m4", -1, 0, 0.5);
    //mygrid.AddDimension("m4", m4min, m4max, m4step);
    //mygrid.AddFixedDimension("ue4", 0);
    mygrid.AddDimension("m4", m4min, m4max, m4step);//0.05 
    //mygrid.AddDimension("um4",-2.0, 0.025, 0.025); //0.05       
    //mygrid.AddDimension("m4", -1.0, 0.05, 0.05);
    mygrid.AddDimension("ue4", -3.0, 0.05, 0.05);
    //mygrid.AddDimension("um4", -2.0, 0.025, 0.025); 
    mygrid.AddFixedDimension("um4", 0);
    std::vector<std::vector<double>> grid = mygrid.GetGrid();
    double ndim0 = mygrid.f_dimensions[0].GetNPoints();
    double ndim1 = mygrid.f_dimensions[1].GetNPoints();
    double ndim2 = mygrid.f_dimensions[2].GetNPoints();
    std::vector<double> delta_ms;
    for(int i=0; i<=ndim0; i++){
        delta_ms.push_back(m4min + i*m4step);
    }
    //Print the grid interesting bits
    mygrid.Print();
    
    if(mode_option == "gen"){
        if(m4index>=0){
          NGrid smaller_grid;
          smaller_grid.AddFixedDimension("m4", delta_ms[m4index]);
          smaller_grid.AddFixedDimension("ue4", -1.0);
          smaller_grid.AddFixedDimension("um4", 0);
          smaller_grid.Print();
          std::cout<<"m4 passed for index "<<m4index<<" is "<<delta_ms[m4index]<<std::endl;
          SBNfeld myfeld(smaller_grid,tag,xml);
          myfeld.GenerateOscillatedSpectra();
          myfeld.SetCoreSpectrum(tag+"_CV.SBNspec.root");
          if(m4index==0){
            myfeld.GenerateBackgroundSpectrum();
	  }
        }else{
          SBNfeld myfeld(mygrid,tag,xml);
          myfeld.GenerateOscillatedSpectra();
          myfeld.SetCoreSpectrum(tag+"_CV.SBNspec.root");
          myfeld.GenerateBackgroundSpectrum(); 
        }

    }else if(mode_option == "genbkg"){
        SBNfeld myfeld(mygrid,tag,xml);
        myfeld.SetCoreSpectrum(tag+"_CV.SBNspec.root");
        myfeld.GenerateBackgroundSpectrum();

    }else if(mode_option == "feldman"){
     
        SBNfeld myfeld(mygrid,tag,xml);

        std::cout<<"Begininning a full Feldman-Cousins analysis for tag : "<<tag<<std::endl;

        myfeld.SetCoreSpectrum(tag+"_BKG_ONLY.SBNspec.root");
        myfeld.SetFractionalCovarianceMatrix(tag+".SBNcovar.root","frac_covariance");

        std::cout<<"Setting random seed "<<random_number_seed<<std::endl;
        myfeld.SetRandomSeed(random_number_seed);
        std::cout<<"Loading precomputed spectra"<<std::endl;
        myfeld.LoadPreOscillatedSpectra();
        std::cout <<"DONE loading precomputed spectra at : " << difftime(time(0), start_time)/60.0 << " Minutes.\n";
        myfeld.LoadBackgroundSpectrum();

        myfeld.SetNumUniverses(number);

        std::cout<<"Calculating the necessary SBNchi objects"<<std::endl;
        myfeld.CalcSBNchis();
        std::cout <<"DONE calculating the necessary SBNchi objects at : " << difftime(time(0), start_time)/60.0 << " Minutes.\n";

        if(grid_pt ==-1){
            std::cout<<"Beginning to peform FullFeldmanCousins analysis"<<std::endl;
            myfeld.FullFeldmanCousins();
        }else if(grid_pt>=0){
             std::cout<<"Beginning to peform Single Grid PointFeldmanCousins analysis on pt: "<<grid_pt<<std::endl;
             myfeld.PointFeldmanCousins((size_t)grid_pt);
        }

    }else if(mode_option == "test"){
//        std::string chi_file = "/nevis/riverside/data/isafa/guanqun/whipping_star/build/bin/MCP2022A_NuMu/chisqs/"+chi_name+"_"+std::to_string(m4index)+"_"+std::to_string(u4index)+".txt";
        std::string chi_file = "/nevis/riverside/data/isafa/guanqun/whipping_star/build/bin/chisqs_0907/"+chi_name+"_"+std::to_string(m4index)+"_"+std::to_string(u4index)+".txt";
 
        if(m4index>=0){
          NGrid smaller_grid;
          smaller_grid.AddFixedDimension("m4", delta_ms[m4index]);
          //smaller_grid.AddDimension("m4", m4min, m4max, m4step);//0.05 
          //std::vector<double> ue4s;
          //for(int i=0; i<ndim1; i++){
          //  ue4s.push_back(-3.0 + i*0.05);}          
          //smaller_grid.AddFixedDimension("ue4", ue4s[u4index]); 
          //smaller_grid.AddFixedDimension("ue4", 0);
          smaller_grid.AddDimension("ue4", -3, 0.05, 0.05);
          //smaller_grid.AddDimension("um4", -2.0, 0.025, 0.025); 
          smaller_grid.AddFixedDimension("um4", 0); 
          smaller_grid.Print();
          std::cout<<"m4 passed for index "<<m4index<<" is "<<delta_ms[m4index]<<std::endl;

          SBNfeld myfeld(smaller_grid,tag,xml);
          myfeld.SetCoreSpectrum(tag+"_BKG_ONLY.SBNspec.root");

          if(bool_stat_only){
            myfeld.SetEmptyFractionalCovarianceMatrix();
            myfeld.SetStatOnly();
            std::cout<<"RUNNING Stat Only!"<<std::endl;
          }else{
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
          std::cout<<"Loading precomputed spectra"<<std::endl;
          myfeld.LoadPreOscillatedSpectra();
          myfeld.LoadBackgroundSpectrum();

          std::cout<<"Calculating the necessary SBNchi objects"<<std::endl;
          myfeld.CalcSBNchis();
          std::cout<<"Beginning to peform a globalScan analysis"<<std::endl;
          
        //if(input_data){
        //    SBNspec * observed_spectrum = new SBNspec(data_filename, xml, false);
            //myfeld.GlobalScan(observed_spectrum);   //with fixed covariance matrix during each iteration
      //    myfeld.GlobalScanVary(observed_spectrum); 	//with covariance matrix varying at each grid point
            //myfeld.GlobalScanNeyman(observed_spectrum);
      //  }else{

         if(grid_pt==0){
             myfeld.GlobalScan(chi_file);//1503
         }else{
             myfeld.GlobalScan(grid_pt, chi_file);
             }
            
      //  }

       }else{
            SBNfeld myfeld(mygrid,tag,xml);
            myfeld.SetCoreSpectrum(tag+"_BKG_ONLY.SBNspec.root");
            if(bool_stat_only){
                myfeld.SetEmptyFractionalCovarianceMatrix();
                myfeld.SetStatOnly();
                std::cout<<"RUNNING Stat Only!"<<std::endl;
            }else{
  
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
            std::cout<<"Loading precomputed spectra"<<std::endl; 
            myfeld.LoadPreOscillatedSpectra(); 
            myfeld.LoadBackgroundSpectrum();
            std::cout<<"Finished loading everything <<<<<<<<<<<<<<<<<<<<<<<"<<std::endl;
            std::cout<<"Calculating the necessary SBNchi objects"<<std::endl;
            myfeld.CalcSBNchis();

            std::cout<<"Beginning to peform a globalScan analysis"<<std::endl;
         
            if(grid_pt==0){
                myfeld.GlobalScan(chi_file);//1503
            }else{
                myfeld.GlobalScan(grid_pt, chi_file);
                }
            
            }    
            
      
    }else if(mode_option == "plot"){
        SBNfeld myfeld(mygrid,tag,xml);

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

