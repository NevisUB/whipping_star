#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <vector>
#include <unistd.h>
#include <getopt.h>
#include <float.h>
#include <cstring>
#include <gsl/gsl_math.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TH3.h"
#include "TMultiGraph.h"
#include "TString.h"
#include "TNtuple.h"
#include "TChain.h"
#include "TMath.h"
#include "TSystem.h"
#include "TMatrixT.h"
#include "TObjArray.h"
#include "TList.h"
#include "TRandom.h"
#include "TError.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TROOT.h"
#include "TMarker.h"

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
log_level_t GLOBAL_LEVEL = LOG_DEBUG;



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
        {"montecarlo", 		required_argument,	0, 'm'},
        {"data", 		required_argument,	0, 'd'},
        {"interpolation",       required_argument,      0, 'i'},
        {"covarmatrix",         required_argument,      0, 'c'},
        {"flat",        required_argument, 0 ,'f'},
        {"randomseed",        required_argument, 0 ,'r'},
        {"help", 		no_argument,	0, 'h'},
        {0,			    no_argument, 		0,  0},
    };

    int iarg = 0;
    opterr=1;
    int index;

    //a tag to identify outputs and this specific run. defaults to EXAMPLE1
    std::string tag;
    bool bool_stat_only = false;
    int interpolation_number = -99;  //number of points for chi2 value interpolation
    double random_number_seed = -1;

    bool input_data = false;
    std::string data_filename;  //root file containing data/mc spectra
    std::string mc_filename;
    std::string covmatrix_file;  //root file containing covariance matrix

    bool bool_flat_det_sys = false;
    double flat_det_sys_percent = 0.0;

    while(iarg != -1)
    {
        iarg = getopt_long(argc,argv, "c:d:x:m:r:p:i:f:sh", longopts, &index);

        switch(iarg)
        {
            case 'x':
                xml = optarg;
                break;
            case 'd':
                input_data = true;
                data_filename = optarg;
                break;
            case 'm':
                mc_filename = optarg;
                break;
            case 'i':
                interpolation_number = (int)strtod(optarg, NULL);
                break;
            case 'c':
                covmatrix_file = optarg;
                break;
            case 'f':
                bool_flat_det_sys = true;
                flat_det_sys_percent = (double)strtod(optarg,NULL);
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
                std::cout<<"sbnfit_fraction_fit is a work in progress."<<std::endl;
                std::cout<<"---------------------------------------------------"<<std::endl;
                std::cout<<"--- Required arguments: ---"<<std::endl;
                std::cout<<"\t-x\t--xml\t\tInput configuration .xml file for SBNconfig"<<std::endl;
                std::cout<<"\t-d\t--data\t\tInput observed data for a global scan"<<std::endl;
                std::cout<<"\t-m\t--montecarlo\t\tInput monte carlo for global scan"<<std::endl;
                std::cout<<"\t-i\t--interpolation\t\tInput number of points for interpolation"<< std::endl;
                std::cout<<"\t-c\t--covariance matrix\t\tInput syst covariance matrix for global scan"<< std::endl;
                std::cout<<"\t-f\t--flat\t\t Input flat systematic fractional covariance matrix"<<std::endl;
                std::cout<<"--- Optional arguments: ---"<<std::endl;
                std::cout<<"\t-s\t--stat\t\tStat only runs"<<std::endl;
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

    std::cout<<"Begining Fraction Fit"<<std::endl;

    NGrid mygrid;

    //now only available for 2 subchannels only
    mygrid.AddDimension("NCPi0Coh", 1.001, 1.002, 0.001);   //0.1 FULL
    //mygrid.AddDimension("NCPi0NotCoh", 0., 2, 0.05);   //0.1 FULL
    mygrid.AddDimension("NCDeltaRadOverlaySM", 0., 5, 0.05);   //0.1 FULL


    std::cout << "Fraction Fit|| "<< "\tStart initializing MC and data spectrum" << std::endl;

    //initialize the MC spectrum
    SBNspec mc_spec(mc_filename, xml);

    //initlaize the data spectrum
    SBNspec data_spec(data_filename, xml);
    data_spec.CollapseVector();  //collapse full vector

    //compare plots before the fit
    tag = "before_fit";
    TMatrixT<double> tt(data_spec.num_bins_total_compressed,data_spec.num_bins_total_compressed);
    tt.Zero();
    mc_spec.CompareSBNspecs(tt,&data_spec, tag);
    tag.clear();


    std::cout << "Fraction Fit||" <<  "\tInitialize fractional systematric covariance matrix" << std::endl;
    //initialize covariance matrix
    TMatrixT<double> frac_syst_matrix(mc_spec.num_bins_total, mc_spec.num_bins_total);


    if(bool_stat_only){
        frac_syst_matrix.Zero();
        std::cout<<"RUNNING Stat Only!"<<std::endl;
    }else if(bool_flat_det_sys){
        std::cout << "RUNNING with flat systematics: " << flat_det_sys_percent << "%!" << std::endl;
        frac_syst_matrix.Zero();
        //set up flat fractional syst covariance
        for(int i=0 ; i< mc_spec.num_bins_total; i++)
            for(int j=0;j<mc_spec.num_bins_total; j++)
                frac_syst_matrix(i,j)=flat_det_sys_percent*flat_det_sys_percent/10000.;
    }
    else{
        std::cout<< "RUNNING with Systematics!" << std::endl;
        frac_syst_matrix.Zero();

        std::cout << "Open covariance matrix root file: " << covmatrix_file << std::endl;
        TFile* f_covar = new TFile(covmatrix_file.c_str(), "read");
        TMatrixT<double>* m_temp;
        m_temp = (TMatrixD*)f_covar->Get("frac_covariance");
        frac_syst_matrix = *m_temp;
        f_covar->Close();	
    }



    std::cout<< "Fraction Fit||" << "\tGrab info of the grid"<< std::endl;
    //check grid size
    std::vector<std::vector<double>> grid = mygrid.GetGrid();
    if(grid.size() != mygrid.f_num_total_points){
        std::cout <<  "the number of points don't match: something wrong with the grid setup!!" << std::endl;
        return 1;
    }


    //collect the name of dimensions: subchannels you want to vary; and the range
    std::vector<std::string> dimension_name;
    const double range_x_low = mygrid.f_dimensions.at(0).f_min;
    const double range_x_up = mygrid.f_dimensions.at(0).f_max;
    const double range_y_low = mygrid.f_dimensions.at(1).f_min;
    const double range_y_up = mygrid.f_dimensions.at(1).f_max;
    //const double range_z_low = mygrid.f_dimensions.at(2).f_min;
    //const double range_z_up = mygrid.f_dimensions.at(2).f_max;
    int nbin_x = mygrid.f_dimensions.at(0).f_N;  //number of point in x axis
    int nbin_y = mygrid.f_dimensions.at(1).f_N;  
    //int nbin_z = mygrid.f_dimensions.at(2).f_N;  

    dimension_name.clear(); 
    for(int i=0; i< mygrid.f_num_dimensions ; i++){
        dimension_name.push_back(mygrid.f_dimensions.at(i).f_name);
    }



    //*********************loop over grid points, calc chi square********************************

    std::cout << "Fraction Fit||"<<  "\tStart GLOBAL SCAN" <<std::endl;
    std::vector<double> chi;  //vector to save chi square values.
    chi.reserve(grid.size());  //reserve the memory

    SBNchi chi_temp(xml);
    for(int i =0; i< grid.size() ;i++){

        std::cout << "on Point " << i << std::endl;
        //set a temperary SBNspec, assume there is already a MC CV root file with corresponding sequence defined in xml	
        SBNspec spec_temp(mc_filename, xml, false);

        //access grid point
        std::vector<double> point = grid[i];

        //check dimension point
        if(point.size() != dimension_name.size()){
            std::cout << "dimension doesn't match: something wrong with grid setup!" << std::endl;
            return 1;
        }

        //scale chosen subchannels
        for(int j=0; j< point.size(); j++ ){
            spec_temp.Scale(dimension_name[j], point[j]);
        }

        spec_temp.CollapseVector();

        //initialize a SBNchi
        //inverted collapsed covariance matrix	
        TMatrixT<double> collapsed_temp(spec_temp.num_bins_total_compressed, spec_temp.num_bins_total_compressed);
        TMatrixT<double> invert_collapsed_temp(spec_temp.num_bins_total_compressed, spec_temp.num_bins_total_compressed);
        collapsed_temp = chi_temp.CalcCovarianceMatrixCNP(frac_syst_matrix, spec_temp.full_vector, spec_temp.collapsed_vector, data_spec.collapsed_vector);
        invert_collapsed_temp = chi_temp.InvertMatrix(collapsed_temp);
        //calculate chi2
        double deter = collapsed_temp.Determinant();
        double chi_value = chi_temp.CalcChi(invert_collapsed_temp, spec_temp.collapsed_vector, data_spec.collapsed_vector)+log(deter);
        chi.push_back(chi_value);

        //print out chi square value
        std::cout << "CHI2: " << i << "/ "<<grid.size()<<" "<< chi.back();
        for(int j=0; j<point.size() ; j++) std::cout << " " << point[j];
        std::cout << std::endl;

    }

    

    std::cout << "Fraction Fit||" << "\tFinished" <<std::endl;
    std::cout << "Total wall time: " << difftime(time(0), start_time)/60.0 << " Minutes.\n";
    return 0;

}

