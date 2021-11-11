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
#include "TH2D.h"
#include "THStack.h"
#include "TLegend.h"
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
#include "TLatex.h"

#include "params.h"
#include "SBNconfig.h"
#include "SBNchi.h"
#include "SBNspec.h"
#include "SBNosc.h"
#include "SBNfit.h"
#include "SBNfit3pN.h"
#include "SBNcovariance.h"
#include "SBNsinglephoton.h"

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


    /*************************************************************
     *************************************************************
     *		Command Line Argument Reading
     ************************************************************
     ************************************************************/
    const struct option longopts[] =
    {
        {"xml", 		required_argument, 	0, 'x'},
        {"tag", 		required_argument, 	0, 't'},
        {"stat", 		no_argument, 		0, 's'},
        {"modifycv", 		required_argument, 	0, 'p'},
        {"mode", 		required_argument,	0, 'm'},
        {"data", 		required_argument,	0, 'd'},
	{"asimov",		required_argument,      0, 'a'},
        {"covariancematrix",    required_argument,      0, 'c'},
        {"geniematrix",         required_argument,      0, 'g'},
        {"flat",                required_argument,      0,'f'},
        {"edependent",          no_argument,            0,'e'},
        {"zeroout", 		no_argument,            0, 'z'},
        {"interpolation",       required_argument,      0, 'i'},
        //{"randomseed",        required_argument, 0 ,'r'},
        {"help", 		no_argument,	0, 'h'},
        {0,			    no_argument, 		0,  0},
    };

    int iarg = 0;
    opterr=1;
    int index;

    std::string xml = "Whatt.xml";
    //a tag to identify outputs and this specific run. defaults to EXAMPLE1
    std::string tag = "NoTag";
    std::string mode = "fit";
    int interpolation_number = -99;  //number of points for chi2 value interpolation
    //double random_number_seed = -1;

    bool bool_stat_only = false;
    bool bool_shape_fit = false;
    bool bool_zero_correlation = false;  //zero out correlation between subchannels
    bool bool_edependent = false;   //energy/momentum dependent fit or not
    bool input_data = false;
    bool bool_modify_genie_cv = false;   //modify genie CV before fitting
    bool bool_gen_sensitivity_curve = false;      //if the plot is intend to be a sensitivity curve
    std::string data_filename;
    std::string asimov_dataset;  //string containing information about the Asimov data
    std::string covmatrix_file;  //root file containing total covariance matrix
    std::string genie_matrix_file;  //root file containing flux/XS covariance matrix
    std::string det_matrix_file; //root file containing each det syst covar matrix;

    bool bool_flat_sys = false;
    double flat_sys_percent = 0.0;
    double delta_scaling = 1.0;

    while(iarg != -1)
    {
        iarg = getopt_long(argc,argv, "x:m:t:d:c:g:i:p:a:f:szeh", longopts, &index);

        switch(iarg)
        {
            case 'x':
                xml = optarg;
                break;
            case 't':
                tag = optarg;
                break;
            case 'd':
                input_data = true;
                data_filename = optarg;
                break;
            case 'm':
	    {
            	std::string plot_argument = optarg, delimiter=",";
	        mode = plot_argument.substr(0,plot_argument.find(delimiter));

		// if we are in plot mode, then we process the second argument
		if(mode == "plot" && (plot_argument.find(delimiter) != std::string::npos) && (plot_argument.substr(plot_argument.find(delimiter)+delimiter.length()) == "sensitivity"))
		    bool_gen_sensitivity_curve = true;	
		
                break;
	    }
            case 'c':
                covmatrix_file = optarg;
                break;
            case 'g':
		bool_shape_fit = true;
                genie_matrix_file = optarg;
                break;
            case 'f':
                bool_flat_sys = true;
                flat_sys_percent = (double)strtod(optarg,NULL);
                break;
                //case 'r':
                //    random_number_seed = (double)strtod(optarg,NULL);
                //    std::cout<<"Reading in random seed argument: "<<random_number_seed<<std::endl;
                //    break;
            case 'p':
                bool_modify_genie_cv= true;
                delta_scaling = (double)strtod(optarg,NULL);
                break;
            case 's':
                bool_stat_only = true;
                break;
 	    case 'z':
		bool_zero_correlation = true;
		break;
            case 'e':
                bool_edependent = true;
                break;
            case 'i':
                interpolation_number = (int)strtod(optarg, NULL);
                break;
	    case 'a':
		asimov_dataset = optarg;
		break;
            case '?':
            case 'h':
                std::cout<<"---------------------------------------------------"<<std::endl;
                std::cout<<"sbnfit_fraction_fit is a work in progress."<<std::endl;
                std::cout<<"---------------------------------------------------"<<std::endl;
                std::cout<<"--- Required arguments: ---"<<std::endl;
                std::cout<<"\t-x\t--xml\t\tInput configuration .xml file for SBNconfig"<<std::endl;
                std::cout<<"\t-d\t--data\t\tInput observed data for fit to real data"<<std::endl;
                std::cout<<"\t-m\t--mode\t\tInput running mode, default is 'fit' mode"<<std::endl;
                std::cout<<"\t-m\t--option\t\t gen --- generate energy/momentum dependent pre-scaling file "<<std::endl;
                std::cout<<"\t-m\t--option\t\t fit --- perform the fit to extract BF value "<<std::endl;
                std::cout<<"\t-m\t--option\t\t plot --- read files and draw plots "<<std::endl;
                std::cout<<"\t-m\t--option\t\t plot,sensitivity --- read files and draw plots with sensitivity banner overlaid "<<std::endl;
		std::cout<<"\t-m\t--option\t\t public --- print out event prediction, uncertainty for public data release" << std::endl;
		std::cout<<"\t-c\t--covariancematrix\t\tInput fractional covariance matrix"<< std::endl;
		std::cout<<"\t-g\t--geniematrix\t\tObsolete: Input FluxXS covariance matrix to remove genie uncertainty"<< std::endl;
                std::cout<<"--- Optional arguments: ---"<<std::endl;
		std::cout<<"\t-t\t--tag\t\t Name tag for output [default to: NoTag]" << std::endl;
                std::cout<<"\t-f\t--flat\t\t Input flat systematic fractional covariance matrix"<<std::endl;
		std::cout<<"\t-a\t--asimov\t\t Input configurations for Asimov dataset, in the form of, eg. NCDelta,3,NCpi0,2" << std::endl;
                std::cout<<"\t-s\t--stat\t\tStat only runs"<<std::endl;
                std::cout<<"\t-e\t--edependent\t\tRun in energy-dependent mode, for extraction of NCpi0 non-coherent"<<std::endl;
		std::cout<<"\t-z\t--zeroout\t\tZero out correlations between subchannels" << std::endl;
                std::cout<<"\t-i\t--interpolation\t\tInput number of points for interpolation"<< std::endl;
                std::cout<<"\t-p\t--modifycv\t\tModify CV: Input scaling factor for NCdelta, and scaling genie CV using result from NCpi0"<< std::endl;
                //std::cout<<"\t-r\t--randomseed\t\tRandomNumber Seed (default from machine)"<<std::endl; 
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

    std::cout<<"Begining Single Photon module"<<std::endl;

    NGrid mygrid, poly_grid;

    // flat fit to all subchannels
    //mygrid.AddConstrainedDimension("All", 0.5, 1.5, 0.01, 1.19);   //0.1 FULL
   // NCpi0 cos(theta_pi0) flat fit
   // TECHNOTE V6 version
    //mygrid.AddConstrainedDimension("NCPi0NotCoh", 0.5, 1.25, 0.01, 1.0);   //0.1 FULL
    //mygrid.AddConstrainedDimension("NCPi0Coh", 0, 5, 0.05, 1.0); //0.1full
    // NCpi0 flat normalization fit
    //mygrid.AddConstrainedDimension("NCPi0NotCoh", 0.2, 1.55, 0.01, 1.0);   //0.1 FULL
    //mygrid.AddConstrainedDimension("NCPi0Coh", 0, 7, 0.05, 1.0); //0.1full
    // NCpi0 momentum, momentum-dependent fit
    //mygrid.AddConstrainedDimension("NCPi0NotCoh", 0.5, 1.25, 0.02, 1.0);   //0.1 FULL
    //mygrid.AddConstrainedDimension("NCPi0NotCoh", 0.5, 1.25, 0.05, 1.0);   //coarser grid
    //mygrid.AddConstrainedDimension("NCPi0Coh", 0, 8, 0.2, 1.0); 
    //mygrid.AddFixedDimension("NCPi0NotCoh", 1.19);   //fixed
    mygrid.AddDimension("NCDelta", 0, 6, 3.18);
    //mygrid.AddDimension("NCDeltaLEE", 0, 5, 0.01 );
    //mygrid.AddDimension("NCDeltaLEE", 0, 2.5, 0.005 );


    poly_grid.AddConstrainedDimension("NCPi0NotCoh", -2.6, 1.0, 0.2, 1);  //zoomed in first order
    //poly_grid.AddConstrainedDimension("NCPi0NotCoh", -4.0, 2.0, 0.4, -1.1);  //first order
    //poly_grid.AddFixedDimension("NCPi0NotCoh", -1.05); // second order 

    if(mode == "gen"){
        if(bool_edependent){
            SBNsinglephoton sp(xml, tag, mygrid, poly_grid, true);
            sp.GeneratePreScaledSpectra();	
        }else{
            SBNsinglephoton sp(xml, tag, mygrid);
            sp.GeneratePreScaledSpectra();
        }
    }
    else if(mode =="fit"|| mode =="plot"|| mode == "calc"){
	SBNsinglephoton sp(xml, tag, mygrid);
	if(bool_edependent){   //setup poly grid if it's an energy-dependent fit
	    sp.SetPolyGrid(poly_grid);
	}

	if(mode == "fit"){
		//load CV, data and prescaled spectra
		sp.LoadCV();
		if(input_data) sp.LoadData(data_filename);
		else sp.SetupAsimovDataset(asimov_dataset);

		//setup systematic fractional covariance matrix
		if(bool_stat_only) sp.SetStatOnly();
		else if(bool_flat_sys) sp.SetFlatFullFracCovarianceMatrix(flat_sys_percent);
		//else  sp.SetFullFractionalCovarianceMatrix(covmatrix_file, "updated_frac_covariance");
		else  sp.SetFullFractionalCovarianceMatrix(covmatrix_file, "frac_covariance");

                /*
		//add another covariance matrix
		if( !bool_stat_only && bool_shape_fit){
		    sp.AddCovarianceMatrix(genie_matrix_file, "frac_covariance");
 	        }

		//zero out correlation if necessary
	   	if(bool_zero_correlation) sp.ZeroOutOffDiagonal();;
 		*/

	        //Obsolete
		//fit with normalization error removed  needs extra flux+XS syst covar matrix
		if(!bool_stat_only && bool_shape_fit){
		    sp.SetGenieFractionalCovarianceMatrix(genie_matrix_file);
		    sp.CalcFullButGenieFractionalCovarMatrix();
		    //sp.ZeroOutGenieCorrelation("NCDeltaLEE");
		}

		//Obsolete
		//if we want to modify NCpi0 to match the result from NCpi0 normalization fit before performing a combined fit
		if(bool_modify_genie_cv){
		  sp.ModifyCV(delta_scaling);
		  //sp.ModifyCV(delta_scaling, {1.0, 1.0});
		}


		sp.LoadSpectraApplyFullScaling(); 
		//sp.CalcChiGridScan();
		sp.CalcChiGridScanVaryMatrices();
	}else if(mode == "plot"){
		sp.GrabFitMap();
		if(interpolation_number != -99) sp.SetInterpolationNumber(interpolation_number);
		sp.SaveHistogram(bool_gen_sensitivity_curve);	
	}
	else{
		if(bool_modify_genie_cv) sp.ModifyCV(delta_scaling);
		double chi2=sp.CalcChi(true);
		std::cout << "SBNsinglephoton || chi2 value between (modified) CV and data is " << chi2 << std::endl;
	}
    }
    else if(mode == "public"){
	SBNsinglephoton sp(xml, tag);

        //load CV, data and prescaled spectra
        sp.LoadCV();
        sp.LoadData(data_filename);
        //setup systematic fractional covariance matrix
        sp.SetFullFractionalCovarianceMatrix(covmatrix_file, "frac_covariance");

	sp.PublicDataPrintOut({"NCDelta"}, {0});
	sp.PublicDataPrintOut({"NCDelta"}, {1});
	sp.PublicDataPrintOut({"NCDelta"}, {3.18});
    }
    else{
        std::cout << "Mode input is not identified, please try with a valid input.." << std::endl;
        std::cout << "Mode options: 'fit'  || Perform fit" << std::endl;
        std::cout << "Mode options: 'plot' || Grab info from root file and plot variables" << std::endl;
        std::cout << "Mode options: 'gen'  || Generate energy/momentum dependent pre-scaling root files" << std::endl;
        std::cout << "Mode options: 'calc' || Calculate the chi2 value of (corrected) CV and data" << std::endl;
	std::cout << "Mode options: 'public' || Do some simple prints out for public data release" << std::endl;
    }
    std::cout << "Single Photon module||" << "\tFinished" <<std::endl;
    std::cout << "Total wall time: " << difftime(time(0), start_time)/60.0 << " Minutes.\n";
    return 0;

}

