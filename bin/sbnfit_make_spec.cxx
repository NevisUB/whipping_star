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
#include "prob.h"
#include "SBNconfig.h"
#include "SBNchi.h"
#include "SBNspec.h"
#include "SBNosc.h"
#include "SBNfit.h"
#include "SBNfit3pN.h"
#include "SBNgenerate.h"

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

    std::string xml = "example.xml";
    bool print_mode = false;

    /*************************************************************
     *************************************************************
     *		Command Line Argument Reading
     ************************************************************
     ************************************************************/
    const struct option longopts[] =
    {
        {"xml", 		required_argument, 	0, 'x'},
	{"compare",             required_argument,      0, 'c'},
	{"covarmatrix",             required_argument,      0, 'm'},
	{"nooscillate",             required_argument,      0, 'n'},
        {"printall", 		no_argument, 		0, 'p'},
        {"tag", 		required_argument,	0, 't'},
        {"help", 		no_argument,	0, 'h'},
        {0,			    no_argument, 		0,  0},
    };

    int iarg = 0;
    opterr=1;
    int index;
    bool compare_spec = false;
    bool do_oscillation = true;
    bool covar_matrix= false; // use covariance matrix or not

    //a tag to identify outputs and this specific run. defaults to EXAMPLE1
    std::string tag = "Central_Value"; //meaning Best Fit Point
    std::string data_file;  //data file name
    std::string ref_file;  //reference root file
    std::string covar_file; //covariance matrix root file

    while(iarg != -1)
    {
        iarg = getopt_long(argc,argv, "x:t:c:m:n:dph", longopts, &index);

        switch(iarg)
        {
            case 'x':
                xml = optarg;
                break;
	    case 'c':
		compare_spec = true;
	        data_file= optarg;
		break;
	    case 'm':
		covar_matrix= true;
		covar_file = optarg;
		break;
	    case 'n':
		do_oscillation = false;
		ref_file = optarg;
		break;
            case 'p':
                print_mode=true;
                break;
            case 't':
                tag = optarg;
                break;
            case '?':
            case 'h':
                std::cout<<"---------------------------------------------------"<<std::endl;
                std::cout<<"sbnfit_make_covariance allows for the building of covariance matricies from input root files containing reconstructed variables and the EventWeight class std::map<std::string, std::vector<double>>."<<std::endl;
                std::cout<<"---------------------------------------------------"<<std::endl;
                std::cout<<"--- Required arguments: ---"<<std::endl;
                std::cout<<"\t-x\t--xml\t\tInput configuration .xml file for SBNconfig"<<std::endl;
                std::cout<<"\t-t\t--tag\t\tA unique tag to identify the outputs [Default to BF]"<<std::endl;
                std::cout<<"--- Optional arguments: ---"<<std::endl;
                std::cout<<"\t-p\t--printall\tRuns in BONUS print mode, making individual spectra plots for ALLVariations. (warning can take a while!) "<<std::endl;
                std::cout<<"\t-h\t--help\t\tThis help menu."<<std::endl;
                std::cout<<"---------------------------------------------------"<<std::endl;

                return 0;
        }
    }




    //std::string dict_location = "../libio/libEventWeight.so";
    //std::cout<<"Trying to load dictionary: "<<dict_location<<std::endl;
    //gSystem->Load(  (dict_location).c_str());

    /*************************************************************
     *************************************************************
     *			Main Program Flow
     ************************************************************
     ************************************************************/
    time_t start_time = time(0);


    if(!compare_spec){
	    std::cout<<"Begining building SBNspec for tag: "<<tag<<std::endl;
	    //now only using gen(xml, NeutrinoModel) can avoid closing SBNgenerate before writing spectrums.
	    NeutrinoModel nullModel(0,0,0);
	    //NeutrinoModel nullModel(sqrt(2), pow(10, 0.0), sqrt((1+sqrt(1-3e-3))/2));

	    //initialize SBNgenerate, which will generate SBNspec and fill the hisotgrams
	    SBNgenerate gen_cv(xml, nullModel);

	    //write out the SBNspec in root files
	    gen_cv.WriteCVSpec(tag);
	    //gen_cv.WritePrecomputedOscSpecs("NuMuDis");
    }
    else{
	SBNspec data_spec(data_file, xml);
	if(do_oscillation == false){
	    SBNspec ref_spec(ref_file, xml);
	    ref_spec.CalcFullVector();

        	//std::cout << "check 1" << std::endl;
	    //ref_spec.ScaleAll(5.23537/66);

	    if(covar_matrix){
		TFile* f_cov = new TFile(covar_file.c_str(), "read");
		TMatrixT<double>* p_covar = (TMatrixT<double>*)f_cov->Get("frac_covariance");
		TMatrixT<double> full_covar(ref_spec.num_bins_total, ref_spec.num_bins_total);
		TMatrixT<double> collapse_covar(ref_spec.num_bins_total_compressed, ref_spec.num_bins_total_compressed);
		SBNchi chi_temp(xml);
		
                full_covar = chi_temp.FillSystMatrix(p_covar, ref_spec.full_vector);  //systematic covar matrix only
		//full_covar = chi_temp.CalcCovarianceMatrix(p_covar, ref_spec.full_vector);
		chi_temp.CollapseModes(full_covar, collapse_covar);
		ref_spec.CompareSBNspecs(collapse_covar, &data_spec, tag);
	    }
	    else ref_spec.CompareSBNspecs(&data_spec, tag);
	}
	else{
	   //if we want to oscillate the spectrum and compare it with certain data spectrum
	    NeutrinoModel oscModel(pow(10, 0.15), pow(10, 0.0), pow(10, -0.15));
	    //NeutrinoModel oscModel(sqrt(2), pow(10, 0.0), sqrt((1+sqrt(1-3e-3))/2));
	    //NeutrinoModel oscModel(pow(10, 0.0603), pow(10, 0.0),pow(10, -0.8));
	    //NeutrinoModel oscModel(pow(10, 0.1), pow(10, 0.0),pow(10, -0.925));

	    //oscillaton based on the model
	    //SBNgenerate gen_osc(xml, oscModel);

	    //write out pre-oscillated spectrum
	    //gen_osc.WritePrecomputedOscSpecs(tag);

	    //construct background spectrum
	    //gen_osc.spec_central_value.Scale("fullosc",0.0);
	    //gen_osc.spec_central_value.WriteOut(tag+"_BKG_ONLY");

	    //calculate fully oscillated spectrum
	    SBNosc osc(tag+"_BKG_ONLY.SBNspec.root", xml);
	    osc.LoadModel(oscModel);   
	    std::vector<double> ans = osc.Oscillate(tag, false);
	    SBNspec osc_spec(ans, xml);  //oscillated SBNspec
	    osc_spec.ScaleAll(5.23537/66.0);
	    //osc_spec.ScaleAll(66.0/5.81731);
	    //data_spec.ScaleAll(66.0/5.81731);
	    //osc_spec.WriteOut("Oscillation_sinsq_3e-3_dmsq_2");
	    tag = "BF_vs_Data";
	    if(covar_matrix){
                TFile* f_cov = new TFile(covar_file.c_str(), "read");
                TMatrixT<double>* p_covar = (TMatrixT<double>*)f_cov->Get("frac_covariance");
                TMatrixT<double> full_covar(osc_spec.num_bins_total, osc_spec.num_bins_total);
                TMatrixT<double> collapse_covar(osc_spec.num_bins_total_compressed, osc_spec.num_bins_total_compressed);
                SBNchi chi_temp(xml);

                full_covar = chi_temp.FillSystMatrix(p_covar, osc_spec.full_vector);  //systematic covar matrix only
                //full_covar = chi_temp.CalcCovarianceMatrix(p_covar, osc_spec.full_vector);   //full stat+syst covar matrix
                chi_temp.CollapseModes(full_covar, collapse_covar);
                osc_spec.CompareSBNspecs(collapse_covar, &data_spec, tag);
            }
            else osc_spec.CompareSBNspecs(&data_spec, tag);
            //else data_spec.CompareSBNspecs(&osc_spec, tag);
	  
	}
    }

    std::cout << "Total wall time: " << difftime(time(0), start_time)/60.0 << " Minutes.\n";
    return 0;

}
