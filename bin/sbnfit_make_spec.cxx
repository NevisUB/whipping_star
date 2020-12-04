#include <iostream>
#include <iomanip>
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
	{"covarmatrix",         required_argument,      0, 'm'},
	{"nooscillate",         required_argument,      0, 'n'},
        {"printall", 		no_argument, 		0, 'p'},
	{"error",               required_argument,            0, 'e'},
        {"tag", 		required_argument,	0, 't'},
        {"help", 		no_argument,	0, 'h'},
        {0,			    no_argument, 		0,  0},
    };

    int iarg = 0;
    opterr=1;
    int index;
    bool compare_spec = false;
    bool print_error = false;
    bool do_oscillation = true;
    bool covar_matrix= false; // use covariance matrix or not

    //a tag to identify outputs and this specific run. defaults to EXAMPLE1
    std::string tag = "Central_Value"; //meaning Best Fit Point
    std::string data_file;  //data file name
    std::string ref_file;  //reference root file
    std::string covar_file; //covariance matrix root file

    while(iarg != -1)
    {
        iarg = getopt_long(argc,argv, "x:t:c:m:e:n:dph", longopts, &index);

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
	    case 'e':
		print_error = true;
		ref_file = optarg;
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


    if(print_error){
	    SBNspec cv_spec(ref_file, xml);
	    cv_spec.CalcFullVector();
	    cv_spec.CalcErrorVector();
	    
	    //loop over each channel plot
	    int overall_index = 0;
	    for(int im = 0; im <cv_spec.mode_names.size(); im++){
                for(int id = 0; id <cv_spec.detector_names.size(); id++){
                    for(int ic = 0; ic <cv_spec.channel_names.size(); ic++){

                         std::string canvas_name = cv_spec.mode_names.at(im)+"_"+cv_spec.detector_names.at(id)+"_"+cv_spec.channel_names.at(ic);
                         std::cout << "On Plot: " << canvas_name << std::endl;
			 //std::cout.width(15);  //set up the width of print out
			 std::cout <<  std::setw(15) << "Bin"; 
			 std::cout <<  std::setw(15) << "Total Error"; 
			 std::cout <<  std::setw(15) << "numu CC"; 
			 std::cout <<  std::setw(15) << "numu NC"; 
			 std::cout <<  std::setw(15) << "nue" << std::endl; 
			 for(int ibin = 0; ibin < cv_spec.num_bins[ic]; ibin++){
			      double total_bin_err = cv_spec.collapsed_err_vector.at(overall_index);
			      std::cout<<  std::setw(15) << ibin+1  <<  std::setw(15) << total_bin_err;
			      for(auto const& h: cv_spec.hist){
				  std::string test = h.GetName();
                                  if(test.find(canvas_name)!=std::string::npos ){	
					if(test.find("numu_cc")!=std::string::npos)
						std::cout << std::setw(15)<<h.GetBinError(ibin+1)/total_bin_err;
					if(test.find("numu_nc")!=std::string::npos)
						std::cout << std::setw(15)<<h.GetBinError(ibin+1)/total_bin_err;
					if(test.find("numu_nue")!=std::string::npos){
						std::cout << std::setw(15)<<h.GetBinError(ibin+1)/total_bin_err << std::endl;
					}
				  }
	        	      }

			      overall_index++;
			 }
		     }
		}
	    }

    }
    else if(!compare_spec){
	    std::cout<<"Begining building SBNspec for tag: "<<tag<<std::endl;
	    //now only using gen(xml, NeutrinoModel) can avoid closing SBNgenerate before writing spectrums.
	    //NeutrinoModel nullModel(0,0,0);
	    NeutrinoModel nullModel(sqrt(1.32), pow(10, 0.0), sqrt((1+sqrt(1-7e-2))/2));

	    //initialize SBNgenerate, which will generate SBNspec and fill the hisotgrams
	    SBNgenerate gen_cv(xml, nullModel);

	    //write out the SBNspec in root files
	    //gen_cv.WriteCVSpec(tag);
	    gen_cv.WritePrecomputedOscSpecs("NuMuDis");
    }
    else{
	SBNspec data_spec(data_file, xml);
	if(do_oscillation == false){
	    SBNspec ref_spec(ref_file, xml);
	    ref_spec.CalcFullVector();
	    ref_spec.CalcErrorVector();

        	//std::cout << "check 1" << std::endl;
	    //ref_spec.ScaleAll(5.23537/66);

	    if(covar_matrix){
		TFile* f_cov = new TFile(covar_file.c_str(), "read");
		TMatrixT<double>* p_covar = (TMatrixT<double>*)f_cov->Get("frac_covariance");
		TMatrixT<double> full_covar(ref_spec.num_bins_total, ref_spec.num_bins_total);
		TMatrixT<double> collapse_covar(ref_spec.num_bins_total_compressed, ref_spec.num_bins_total_compressed);
		SBNchi chi_temp(xml);
		
                full_covar = chi_temp.FillSystMatrix(p_covar, ref_spec.full_vector, ref_spec.full_err_vector);  //systematic covar matrix only
		//full_covar = chi_temp.CalcCovarianceMatrix(p_covar, ref_spec.full_vector);
		chi_temp.CollapseModes(full_covar, collapse_covar);
		ref_spec.CompareSBNspecs(collapse_covar, &data_spec, tag);
	    }
	    else ref_spec.CompareSBNspecs(&data_spec, tag);
	}
	else{
	   //if we want to oscillate the spectrum and compare it with certain data spectrum
	    //NeutrinoModel oscModel(pow(10, dm grid), pow(10, e4 grid), pow(10, u4 grid));
	    // numu disappearance
	    NeutrinoModel oscModel(sqrt(1.32), pow(10, 0.0), sqrt((1+sqrt(1-7e-2))/2));
	    // nue appearance
	    //NeutrinoModel oscModel(sqrt(1.20), sqrt(3e-3)/2, 1.0);
	    // nue disappearance
	    //NeutrinoModel oscModel(sqrt(1.32), sqrt((1+sqrt(1-7e-2))/2), 1.0);

	    //oscillaton based on the model
	    SBNgenerate gen_osc(xml, oscModel);

	    //write out pre-oscillated spectrum
	    gen_osc.WritePrecomputedOscSpecs(tag);
	    //construct background spectrum
	    //gen_osc.spec_central_value.Scale("fullosc",0.0);
	    //gen_osc.spec_central_value.WriteOut(tag+"_BKG_ONLY");

	    //calculate fully oscillated spectrum
	    SBNosc osc(tag+"_BKG_ONLY.SBNspec.root", xml);
	    //std::cout << "check " << __LINE__ <<std::endl;
	    osc.LoadModel(oscModel);   
	    //std::cout << "check " << __LINE__ <<std::endl;
	    std::vector<std::vector<double>> ans = osc.Oscillate(tag, false);
	    //std::cout << "check " << __LINE__ <<std::endl;
	    SBNspec osc_spec(ans[0], ans[1], xml);  //oscillated SBNspec
	    //std::cout << "check " << __LINE__ <<std::endl;
	    //osc_spec.ScaleAll(5.23537/66.0);
	    //osc_spec.ScaleAll(66.0/5.81731);
	    //data_spec.ScaleAll(5.23537/66.0);
	    osc_spec.WriteOut("InjectedPoint_sinsq_7e-2_dmsq_1.32");
	    tag = "ExampleOscillatedSpectra_vs_Data";
	   /* if(covar_matrix){
                TFile* f_cov = new TFile(covar_file.c_str(), "read");
                TMatrixT<double>* p_covar = (TMatrixT<double>*)f_cov->Get("frac_covariance");
                TMatrixT<double> full_covar(osc_spec.num_bins_total, osc_spec.num_bins_total);
                TMatrixT<double> collapse_covar(osc_spec.num_bins_total_compressed, osc_spec.num_bins_total_compressed);
                SBNchi chi_temp(xml);

                full_covar = chi_temp.FillSystMatrix(p_covar, osc_spec.full_vector, osc_spec.full_err_vector);  //systematic covar matrix only
                //full_covar = chi_temp.CalcCovarianceMatrix(p_covar, osc_spec.full_vector);   //full stat+syst covar matrix
                chi_temp.CollapseModes(full_covar, collapse_covar);
                osc_spec.CompareSBNspecs(collapse_covar, &data_spec, tag);
            }
            //else osc_spec.CompareSBNspecs(&data_spec, tag);
            else data_spec.CompareSBNspecs(&osc_spec, tag);
	  */ 
	}
    }

    std::cout << "Total wall time: " << difftime(time(0), start_time)/60.0 << " Minutes.\n";
    return 0;

}
