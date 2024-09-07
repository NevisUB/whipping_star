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
	{"reference",           required_argument,      0, 'r'},
        {"printpdf", 		required_argument, 	0, 'p'},
	{"m4", 		        required_argument, 	0, 'i'},
        {"ue4", 		required_argument, 	0, 'j'},
        {"um4", 		required_argument, 	0, 'k'},
        {"channel", 		required_argument, 	0, 'o'},
        {"error",               required_argument,      0, 'e'},
        {"tag", 		required_argument,	0, 't'},
        {"help", 		no_argument,	0, 'h'},
        {0,			    no_argument, 		0,  0},
    };

    int iarg = 0;
    opterr=1;
    int index;
    bool compare_spec = false;
    bool print_error = false;
    bool covar_matrix= false; // use covariance matrix or not
    double m4=0;
    double ue4=0;    
    double um4=0;
    //a tag to identify outputs and this specific run. defaults to EXAMPLE1
    std::string tag = "Central_Value"; //meaning Best Fit Point
    std::string data_file;  //data file name
    std::string ref_file;  //reference root file
    std::string covar_file; //covariance matrix root file
    std::string channel = "null";

    while(iarg != -1)
    {
        iarg = getopt_long(argc,argv, "x:t:c:m:e:r:p:dh", longopts, &index);

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
	    case 'r':
		ref_file = optarg;
		break;
            case 'p':
                print_mode=true;
		ref_file = optarg;
                break;
	    case 'e':
		print_error = true;
		ref_file = optarg;
		break;
            case 't':
                tag = optarg;
                break;
            case 'i':
                m4 = (double)strtod(optarg,NULL);
                break;
            case 'j':
                ue4 = (double)strtod(optarg,NULL);
                break;
            case 'k':
                um4 = (double)strtod(optarg,NULL);
                break;
            case 'o':
                channel = optarg;
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

    if(print_mode){
 	    SBNspec cv_spec(ref_file, xml);
            cv_spec.CalcFullVector();
            cv_spec.CalcErrorVector();
                
    //	    if(covar_matrix){
//	        TMatrixT<double> collapse_covar(cv_spec.num_bins_total_compressed, cv_spec.num_bins_total_compressed);
  //              TFile* f_cov = new TFile(covar_file.c_str(), "read");
    //            TMatrixT<double>* p_covar = (TMatrixT<double>*)f_cov->Get("frac_covariance");
      //          TMatrixT<double> full_covar(cv_spec.num_bins_total, cv_spec.num_bins_total);
        //        SBNchi chi_temp(xml);

          //      full_covar = chi_temp.FillSystMatrix(p_covar, cv_spec.full_vector, cv_spec.full_err_vector);  //systematic covar matrix only
            //    chi_temp.CollapseModes(full_covar, collapse_covar);
	//	f_cov->Close();
		//cv_spec.DrawSpecs(tag, &collapse_covar);
       //    }else{
		//cv_spec.DrawSpecs(tag);
    }
    else if(print_error){
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
	 
   	    std::string write_out_tag;
 
	    //if we want to oscillate the spectrum and compare it with certain data spectrum
	    //NeutrinoModel oscModel(pow(10, dm grid), pow(10, e4 grid), pow(10, u4 grid));
	    //**** numu disappearance****
	    //NeutrinoModel oscModel(sqrt(1.32), pow(10, 0.0), sqrt((1+sqrt(1-7e-2))/2));
 	    //write_out_tag = "InjectedPoint_NumuDisappear_sinsq_7e-2_dmsq_1.32";
	    //**** nue appearance ****
	    //NeutrinoModel oscModel(sqrt(1.32), sqrt(3e-3)/2, 1.0);
 	    //write_out_tag = "InjectedPoint_NueAppear_sinsq_3e-3_dmsq_1.32";
	    //**** nue disappearance ****
            NeutrinoModel oscModel(std::pow(10.0, m4), std::pow(10.0, ue4), std::pow(10.0, um4));
            //NeutrinoModel oscModel(sqrt(3), sqrt((1+sqrt(1-0.4))/2), 1.0);
 	    write_out_tag = "spectrum_"+std::to_string(m4)+"_"+std::to_string(ue4)+"_"+std::to_string(um4)+"_"+channel;

	    // for direct accessible Ue4, Uu4.
	    //NeutrinoModel oscModel(dm, Ue4, Uu4);
	    //oscillaton based on the model
	    SBNgenerate gen_osc(xml, oscModel);

	    //write out pre-oscillated spectrum
	    gen_osc.WritePrecomputedOscSpecs(tag);
	    //construct background spectrum
	    //gen_osc.spec_central_value.Scale("fullosc",0.0);
	    //gen_osc.spec_central_value.WriteOut(tag+"_BKG_ONLY");

	    //calculate fully oscillated spectrum
	    SBNosc osc(tag+"_BKG_ONLY.SBNspec.root", xml);
	    osc.LoadModel(oscModel);   
	    //std::cout << "check " << __LINE__ <<std::endl;
	    std::vector<std::vector<double>> ans = osc.Oscillate(tag, false);
	    SBNspec osc_spec(ans[0], ans[1], xml);  //oscillated SBNspec
	    //osc_spec.ScaleAll(5.23537/66.0);
	    //osc_spec.ScaleAll(66.0/5.81731);
	    //data_spec.ScaleAll(5.23537/66.0);
	    osc_spec.WriteOut(write_out_tag);
    }
    else{
	    SBNspec data_spec(data_file, xml);
	    SBNspec ref_spec(ref_file, xml);
	    ref_spec.CollapseVector();
	    ref_spec.CalcErrorVector();

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
		for(int i = 0; i != collapse_covar.GetNcols(); ++i){
		    std::cout << "bin: " << i << ", syst error: " <<100*sqrt(collapse_covar(i,i))/ref_spec.collapsed_vector[i] << "%, stats error: " << 100/sqrt(ref_spec.collapsed_vector[i]) <<"%" << std::endl;
		}
//		ref_spec.CompareSBNspecs(collapse_covar, &data_spec, tag);
	    }
	    else ref_spec.CompareSBNspecs(&data_spec, tag);
            //else data_spec.CompareSBNspecs(&ref_spec, tag);

    }

    std::cout << "Total wall time: " << difftime(time(0), start_time)/60.0 << " Minutes.\n";
    return 0;

}
