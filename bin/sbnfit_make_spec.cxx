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
	{"mcfile", required_argument, 0, 'r'},
	{"covarmatrix",             required_argument,      0, 'm'},
        {"printall", 		no_argument, 		0, 'p'},
        {"tag", 		required_argument,	0, 't'},
        {"help", 		no_argument,	0, 'h'},
        {0,			    no_argument, 		0,  0},
    };

    int iarg = 0;
    opterr=1;
    int index;
    bool compare_spec = false;
    bool covar_matrix= false; // use covariance matrix or not

    //a tag to identify outputs and this specific run. defaults to EXAMPLE1
    std::string tag = "Central_Value"; //meaning Best Fit Point
    std::string data_file;  //data file name
    std::string ref_file;  //reference root file
    std::string covar_file; //covariance matrix root file

    while(iarg != -1)
    {
        iarg = getopt_long(argc,argv, "x:t:c:m:n:r:dph", longopts, &index);

        switch(iarg)
        {
            case 'x':
                xml = optarg;
                break;
	    case 'c':
		compare_spec = true;
	        data_file= optarg;
		break;
	    case 'r':
		ref_file=optarg;
		break;
	    case 'm':
		covar_matrix= true;
		covar_file = optarg;
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




     /*************************************************************
     *			Main Program Flow
     ************************************************************
     ************************************************************/
    time_t start_time = time(0);


    if(!compare_spec){
	    std::cout<<"Begining building SBNspec for tag: "<<tag<<std::endl;
	    //now only using gen(xml, NeutrinoModel) can avoid closing SBNgenerate before writing spectrums.
	    NeutrinoModel nullModel(0,0,0);

	    //initialize SBNgenerate, which will generate SBNspec and fill the hisotgrams
	    SBNgenerate gen_cv(xml, nullModel);

	    //write out the SBNspec in root files
	    gen_cv.WriteCVSpec(tag);
    }
    else{
	SBNspec data_spec(data_file, xml);
	SBNspec ref_spec(ref_file, xml);
	//ref_spec.Scale("NCDeltaRadOverlayLEE", 2.3);
	//ref_spec.Scale("NCPi0NotCoh", 1.1);
	//ref_spec.Scale("NCPi0Coh", 2.8);
	ref_spec.CollapseVector();
	ref_spec.CalcErrorVector();
	data_spec.CollapseVector();

        	//std::cout << "check 1" << std::endl;
	SBNchi chi_temp(xml);
	chi_temp.is_stat_only = false;

	    if(covar_matrix){
		TFile* f_cov = new TFile(covar_file.c_str(), "read");
		TMatrixT<double>* p_covar = (TMatrixT<double>*)f_cov->Get("frac_covariance");
		//TMatrixT<double> full_covar(ref_spec.num_bins_total, ref_spec.num_bins_total);
		TMatrixT<double> collapse_covar(ref_spec.num_bins_total_compressed, ref_spec.num_bins_total_compressed);
		
                collapse_covar = chi_temp.FillSystMatrix(*p_covar, ref_spec.full_vector, ref_spec.full_err_vector, true);  //systematic covar matrix only
		//full_covar = chi_temp.CalcCovarianceMatrix(p_covar, ref_spec.full_vector);
		//chi_temp.CollapseModes(full_covar, collapse_covar);
		//calculate chi2 between data and MC
		TMatrixT<double> full_covar = chi_temp.AddStatMatrixCNP(&collapse_covar, ref_spec.collapsed_vector, data_spec.collapsed_vector);
	    	TMatrixT<double> inverted_covar = chi_temp.InvertMatrix(full_covar);
		double chi = chi_temp.CalcChi(inverted_covar, ref_spec.collapsed_vector, data_spec.collapsed_vector);
		std::cout << "Chi2 between data and MC is " << chi << std::endl;
		chi_temp.DrawComparisonIndividual(ref_spec, data_spec, collapse_covar, tag, true);
	    }
	    else chi_temp.DrawComparisonIndividual(ref_spec, data_spec, tag);
   }
    

    std::cout << "Total wall time: " << difftime(time(0), start_time)/60.0 << " Minutes.\n";
    return 0;

}
