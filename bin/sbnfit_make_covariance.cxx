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
        {"detsys",              no_argument,            0, 'd'},
        {"selectchannel",       required_argument,      0, 's'},
    	{"covar",               required_argument,      0, 'c'},
        {"printall", 		no_argument, 		0, 'p'},
        {"tag", 		required_argument,	0, 't'},
        {"minimal",     no_argument,    0, 'm'},
        {"writesignal",         required_argument,      0, 'w'},
        {"help", 		no_argument,	0, 'h'},
    	{"constrain",		no_argument,    0, 'c'},
        {0,			    no_argument, 		0,  0},
    };

    int iarg = 0;
    opterr=1;
    int index;
    bool bool_use_universe = true;
    bool use_existing_covar=false;
    bool constrain_mode=false;
    bool minimal_matrix_plots = false;

    //a tag to identify outputs and this specific run. defaults to EXAMPLE1
    std::string tag = "DEFAULT_TAG";
    std::string signal_def = "DEFAULT_SIGNAL_DEF";
    std::string comma_separated_covar_info = "NONE";
    std::string comma_separated_channels;

    while(iarg != -1)
    {
        iarg = getopt_long(argc,argv, "x:s:c:t:w:mdph", longopts, &index);

        switch(iarg)
        {
            case 'x':
                xml = optarg;
                break;
            case 'p':
                print_mode=true;
                break;
            case 'd':
                bool_use_universe=false;
                break;
	    case 's':
                comma_separated_channels = optarg; 
		break;
            case 't':
                tag = optarg;
                break;
        case 'm':
                minimal_matrix_plots = true;
                break;
	    case 'c':
	      constrain_mode=true;
	      break;
	    case 'c':
	        use_existing_covar=true;
		comma_separated_covar_info = optarg;
		break;
	    case 'w':
		signal_def = optarg;
		break;
            case '?':
            case 'h':
                std::cout<<"---------------------------------------------------"<<std::endl;
                std::cout<<"sbnfit_make_covariance allows for the building of covariance matricies from input root files containing reconstructed variables and the EventWeight class std::map<std::string, std::vector<double>>."<<std::endl;
                std::cout<<"---------------------------------------------------"<<std::endl;
                std::cout<<"--- Required arguments: ---"<<std::endl;
                std::cout<<"\t-x\t--xml\t\tInput configuration .xml file for SBNconfig"<<std::endl;
                std::cout<<"\t-t\t--tag\t\tA unique tag to identify the outputs [Default to DEFAULT_TAG]"<<std::endl;
		std::cout<<"\t-c\t--covar\t\tUse with option -s.\tInput: covariance root file,matrix name. Given selected channels, generate a new covariance matrix for these channels from submatrices of the given matrix"<< std::endl;
		std::cout<<"\t-s\t--selectchannel\t\tUse together with option -c. \tInput: comma separated channel names" << std::endl;
		std::cout<<"\t-w\t--writeout\t\tWrite out spectra with signal and background separated. \tInput: signal definition - any subchannel whose name include signal definition will be summed into signal, otherwise they will be added into background" << std::endl;
                std::cout<<"--- Optional arguments: ---"<<std::endl;
                std::cout<<"\t-d\t--detsys\t use root files with systematically varied histograms (detsys) to build the covariance matrix" << std::endl;
                std::cout<<"\t-p\t--printall\tRuns in BONUS print mode, making individual spectra plots for ALLVariations. (warning can take a while!) "<<std::endl;
                std::cout<<"\t-m\t--minimal\tDoesn't print individual variations in matrix_plots"<<std::endl;
                std::cout<<"\t-h\t--help\t\tThis help menu."<<std::endl;
                std::cout<<"---------------------------------------------------"<<std::endl;

                return 0;
        }
    }

    if(print_mode==true && !bool_use_universe){
        std::cout<<"Error! Sorry, cant use printall and detsys together yet!"<<std::endl;
        return 1;
    }

    gSystem->Load("/uboone/app/users/markrl/Hive_v3.0/hellstroms_hive/hive/root_linkdefs/loc/SL7/denan_cxx.so");

    //std::string dict_location = "../libio/libEventWeight.so";
    //std::cout<<"Trying to load dictionary: "<<dict_location<<std::endl;
    //gSystem->Load(  (dict_location).c_str());

    /*************************************************************
     *************************************************************
     *			Main Program Flow
     ************************************************************
     ************************************************************/
    time_t start_time = time(0);

    std::cout<<"Begining Covariance Calculation for tag: "<<tag<<std::endl;

    
    if(use_existing_covar){
	//grab submatrices from existing covariance matrix

	SBNcovariance example_covar(xml, tag);	


	//grab input covariance matrix info
	std::string delimiter=",";
	std::string input_file = comma_separated_covar_info.substr(0, comma_separated_covar_info.find(delimiter));
	std::string matrix_name = comma_separated_covar_info.substr(comma_separated_covar_info.find(delimiter)+delimiter.length());

	//grab the channel info
	std::cout << comma_separated_channels << std::endl;
	std::vector<std::string> channels;
	size_t pos_front = 0, pos_end = comma_separated_channels.find(delimiter);
	while(pos_end != std::string::npos){
	    channels.push_back(comma_separated_channels.substr(pos_front, pos_end-pos_front));
	    pos_front = pos_end + delimiter.length();
	    pos_end = comma_separated_channels.find(delimiter, pos_front);
	}
	channels.push_back(comma_separated_channels.substr(pos_front));

	//form submatrix and save!
	example_covar.GrabSubMatrix(input_file, matrix_name, channels);
    }
    //Create a SBNcovariance object initilizing with the inputted xml
    //This will load all the files and weights as laid out
    else if(bool_use_universe){
        SBNcovariance example_covar(xml);

	//Write out variations in histograms
	example_covar.WriteOutVariation(signal_def);

        //Form the covariance matrix from loaded weights and MC events
        example_covar.FormCovarianceMatrix(tag);

        //and make some plots of the resulting things
        //Will be outputted in the form: SBNfit_covariance_plots_TAG.root
        example_covar.PrintMatricies(tag,minimal_matrix_plots);

        //Constraint will be patched in shortly: mark
	std::ofstream ratio_con;
	std::vector<double> average_ratio;
	 ratio_con.open("ratio_list_"+tag+".txt",std::ofstream::trunc);
	 if(constrain_mode){  
           for (int i=0;i<example_covar.variations.size();i++){
	     std::cout<<"var " <<example_covar.variations.at(i)<<std::endl;
	     average_ratio=example_covar.DoConstraint(0,1,tag,i);
	     //ratio_con<<i<<" " <<average_ratio<<std::endl;
	     ratio_con<<"var "<<i<<" name "<<example_covar.variations.at(i)<<" events "<<average_ratio[0]<<" uncon "<<average_ratio[1]<<" con "<<average_ratio[2]<<" ratio "<<average_ratio[3]<<std::endl;
        }
	   
	   example_covar.DoConstraint_test(0,1,tag);
	 }  
        if(print_mode){
            //This takes a good bit longer, and prints every variation to file. 
            example_covar.PrintVariations(tag);
            example_covar.PrintVariations_2D(tag);
        }

    }else{
        SBNcovariance example_covar(xml, bool_use_universe);

	//Write out variations in histograms
	example_covar.WriteOutVariation(signal_def);
	
        //Form the covariance matrix from loaded weights and MC events
        example_covar.FormCovarianceMatrix(tag);

        //and make some plots of the resulting things
        //Will be outputted in the form: SBNfit_covariance_plots_TAG.root
        example_covar.PrintMatricies(tag,minimal_matrix_plots);

        //Constraint will be patched in shortly: mark
        //example_covar.DoConstraint(0,1);

        if(print_mode){
            //This takes a good bit longer, and prints every variation to file. 
            example_covar.PrintVariations(tag);
        }
    }


    std::cout << "Total wall time: " << difftime(time(0), start_time)/60.0 << " Minutes.\n";
    return 0;

}
