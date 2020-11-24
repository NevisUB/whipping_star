#ifndef SBNCONFIG_H_
#define SBNCONFIG_H_

#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <map>
#include <time.h>
#include <numeric>
#include <algorithm>
#include <iomanip>
#include "TH1D.h"
#include "TFile.h"
#include "TTreeFormula.h"

#include "tinyxml.h"
#include "branch_variable.h"

template <typename T>
std::string to_string_prec(const T a_value, const int n = 6)
{
	std::ostringstream out;
	out <<std::fixed<< std::setprecision(n) << a_value;
	return out.str();
}
//#define TYPE_FLOAT
#ifdef TYPE_FLOAT  
    typedef float eweight_type;
#else
    typedef double eweight_type;
#endif

namespace sbn{
// All declarations are within the namespace scope.

//Order is important, "xml order" means that we loop over all modes,detectors channels and subchannels, But finished each mode,det and channels before moving on.
// if we have modes m1 and m2, detectors d1 and d2, channels c1 and c2 and subchannels c1s1 c1s2 and c2s1 then the order of spectra would be..
// m1_d1_c1_c1s1
// m1_d1_c1_c1s2
// m1_d1_c2_c2s1
// m1_d2_c1_c1s1
// m1_d2_c1_c1s2
// m1_d2_c2_c2s1

// m2_d1_c1_c1s1
// m2_d1_c1_c1s2
// m2_d1_c2_c2s1
// m2_d2_c1_c1s1
// m2_d2_c1_c1s2
// m2_d2_c2_c2s1

class SBNconfig {

	protected:
	
	public:
	
	//Constructors
	SBNconfig(std::string, bool, bool); //read xml and do configuration. first 'bool': verbose or not, second 'bool': is the eventweights of different universes used to build the covariance matrix, or do we feed into root file with histograms from different systematic variations to build the covariance matrix.
	SBNconfig(std::string,bool);
	SBNconfig(std::string);
	SBNconfig(){};
	SBNconfig(std::vector<std::string>, std::vector<std::string>, std::vector<std::string>, std::vector<std::vector<std::string>>, std::vector<std::vector<double>>);
	//This is going to be a manual Setup thing

	//Some stringsteam stuff
	std::string otag;

	// Fullnames is kinda important, it contains all the concatanated names of all individual histograms that have been configured with the "use=1" attribute
	// The order is IMPORTANT its the same as defined in xml
	std::vector<std::string> fullnames;


	//Bools to contain what is and is not in the xml
	bool has_oscillation_patterns;
	bool is_verbose;
	bool use_universe;

	int num_detectors;
	int num_detectors_xml;
	int num_channels;
	int num_channels_xml;
	int num_modes;
	int num_modes_xml;

    double plot_pot = 1.0;

	//vectors of length num_channels
	std::vector<int> num_subchannels; 
	int* a_num_subchannels;
	std::vector<int> num_subchannels_xml; 
	std::vector<int> num_bins;
	int* a_num_bins;

	std::string xmlname;	

	std::string data_path;

	int num_bins_detector_block;
	int num_bins_mode_block;
	int num_bins_total;

	int num_bins_detector_block_compressed;
	int num_bins_mode_block_compressed;
	int num_bins_total_compressed;

	std::string correlation_matrix_rootfile;
	std::string correlation_matrix_name;

	
	//the xml names are the way we track which channels and subchannels we want to use later
	std::vector<std::string> mode_names; 			
	std::vector<std::string> detector_names; 		
	std::vector<std::string> channel_names; 		
	std::vector<std::string> channel_units; 		
	std::vector<std::vector<std::string >> subchannel_names; 

    std::vector<std::string> mode_plotnames; 			
	std::vector<std::string> detector_plotnames; 		
	std::vector<std::string> channel_plotnames; 		
	std::vector<std::vector<std::string >> subchannel_plotnames; 

    


	// vector Bools for turning on and off certain modes/detectors..etc..
	std::vector<bool> mode_bool; 
	std::vector<bool> detector_bool; 
	std::vector<bool> channel_bool; 
	std::vector<std::vector<bool >> subchannel_bool; 

	std::vector<int> channel_used;
	std::vector<int> detector_used;
	std::vector<int> mode_used;
	std::vector<std::vector<int>> subchannel_used;

	//An oscillation pattern, for oscillations
	std::vector<std::vector<int> > subchannel_osc_patterns; 

	//self explanatory
	std::vector<std::vector<double> > bin_edges;
	std::vector<std::vector<double> > bin_widths;

	//Given a string e.g "nu_ICARUS_elike_intrinisc" this map returns the index of the corresponding covariance matrix. Not really used.
	std::map <std::string, std::vector<int> > map_tag_to_covariance_index;
    
    std::map<std::string, std::string> map_subchannel_plotnames;

	// If you have a large covariance matrix/spectrum (say containing nu and nubar mode) but only want to run with nu-mode (so you Set use=0 in nubarmode) the vector used_bins contains all the bins that are actually in use. 
	std::vector<int> used_bins; 

	//For generating a covariance matrix from scratch, this contains the number of montecarlos (weights in weight vector) and their names.
	// For some reason I have decided that the first montecarlo, weight[0] must be the central value, =1
	int num_montecarlo_files;
	//std::vector<int> num_montecarlo;
	std::vector<std::string> montecarlo_name;	 //name means treenae here
	std::vector<std::string> montecarlo_file;	
    std::vector<std::string> montecarlo_additional_weight_names;
    std::vector<std::string> montecarlo_eventweight_branch_names;
    std::vector<bool> montecarlo_additional_weight_bool;
    std::vector<double> montecarlo_additional_weight;
    std::vector<TTreeFormula*> montecarlo_additional_weight_formulas;

    std::vector<std::string> weightmaps_formulas;
    std::vector<std::string> weightmaps_uses;
    std::vector<std::string> weightmaps_patterns;
    std::vector<std::string> weightmaps_mode;

    std::map<std::string,bool> variation_whitelist;
    std::map<std::string,bool> variation_blacklist;

    //A map between a MC file and its friends
    std::map<std::string,std::vector<std::string>> montecarlo_file_friend_map;
    std::map<std::string,std::vector<std::string>> montecarlo_file_friend_treename_map;


	std::vector<int> montecarlo_maxevents;	
	std::vector<double> montecarlo_scale;	
	std::vector<double> montecarlo_pot;	
	std::vector<bool> montecarlo_fake;
	
	std::vector<double> pot_scaling;
	std::vector<double> pot;
	
	std::vector<std::vector<std::string>> parameter_names;	//obsolete code
	std::vector<std::vector<BranchVariable*>> branch_variables;

 	/**********************created for single photon****************************/
        //systematics root files provided correspond to
        std::vector<std::string> systematic_name;
	/*********************************** Member Functions ********************************/	

	int CalcTotalBins();
	
};

}


#endif
