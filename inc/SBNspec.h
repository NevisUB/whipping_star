#ifndef SBNSPEC_H_
#define SBNSPEC_H_

#include <cmath>
#include <map>
#include <vector>
#include <iostream>
#include "SBNconfig.h"
#include <TH1D.h>
#include <string>
#include <TF1.h>
#include <THStack.h>
#include <TMatrixT.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TColor.h>
#include <TStyle.h>
#include <TText.h>

//#include <TROOT.h>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <numeric>

#include <ctime>
#include <TFile.h>
#include "params.h"
#include <TRandom3.h>


template <typename T> 
std::vector<size_t> SortIndexes(const std::vector<T> &v) { 

	// initialize original index locations 
	std::vector<size_t> idx(v.size()); 
	iota(idx.begin(), idx.end(), 0); 

	// sort indexes based on comparing values in v 
	sort(idx.begin(), idx.end(), 
			[&v](size_t i1, size_t i2) {return v[i1] > v[i2];}); 

	return idx; 
}

namespace sbn{
	//This is the basic class that holds all spectral information in whatever reco or true variable you have decided you want in the xml files.
	// Inherits from SBNconfig as thats how its all configured/kept equal! :



	class SBNspec : public SBNconfig{

		public:


			SBNspec() {};
			SBNspec(std::string); //Load in config file EMPTY hists
			SBNspec(std::string, int); //Load in config file, create EMPTY hists, with optional numbering (e.g for multisims!) 
			SBNspec(std::string, int, bool); //Load in config file, create EMPTY hists, with optional numbering (e.g for multisims!)
			SBNspec(std::string, int, bool, bool); //Load in config file, plus the choice of using root files with systematic variations already applied (the second 'bool'). 

			SBNspec(std::string, std::string);
			SBNspec(std::string, std::string, bool);

			SBNspec(std::vector<double> input_full_vec, std::string whichxml);
			SBNspec(std::vector<double> input_full_vec, std::string whichxml, bool isverbose);
            SBNspec(std::vector<double> input_full_vec, std::string whichxml, int universe, bool isverbose);


			// this vector of hists contains all spectra used.
			// The order of filling is the same as the order defined in xml file!
			std::vector<TH1D > hist;
			std::map<std::string, int> map_hist;

			//This is the full concatanated vector (in xml order)	
			std::vector<double > full_vector;
			//This is the compessed vector, collapsing all subchannels down to a single channel
			std::vector<double > collapsed_vector;
			std::vector<float > f_collapsed_vector;

			//error vector for full subchannels, and for the collapsed channel
			//would be instrisinc statistical error
			std::vector<double> full_err_vector;
			std::vector<double> collapsed_err_vector;

			//need to store a history of the scales for oscillation purposes.  FIX THIS
			std::string scale_hist_name;
			double scale_hist_val;
			bool has_been_scaled;



			/*************************** Member Functions ****************/


			int ScaleRandom(); //mainly a debugging function, just randomly scales each hist by 0-2
			int ScalePoisson(); //Scales every histogram by a poissonian random number
			int ScalePoisson(TRandom3*); //Scales every histogram by a poissonian random number


			int SetAsGaussian(double mean, double sigma, int n);
			int SetAsFlat(double val);

			//Scales all vectors in hist by double
			int ScaleAll(double);
			//Scales all vectors whose xml names contains the string name. so you can scale all nu mode at one with Scale("nu");
			int Scale(std::string name, double val);
			//Same as above, but applies a bin-center dependant fucntion to the selected histograms
			int Scale(std::string name, TF1 *);
			//Same as above but normalises to value rather than scales
			int NormAll(double);
			int Norm(std::string name, double val);
			int Keep(std::string name, double val);  // only keep histograms whose name contains string name, and scale it by val
			int Clear();

			//Addes two SBNspec toGether. must have same xml!
			int Add(SBNspec*);

			//Addes histin to the mode_det_channel_subchannel 
			int Add(std::string which_hist, TH1 * histin);


			//Recaculates the full_vector and collapsed_vector's
			int CalcFullVector();
			int CollapseVector();

			//calculate the error vector
			int CalcErrorVector();

			double GetTotalEvents();

			int GetGlobalBinNumber(double invar, int which_hist);
			int GetGlobalBinNumber(int local_bin, std::string histname);
			int GetLocalBinNumber(double invar, int which_hist);

			int GetHistNumber(int f);


			//Just some debugging/checking
			int PrintFullVector();
			int PrintCollapsedVector();
			//WriteOut saves all to an externam rootfile, each individual subchannel and a stacked channel plot.
			int WriteOut(std::string);
			
            std::vector<int> GetIndiciesFromSubchannel(std::string const & subchannel);
			int CompareSBNspecs(SBNspec * compsec, std::string tag);
			int CompareSBNspecs(SBNspec * compsec, bool, std::string tag);
			int CompareSBNspecs(TMatrixT<double> M, SBNspec * compsec, std::string tag);
			int CompareSBNspecs(TMatrixT<double> M, SBNspec * compsec, bool, std::string tag);
			};



}

#endif
