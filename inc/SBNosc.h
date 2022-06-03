#ifndef SBNOSC_H_
#define SBNOSC_H_

#include "SBNspec.h"
#include "prob.h"
#include "params.h"
#include <string>
#include <utility>

/**********************************************
 *	This is a less general 3+N osc spectrum
 * *******************************************/
namespace sbn{
// Note for this to work, you need to have precomputed spectra in data/precomp
// Ability to precompute is not yet included!

// They are labeled SBN_FREQ_MASS_.root where  FREQ is either SIN or SINSQ and MASS is the log10 of the sterile ev^2, e.g -0.04, or 1.20 


class SBNosc : public SBNspec{
	protected:
	public:
	//this is the structure contains 3+N oscillation parameters (find in prob.h)
	struct NeutrinoModel working_model;	
	
	// which_mode to oscillate in  (APP, DIS, etc..) 
	int which_mode;
	double mass_step_size;	//has to be 0.04 for now

	SBNosc(std::string, std::string); //constructor
	SBNosc(std::string, std::string, NeutrinoModel); //constructor
        SBNosc(SBNspec & specin);

	//find all the frequencies! Needs to know if a frequency corresponds to 41 51 54..etc.. so thats the int
	std::vector< std::pair <double, int>> mass_splittings;	

	//Oscillate the contained std::vector<TH1D> hists 
	int OscillateThis(std::string);	
	// Or just oscillate a copy and return the ompressed vector
        std::vector<std::vector<double>> Oscillate(std::string,bool compress); 
        std::vector<std::vector<double>> Oscillate(std::string);
	std::vector<std::vector<double>> Oscillate(std::string, double);
	//std::vector<double> OscillateWithAmp(double amp, double amp_sq);

	int LoadModel(NeutrinoModel);	
	int calcMassSplittings();	

	int PrecomputeSpectra(double dm);


	//Oscillation mode 
	int SetMode(int);
	void SetAppMode();
	void SetDisMode();
	void SetBothMode();
	void SetWierdMode();
	void SetDisEMode();


};

};
#endif
