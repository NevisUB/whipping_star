#ifndef SBNCHI_H_
#define SBNCHI_H_

#include <cmath>
#include <vector>
#include <iostream>
#include "SBNspec.h"
#include "SBNconfig.h"

#include "TH2.h"
#include "TMatrixT.h"
#include "TRandom3.h"
#include "TFile.h"

#include <ctime>
#include "params.h"


#include "TMatrixDEigen.h"
#include "TMatrixDSymEigen.h"

namespace sbn{


class SBNchi : public SBNconfig{
	
	public:
	//This is the background spectra that you are comparing too. Obviously can be any existing injected spectra either.
	SBNspec bkgSpec;
	bool stat_only;

	//always contains the last chi^2 value calculated
	double lastChi;
	std::vector<std::vector<double>> lastChi_vec;
	

	//Either initilize from a SBNspec (and use its .xml file)
	SBNchi(SBNspec); 
	//Either initilize from a SBNspec and another xml file
	SBNchi(SBNspec,std::string); 
	//Either initilize from a SBNspec  a TMatrix you have calculated elsewhere
	SBNchi(SBNspec,TMatrixT<double>);
	//If you want to change background, run after change bkgSpec internally 
	int load_bkg();



	//load up systematic covariabnce matrix from a rootfile, location in xml
	TMatrixT<double> sys_fill_direct(std::string, std::string);
	TMatrixT<double> sys_fill_direct();

	TMatrixT<double> matrix_collapsed;

	void fake_fill(TMatrixT <double>&  M);
	void stats_fill(TMatrixT <double>&  M, std::vector<double> diag);


	// These are the powerhouse of of the SBNchi, the ability to collapse any number of modes,detectors,channels and subchannels down to a physically observable subset
	// layer 1 is the cheif bit, taking each detector and collapsing the subchannels
	void collapse_layer1(TMatrixT <double> & M, TMatrixT <double> & Mc);
	//layer 2 just loops layer_1 over all detectors in a given mode
	void collapse_layer2(TMatrixT <double> & M, TMatrixT <double> & Mc);
	//layer 3 just loops layer_2 over all modes we have setup
	void collapse_layer3(TMatrixT <double> & M, TMatrixT <double> & Mc);

	TMatrixT<double> * getCompressedMatrix();

	//Return chi^2 from eith a SBnspec (RECCOMENDED as it checks to make sure xml compatable)
	double CalcChi(SBNspec sigSpec);
	// Or a vector
	double CalcChi(std::vector<double> );
	//Or you are taking covariance from one, and prediciton from another
	double CalcChi(SBNspec *sigSpec, SBNspec *obsSpec);
	//or a log ratio (miniboone esque)
	double CalcChiLog(SBNspec *sigSpec);


	//Some reason eventually store the reuslt in vectors, I think there was memory issues. 
	std::vector<std::vector<double >> to_vector(TMatrixT <double> McI);
	std::vector<std::vector<double >> vMcI;
	std::vector<std::vector<double >> vMc;

	TH2D getChiogram();

	// Todo:
	std::vector<double >  calc_signal_events(struct neutrinoModel &nuModel);


	//int init_minim();
	//double minim_CalcChi(const double * x);
	//double minimize(neutrinoModel newModel, double ipot, double ipotbar);

};


};
#endif
