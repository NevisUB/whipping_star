#ifndef SBNCHI_H_
#define SBNCHI_H_

#include <cmath>
#include <vector>
#include <iostream>
#include <ctime>
#include <random>

#include "SBNspec.h"
#include "SBNconfig.h"

#include "TH1.h"
#include "TH2.h"
#include "TMatrixT.h"
#include "TMatrixTSym.h"
#include "TVectorT.h"
#include "TRandom3.h"
#include "TText.h"
#include "TFile.h"
#include "TStyle.h"
#include "TLine.h"

#include "params.h"

#include "TDecompChol.h"
#include "TDecompSVD.h"
#include "TMatrixDEigen.h"
#include "TMatrixDSymEigen.h"

namespace sbn{


class SBNchi : public SBNconfig{

	public:

	//Either initilize from a SBNspec (and use its .xml file)
	SBNchi(SBNspec);
	//Either initilize from a SBNspec and another xml file
	SBNchi(SBNspec,std::string);
	//Either initilize from a SBNspec  a TMatrix you have calculated elsewhere
	SBNchi(SBNspec,TMatrixT<double>);
	SBNchi(SBNspec,TMatrixT<double>,bool);
	SBNchi(SBNspec,TMatrixT<double>,std::string, bool);
	SBNchi(SBNspec in, TMatrixT<double> matrix_systematicsin, std::string inxml, bool is_verbose, double random_seed);

    //Initialise a stat_only one;
	SBNchi(SBNspec, bool is_stat_only);
	SBNchi(std::string);
	

	//This is the core spectra that you are comparing too. This is used to calculate covariance matrix and in a way is on the 'bottom' of the chi^2.
	SBNspec core_spectrum;
	bool is_stat_only;

	//always contains the last chi^2 value calculated
	double last_calculated_chi;
	std::vector<std::vector<double>> vec_last_calculated_chi;



	TMatrixT<double> matrix_systematics;
	TMatrixT<double> matrix_fractional_covariance;
	TMatrixT<double> matrix_collapsed;

	//Used in cholosky decompositions
	bool cholosky_performed;
	TMatrixT<float> matrix_lower_triangular;
	std::vector<std::vector<float>> vec_matrix_lower_triangular;

	//Some reason eventually store the reuslt in vectors, I think there was memory issues.
	std::vector<std::vector<double >> vec_matrix_inverted;
	std::vector<std::vector<double >> vec_matrix_collapsed;

    /***** Random Number Generation ****/
    std::random_device random_device_seed;
    std::mt19937 *rangen_twister; //merseinne twister
    std::minstd_rand * rangen_linear;
    std::ranlux24_base * rangen_carry;
    void InitRandomNumberSeeds();
    void InitRandomNumberSeeds(double);
    TRandom3 * rangen;

	/*********************************** Member Functions ********************************/	


	int ReloadCoreSpectrum(SBNspec *bkgin);

	//load up systematic covariabnce matrix from a rootfile, location in xml
	//These are pretty obsolete.
	TMatrixT<double> FillSystematicsFromXML(std::string, std::string);
	TMatrixT<double> FillSystematicsFromXML();
	TMatrixT<double> FillSystMatrix(TMatrixT<double>* M, std::vector<double>& spec, std::vector<double>& err);

	void FakeFillMatrix(TMatrixT <double>&  M);
	void FillStatsMatrix(TMatrixT <double>&  M, std::vector<double> diag);

	// These are the powerhouse of of the SBNchi, the ability to collapse any number of modes,detectors,channels and subchannels down to a physically observable subSet
	// layer 1 is the cheif bit, taking each detector and collapsing the subchannels
	void CollapseSubchannels(TMatrixT <double> & M, TMatrixT <double> & Mc);
	//layer 2 just loops layer_1 over all detectors in a given mode
	void CollapseDetectors(TMatrixT <double> & M, TMatrixT <double> & Mc);
	//layer 3 just loops layer_2 over all modes we have Setup
	void CollapseModes(TMatrixT <double> & M, TMatrixT <double> & Mc);

    TMatrixT<double> InvertMatrix(TMatrixT<double> &M);
    TMatrixT<double> CalcCovarianceMatrix(TMatrixT<double>*M, TVectorT<double>& spec, TVectorT<double>& err);
    TMatrixT<double> CalcCovarianceMatrix(TMatrixT<double>*M, std::vector<double>& spec, std::vector<double>& err);
    TMatrixT<double> CalcCovarianceMatrix(TMatrixT<double>*M, std::vector<double>& spec, std::vector<double>& err, std::vector<double>& datavec);
    TMatrixT<double> CalcCovarianceMatrix(TMatrixT<double>*M, std::vector<double>& spec, std::vector<double>& err, std::vector<float>& datavec);
    TMatrixT<double> SplitCovarianceMatrix(TMatrixT<double>*M, std::vector<double>& spec, std::vector<double>& err, int m);
    TMatrixT<double> ChangeToNeymanStats(TMatrixT<double>*M, std::vector<float>& mcspec, std::vector<float>& datavec);   

	TMatrixT<double> * GetCollapsedMatrix();
	int FillCollapsedCovarianceMatrix(TMatrixT<double>*);
	int FillCollapsedCorrelationMatrix(TMatrixT<double>*);
	int FillCollapsedFractionalMatrix(TMatrixT<double>*);

	//Return chi^2 from eith a SBnspec (RECCOMENDED as it checks to make sure xml compatable)
	//double CalcChi(SBNspec sigSpec);
	double CalcChi(SBNspec *sigSpec);
	// Or a vector
	double CalcChi(std::vector<double> );
	//Or you are taking covariance from one, and prediciton from another
	double CalcChi(SBNspec *sigSpec, SBNspec *obsSpec);
	//or a log ratio (miniboone esque)
	double CalcChiLog(SBNspec *sigSpec);

	double CalcChi(std::vector<double> * sigVec);
	float  CalcChi(std::vector<float> * sigVec);
	double CalcChi(double* sigVec);

	double CalcChi(double ** inv, double *, double *);
	float CalcChi(float ** inv, float *, float *);

	std::vector<std::vector<double >> TMatrixDToVector(TMatrixT <double> McI);
	

	//Cholosky related
	int PerformCholoskyDecomposition(SBNspec *specin);

//	SBNspec SampleCovariance(SBNspec *specin); 
    std::vector<float> SampleCovariance(SBNspec *specin); 

	TH1D SamplePoissonVaryCore(SBNspec *specin, int num_MC);
	TH1D SamplePoissonVaryInput(SBNspec *specin, int num_MC, double maxchi);
	TH1D SamplePoissonVaryInput(SBNspec *specin, int num_MC, std::vector<double>*);
	TH1D SampleCovarianceVaryInput(SBNspec *specin, int num_MC, double maxchi);
	TH1D SampleCovarianceVaryInput(SBNspec *specin, int num_MC, std::vector<double>*);

    double max_sample_chi_val;

	int CollapseVectorStandAlone(std::vector<double> * full_vector, std::vector<double> *collapsed_vector);

	int CollapseVectorStandAlone(double* full_vector, double* collapsed_vector);
	int CollapseVectorStandAlone(float* full_vector, float* collapsed_vector);



    int SingleValueDecomposition(double ** matrix, double ** U, double**V, double *single_values );






		//some plotting things
	TH2D* GetChiogram();
	int PrintMatricies(std::string);
	int PlotMatrix(TMatrixT<double> matrix, std::string tag, bool plot_pdf);
    	int DrawSampleCovariance(std::string);

};


};
#endif
