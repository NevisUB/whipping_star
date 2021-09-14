#ifndef SBNCHI_H_
#define SBNCHI_H_

#include <cmath>
#include <vector>
#include <numeric>
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
#include "TFile.h"
#include "TStyle.h"
#include "TLine.h"
#include "TMath.h"
#include "TLatex.h"
#include "params.h"

#include "TDecompChol.h"
#include "TDecompSVD.h"
#include "TMatrixDEigen.h"
#include "TMatrixDSymEigen.h"
#include <Eigen/Dense>
#include <Eigen/SVD>

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
	bool is_stat_only; //controls whether systematic uncertainty (including MC intrinsic error) will be included

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
    std::normal_distribution<float>* m_dist_normal;

	/*********************************** Member Functions ********************************/	



	double CalcChi_statonlyCNP(std::vector<double>&, std::vector<double>&);
	int ReloadCoreSpectrum(SBNspec *bkgin);

	//load up systematic covariabnce matrix from a rootfile, location in xml
	//These are pretty obsolete.
	TMatrixT<double> FillSystematicsFromXML(std::string, std::string);
	TMatrixT<double> FillSystematicsFromXML();

	void FakeFillMatrix(TMatrixT <double>&  M);
	void FillStatsMatrix(TMatrixT <double>&  M, std::vector<double> diag);
	TMatrixT<double> FillSystMatrix(const TMatrixT<double>& M, const std::vector<double>& spec, const std::vector<double>& spec_err);
	TMatrixT<double> FillSystMatrix(const TMatrixT<double>& M, const std::vector<double>& spec, const std::vector<double>& spec_err, bool);

	// These are the powerhouse of of the SBNchi, the ability to collapse any number of modes,detectors,channels and subchannels down to a physically observable subSet
	// layer 1 is the cheif bit, taking each detector and collapsing the subchannels
	void CollapseSubchannels(TMatrixT <double> & M, TMatrixT <double> & Mc);
	//layer 2 just loops layer_1 over all detectors in a given mode
	void CollapseDetectors(TMatrixT <double> & M, TMatrixT <double> & Mc);
	//layer 3 just loops layer_2 over all modes we have Setup
	void CollapseModes(TMatrixT <double> & M, TMatrixT <double> & Mc);

    TMatrixT<double> InvertMatrix(TMatrixT<double> &M);
    TMatrixT<double> CalcCovarianceMatrix(TMatrixT<double>*M, TVectorT<double>& spec, TVectorT<double>& spec_err);
    TMatrixT<double> CalcCovarianceMatrix(TMatrixT<double>*M, std::vector<double>& spec, std::vector<double>& spec_err);
    TMatrixT<double> SplitCovarianceMatrix(TMatrixT<double>*M, std::vector<double>& spec, std::vector<double>& spec_err, int m);
    TMatrixT<double> CalcCovarianceMatrixCNP(TMatrixT<double> &M, std::vector<double>& spec, std::vector<double>& spec_err,  std::vector<double>& spec_collapse, const std::vector<double>& datavec );
    TMatrixT<double> CalcCovarianceMatrixCNP(TMatrixT<double>* M, std::vector<double>& spec, std::vector<double>& spec_err, std::vector<double>& spec_collapse, const std::vector<float>& datavec );
    TMatrixT<double> CalcCovarianceMatrixCNP(TMatrixT<double>* M, std::vector<double>& spec, const std::vector<float>& datavec );
    TMatrixT<double> CalcShapeOnlyCovarianceMatrix(TMatrixT<double> &M, SBNspec *mc, SBNspec* bkg);    
    TMatrixT<double> CalcShapeMixedCovarianceMatrix(TMatrixT<double> &M, SBNspec *mc, SBNspec* bkg);    
    TMatrixT<double> CalcNeymanCovarianceMatrix(TMatrixT<double>* M, std::vector<double>& spec, std::vector<double>& spec_err, std::vector<double>& data_full_vec);
    TMatrixT<double> AddStatMatrixCNP(TMatrixT<double>* M, std::vector<double>& spec, const std::vector<double>& datavec );
    TMatrixT<double> AddStatMatrix(TMatrixT<double>* M, const std::vector<double>& datavec );

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
	double CalcChi(TMatrixT<double> M, std::vector<double>& spec, std::vector<double>& data);
	double CalcChi(TMatrixT<double> M, std::vector<double>& spec, std::vector<double>& data, bool);
	double CalcShapeChi(TMatrixT<double> M, std::vector<double>& data, std::vector<double>& constrained_bkgd, std::vector<double>& mc, bool is_background);
	std::vector<std::vector<double >> TMatrixDToVector(TMatrixT <double> McI);
	

	//Cholosky related
	int PerformCholoskyDecomposition(SBNspec *specin);

//	SBNspec SampleCovariance(SBNspec *specin); 
    std::vector<float> SampleCovariance(SBNspec *specin); 

    //

	TH1D SamplePoissonVaryCore(SBNspec *specin, int num_MC);
	TH1D SamplePoissonVaryInput(SBNspec *specin, int num_MC, double maxchi);
	TH1D SamplePoissonVaryInput(SBNspec *specin, int num_MC, std::vector<double>*);
	TH1D SampleCovarianceVaryInput(SBNspec *specin, int num_MC, double maxchi);
	TH1D SampleCovarianceVaryInput(SBNspec *specin, int num_MC, std::vector<double>*);


    TH1D SamplePoisson_NP(SBNspec *specin, SBNchi &chi_h0, SBNchi & chi_h1, int num_MC, std::vector<double> *chival,int which_sample);
    TH1D SamplePoisson_NP(SBNspec *specin, SBNchi &chi_h0, SBNchi & chi_h1, int num_MC, double,int which_sample);



    double max_sample_chi_val;

	int CollapseVectorStandAlone(std::vector<double> * full_vector, std::vector<double> *collapsed_vector);

	int CollapseVectorStandAlone(double* full_vector, double* collapsed_vector);
	int CollapseVectorStandAlone(float* full_vector, float* collapsed_vector);



    int SingleValueDecomposition(double ** matrix, double ** U, double**V, double *single_values );

    std::vector<float> GeneratePseudoExperiment();


		//some plotting things
	TH2D* GetChiogram();
	int PrintMatricies(std::string);
        int DrawSampleCovariance(std::string);
	int DrawComparisonIndividualFracMatrix(SBNspec&, SBNspec&, TMatrixT<double>&, std::string);
	int DrawComparisonIndividualFracMatrix(SBNspec&, SBNspec&, TMatrixT<double>&, std::string, bool);
	int DrawComparisonIndividual(SBNspec& sig, SBNspec& data, std::string);
	int DrawComparisonIndividual(SBNspec& sig, SBNspec& data, std::string, bool);
	int DrawComparisonIndividual(SBNspec& sig, SBNspec& data, TMatrixT<double>&, std::string);
	int DrawComparisonIndividual(SBNspec& sig, SBNspec& data, TMatrixT<double>&, std::string, bool);
};


};
#endif
