#ifndef SBNSINGLEPHOTON_H_
#define SBNSINGLEPHOTON_H_

#include <cmath>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include "SBNspec.h"
#include "SBNosc.h"
#include "SBNchi.h"
#include "SBNconfig.h"
#include "SBNgenerate.h"

#include "TH1.h"
#include "TH2.h"
#include "TMatrixT.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TStyle.h"
#include "TLine.h"
#include "TLatex.h"
#include "TText.h"
#include "TList.h"
#include "TMath.h"
#include "TGraph.h"
#include "TSystem.h"
#include "TMarker.h"
#include "TObjArray.h"
#include "TArrow.h"

#include "TMath.h"
#include <ctime>
#include "params.h"

#include "TDecompChol.h"
#include "TDecompSVD.h"
#include "TMatrixDEigen.h"
#include "TMatrixDSymEigen.h"

#include "Math/ProbFunc.h"
#include "Math/DistFunc.h"

#include "prob.h"
#include "ngrid.h"
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>


namespace sbn{
 /***********
 *  A singlephoton class built purely for perform NCpi0 normalization fit, NC delta fit for gLEE analysis
 *  	Guanqun Ge, July 19 2020
 */

    class SBNsinglephoton: public SBNconfig{


	//regular grid, for a flat normalization
        NGrid m_grid;
        int m_flat_total_gridpoints;
        std::vector<std::vector<double>> m_vec_grid;

        // extra grid for polynomial scaling purpose
	NGrid m_poly_grid; 
	int m_poly_total_gridpoints;
	bool m_bool_poly_grid;  //use polynomial grid or not
	std::vector<std::vector<double>> m_vec_poly_grid;

	int m_fit_dimension;   //total num of parameters being fit
	int m_total_gridpoints;  //total number of grid points 

	int m_cv_spec_index = -1;   //global index of cv spectrum
        std::vector<SBNspec> m_scaled_spec_grid;  

        TMatrixT<double> * m_full_fractional_covariance_matrix; //full systematic( include detector) covariance matrix
        TMatrixT<double> * m_full_but_genie_fractional_covariance_matrix;   //full systematic covariance matrix other than genie
        TMatrixT<double> * m_genie_fractional_covariance_matrix;   //genie only covariance matrix

	SBNspec* m_cv_spectrum;   //genie CV spectra
	std::string asimov_dataset;
	SBNspec* m_data_spectrum; // data spectra
        SBNchi*  m_chi;

	double m_cv_delta_scaling;   //NCdelta scaling,should only be applied to FAKE DATA for sensitivity study!
	bool m_bool_modify_cv;  //if CV is modified for systematic covar matrix calculation
	bool m_bool_cv_spectrum_generated;   //if CV spec is generated
	bool m_bool_cv_spectrum_loaded;   //if CV spec is loaded, for fit purposes
	bool m_bool_data_spectrum_loaded;  //if data spec is loaded

        bool m_use_CNP;

        int m_max_number_iterations;
        double m_chi_min_convergance_tolerance;

        double m_random_seed;
        std::string tag;   //tag for the analysis: NCpi0, NCDelta etc

	TFile* fin;  //file that has vector of chi, and BF index in it.
	int m_interpolation_number;
	std::map<int, std::vector<double>> m_map;   // map of best-fit point index and chi2 vector

	//things used only when generating pre-scaled files
	int num_files;
        std::vector<int> nentries;
        std::vector<TFile *> files;
        std::vector<TTree *> trees;

	void ZeroOutCorrelation(TMatrixT<double>*, const std::string& );
	int RemoveNan(TMatrixT<double>*); //remove the nan's from matrix

        public:

	SBNsinglephoton(std::string xmlname, std::string intag);
	SBNsinglephoton(std::string xmlname, std::string intag, NGrid ingrid);
	SBNsinglephoton(std::string xmlname, std::string intag, NGrid ingrid, NGrid in_polygrid, bool has_polygrid);

	// MEMBER FUNCTION//
	
	
        int PreScaleSpectrum(std::string xmlname, double, std::vector<double>& param);	
        int PreScaleSpectrum(std::string xmlname, std::vector<double>& param);  //generate scaled spectrum per set of parameter

	/* basically the same function as PreScaleSpectrum, only difference is it returns the scaled spec */
        int ScaleSpectrum(SBNspec*, double, std::vector<double>& param);	
	int GeneratePreScaledSpectra();    //generate pre-scaled spectra for full polynomial grid
	int LoadSpectraApplyFullScaling();  //load spectra and apply full scaling
	int LoadSpectraOnTheFly();
	//calc scale factors based on event energy and polynomial parameters
        double ScaleFactor(double E, double, std::vector<double>& param);
	//int WriteOutCV(std::string tag); 


	int CalcChiGridScanShapeOnlyFit();  //calculate chi2 surface for NCpi0 fit.
	int CalcChiGridScan();  //simply grid scan with covariance matrix fixed during each iteration
	int CalcChiGridScanVaryMatrices();  //simply grid scan with covariance matrix varied at each grid point
	double CalcChi(bool use_cnp);

	int SetStatOnly();
	int SetFlatFullFracCovarianceMatrix(double );
	int SetFullFractionalCovarianceMatrix(std::string filename, std::string matrix_name);
	int SetGenieFractionalCovarianceMatrix(std::string filename);
	int AddCovarianceMatrix(const std::string &filename, const std::string &covarname);
	int CalcFullButGenieFractionalCovarMatrix();
	void PrintCollapsedFractionalMatrix(const TMatrixT<double> &, const SBNspec &);
	void ZeroOutGenieCorrelation(const std::string& );
	void ZeroOutFullCorrelation(const std::string& );
	void ZeroOutOffDiagonal();

	int LoadCV();
	void SetupAsimovDataset(const std::string &);
	int LoadData(std::string filename);
	int SetPolyGrid(NGrid ingrid);
	
	void PublicDataPrintOut(const std::vector<std::string>& subchannel_names, const std::vector<double>& scale);
	int PrintOutFitInfo(const std::map<int, std::vector<double>>& , std::string tag, bool);
	int WriteOutInfo(std::map<int, std::vector<double>>& );
	//not finished yet, what information do we wnat to print out
	int GrabFitMap();
	void SetInterpolationNumber(int );
	int SaveHistogram(bool sensitivity=false);
	int SaveHistogram(std::map<int, std::vector<double>>& , bool sensitivity=false);
	TH2D Do2DInterpolation(int, const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& value, std::string);
	std::vector<TGraph> FindContour(TH2D&, int n, std::string);
	void DrawContour(TH2D*, std::vector<TGraph>&, std::vector<double>);
	void DrawContour(TH2D*, std::vector<TGraph>&, std::string, std::vector<double>);
	
	int ModifyCV(double factor);
	int ModifyCV(double, std::vector<double> param);

	SBNspec GeneratePointSpectra(int );

	protected:

	int OpenFiles();
	bool m_file_open;
	int CloseFiles();	
	int PoissonFluctuation(SBNspec *);

	private:
	void MaskScaleFactor(double&);
	void DrawMarginalizedChi(std::vector<TH1D*> vec_hist, std::string intag, bool sensitivity=false, std::vector<double> bf={});
	void LocateCVGlobalIndex();
	void CheckCVLoad();
	void CheckNumScaledSpectra();
	void CheckDataLoad();

    };



}
#endif
