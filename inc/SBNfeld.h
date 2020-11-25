#ifndef SBNFELD_H_
#define SBNFELD_H_

#include <cmath>
#include <vector>
#include <iostream>
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
#include "TMath.h"
#include "TGraph.h"

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


namespace sbn{


    class SBNfeld: public SBNconfig{


        NGrid m_grid;
        int m_num_total_gridpoints;

        std::vector<std::vector<double>> m_vec_grid;
        std::vector<SBNspec*> m_cv_spec_grid;
        std::vector<SBNchi*> m_sbnchi_grid;

        TMatrixT<double> * m_full_fractional_covariance_matrix;

        SBNosc *m_core_spectrum;
        SBNosc *m_background_spectrum;
        SBNchi *m_background_chi;
        TVectorT<double> *m_tvec_background_spectrum;
        TVectorT<double> *m_tvec_background_err;

        bool m_bool_core_spectrum_set;
        bool m_bool_background_spectrum_set;
        bool m_bool_stat_only;
        bool m_bool_print_comparasons;

        int m_max_number_iterations;
        double m_chi_min_convergance_tolerance;

        int m_num_universes;
        double m_random_seed;
        std::string tag;

        double global_scale;

        public:

        SBNfeld(NGrid ingrid, std::string intag,  std::string inxmlname) : SBNconfig(inxmlname), m_grid(ingrid), tag(intag) {
            m_vec_grid = m_grid.GetGrid();
            m_num_total_gridpoints = m_grid.f_num_total_points;
            m_bool_core_spectrum_set = false;
            m_bool_background_spectrum_set = false;
            m_bool_stat_only = false;
            m_num_universes = 2500;
            m_bool_print_comparasons = true;
            m_random_seed = -1;
            //m_max_number_iterations = 5;
            m_max_number_iterations = 10;


            global_scale = 1.0;
            //global_scale = 5.23537/66.0;

            m_chi_min_convergance_tolerance = 0.001;
        }


        //Member Functions
        
	//use Pearson chi during the fit
        std::vector<double> PerformIterativeFit(std::vector<float> &datavec, size_t grid_pt, TMatrixT<double>& inverse_background_collapsed_covariance_matrix);
 	// use Neyman chi during the fit
        std::vector<double> PerformIterativeFitNeyman(std::vector<float> &datavec, size_t grid_pt, TMatrixT<double>& inverse_background_collapsed_covariance_matrix);


        
        int FullFeldmanCousins();
        int PointFeldmanCousins(size_t);
        int GlobalScan();
        int GlobalScan(int);
        int GlobalScan(SBNspec *obs);
        int GlobalScanNeyman(SBNspec *obs);
        int RasterScan(); 
        
        int GenerateOscillatedSpectra();
        
        int LoadPreOscillatedSpectrum(int);
        int LoadPreOscillatedSpectra();

        int SetRandomSeed(double);

        int GenerateBackgroundSpectrum(); 
        int LoadBackgroundSpectrum();
        int LoadBackgroundSpectrum(std::string);

        int CalcSBNchis();

        int SetCoreSpectrum(std::string);
        int SetFractionalCovarianceMatrix(TMatrixT<double> *);
        int SetFractionalCovarianceMatrix(std::string, std::string);
        int SetEmptyFractionalCovarianceMatrix();
        int SetNumUniverses(int);
        int SetStatOnly();

        NeutrinoModel convert3p1(std::vector<double> ingrid);
        
        
        int AddFlatDetSystematic(double percent);
        


        int GenerateScaledSpectra();
        std::string m_subchannel_to_scale;

        //This is a stopgap for better SBNchi integration.Hrump, need to fix that wierd float oddity. 
        float CalcChi(std::vector<float>& data, std::vector<double>& prediction, TMatrixT<double> & inverse_covariance_matrix );
	//whther or not include determinant of covariance matrix in chi calculation
        float CalcChi(std::vector<float>& data, std::vector<double>& prediction, TMatrixT<double> & inverse_covariance_matrix, bool det_term);

    };




}
#endif
