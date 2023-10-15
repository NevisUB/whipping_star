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
        //std::vector<SBNspec*> m_cv_spec_grid;//TMP public
        std::vector<SBNchi*> m_sbnchi_grid;

        TMatrixT<double> * m_full_fractional_covariance_matrix;



        SBNosc *m_core_spectrum;
        SBNosc *m_background_spectrum;
        SBNchi *m_background_chi;
        TVectorT<double> *m_tvec_background_spectrum;
        TVectorT<double> *m_tvec_background_mcerr;
        TVectorT<double> *m_tvec_background_err;

        bool m_bool_core_spectrum_set;
        bool m_bool_background_spectrum_set;
        bool m_bool_stat_only;
        bool m_bool_print_comparasons;

        bool m_use_CNP;

        int m_max_number_iterations;
        double m_chi_min_convergance_tolerance;

        int m_num_universes;
        double m_random_seed;
        std::string tag;

        double global_scale;

        public:
        std::vector<SBNspec*> m_cv_spec_grid;



        SBNfeld(NGrid ingrid, std::string intag,  std::string inxmlname) : SBNconfig(inxmlname), m_grid(ingrid), tag(intag) {
            m_vec_grid = m_grid.GetGrid();
            m_num_total_gridpoints = m_grid.f_num_total_points;
            m_bool_core_spectrum_set = false;
            m_bool_background_spectrum_set = false;
            m_bool_stat_only = false;
            m_num_universes = 2500;
            m_bool_print_comparasons = true;
            m_random_seed = -1;
            m_max_number_iterations = 5;
            m_use_CNP = false;

            global_scale = 1;//5.81731e19/6.6e20;
            //global_scale = 4.51931e+19/6.6e20;

            m_chi_min_convergance_tolerance = 0.00001;
        }


        //Member Functions
        
         double UpdateInverseCovarianceMatrixCNP(size_t best_grid_point, const std::vector<float> &datavec, TMatrixT<double>& inverse_collapsed, SBNchi * helper);
         int UpdateInverseCovarianceMatrix(size_t best_grid_point, TMatrixT<double>& inverse_collapsed, SBNchi * helper);
         std::vector<double> PerformIterativeGridFit(const std::vector<float> &datavec, const size_t grid_pt, const TMatrixT<double>& inverse_background_collapsed_covariance_matrix);
         std::vector<double> PerformIterativeGridFit(const std::vector<float> &datavec, const size_t grid_pt, const TMatrixT<double>& inverse_background_collapsed_covariance_matrix,bool);


        int UseCNP(){m_use_CNP = true;};
        int FullFeldmanCousins();
        int CompareToData(SBNspec *datain);
        int CompareToData(SBNspec *datain,std::vector<double>,std::vector<double>);
        int PointFeldmanCousins(size_t);
        std::vector<double> SBNScan();
        std::vector<double> GlobalScan();
        std::vector<double> GlobalScan(int);
        std::vector<double> GlobalScan2(int);
        std::vector<double> GlobalScan(SBNspec *obs);
        int RasterScan(); 
        
        int GenerateOscillatedSpectra();
        
        int LoadPreOscillatedSpectrum(int);
        int LoadPreOscillatedSpectra();

        int SetRandomSeed(double);

        int GenerateBackgroundScaledSpectrum();
        int GenerateBackgroundSpectrum(); 
        int SetBackgroundSpectrum(std::string filein, std::string scale_nam, double val);
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
       

        std::vector<double> getConfidenceRegion(TGraph *gmin, TGraph *gmax,double val);
        std::vector<double> getConfidenceRegion(double high, double low,double val);
        
        int AddFlatDetSystematic(double percent);
        std::vector<TGraph*> LoadFCMaps(std::string filein);
        std::vector<TGraph*> MakeFCMaps(std::string filein);
        std::vector<TGraph*> MakeFCMaps(std::string filein,size_t);
        

        std::vector<double> MakeMedianMaps(std::string filein,size_t i);



        int GenerateScaledSpectra();
        std::string m_subchannel_to_scale;

        //This is a stopgap for better SBNchi integration.Hrump, need to fix that wierd float oddity. 
        float CalcChi(const std::vector<float>& data, const std::vector<double>& prediction, const TMatrixT<double> & inverse_covariance_matrix );

    };


int rootEigenWrapper(std::string const & xmlname, std::string const & tag, std::string const & signal_subchannel, double bkg_scale_value, NGrid & mygrid, int number , bool bool_stat_only);

}
#endif
