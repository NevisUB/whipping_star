#ifndef SBNCHI_H_
#define SBNCHI_H_

#include <cmath>
#include <vector>
#include <numeric>
#include <iostream>
#include <ctime>
#include <random>
#include <sys/stat.h>

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
#include "TText.h"
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

    std::vector<TMatrixT<double>> splitNormShape(TMatrixT<double> & Min);

    struct CLSresult{

        public:
            std::string m_tag;
            TH1D m_pdf;
            float m_min_value;
            float m_max_value;
            std::vector<float> m_values;
            std::vector<double> m_quantiles; 
            std::vector<double> m_nlower; 
            std::vector<float> m_pvalues; 

            CLSresult(){
                m_tag = "Default";
                m_min_value = 9999999;
                m_max_value = -999999;
            }
    };


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
	    ~SBNchi(){
	std::cout << "start free SBNchi " << std::endl;
		 if(rangen_twister)   {delete rangen_twister; rangen_twister = nullptr;}
        	 if(rangen_linear)    {delete rangen_linear; rangen_linear = nullptr;}
        	 if(rangen_carry)     {delete rangen_carry; rangen_carry = nullptr;}
        	 if(rangen)	      {delete rangen; rangen = nullptr;}
  		 if(m_dist_normal)  {delete  m_dist_normal; m_dist_normal = nullptr;}
	std::cout << "end free SBNchi " << std::endl;
	    }

	//This is the core spectra that you are comparing too. This is used to calculate covariance matrix and in a way is on the 'bottom' of the chi^2.
	SBNspec core_spectrum;
	bool is_stat_only; //controls whether systematic uncertainty (including MC intrinsic error) will be included



            //always contains the last chi^2 value calculated
            double last_calculated_chi;
            std::vector<std::vector<double>> vec_last_calculated_chi;



            TMatrixT<double> matrix_systematics;
            TMatrixT<double> m_matrix_systematics_collapsed;
            TMatrixT<double> matrix_fractional_covariance;
            TMatrixT<double> matrix_collapsed;

            //Used in cholosky decompositions
            double m_tolerance;
            bool cholosky_performed;
            TMatrixT<float> matrix_lower_triangular;
            std::vector<std::vector<float>> vec_matrix_lower_triangular;



            //Some reason eventually store the reuslt in vectors, I think there was memory issues.
            std::vector<std::vector<double >> vec_matrix_inverted;
            std::vector<std::vector<double >> vec_matrix_collapsed;

            /***** Random Number Generation ****/
            std::random_device random_device_seed;
            std::mt19937 *rangen_twister = nullptr; //merseinne twister
            std::minstd_rand * rangen_linear = nullptr;
            std::ranlux24_base * rangen_carry = nullptr;
            void InitRandomNumberSeeds();
            void InitRandomNumberSeeds(double);
            TRandom3 * rangen = nullptr;
            std::normal_distribution<float>* m_dist_normal = nullptr;

            /*********************************** Member Functions ********************************/	

            double m_cmin;
            double m_cmax;
            int plot_one(TMatrixD matrix, std::string tag, TFile *fin,bool,bool,bool);



            int ReloadCoreSpectrum(SBNspec *bkgin);

            //load up systematic covariabnce matrix from a rootfile, location in xml
            //These are pretty obsolete.
            TMatrixT<double> FillSystematicsFromXML(std::string, std::string);
            TMatrixT<double> FillSystematicsFromXML();


            // These are the powerhouse of of the SBNchi, the ability to collapse any number of modes,detectors,channels and subchannels down to a physically observable subSet
            // layer 1 is the cheif bit, taking each detector and collapsing the subchannels
            void CollapseSubchannels(TMatrixT <double> & M, TMatrixT <double> & Mc);
            //layer 2 just loops layer_1 over all detectors in a given mode
            void CollapseDetectors(TMatrixT <double> & M, TMatrixT <double> & Mc);
            //layer 3 just loops layer_2 over all modes we have Setup
            void CollapseModes(TMatrixT <double> & M, TMatrixT <double> & Mc);

            /*
            TMatrixT<double> InvertMatrix(TMatrixT<double> &M);
            TMatrixT<double> CalcCovarianceMatrix(TMatrixT<double>*M, TVectorT<double>& spec);
            TMatrixT<double> CalcCovarianceMatrix(TMatrixT<double>*M, TVectorT<double>& spec, TVectorT<double> &err);
            TMatrixT<double> CalcCovarianceMatrix(TMatrixT<double>*M, TVectorT<double>& spec, bool);
            TMatrixT<double> CalcCovarianceMatrix(TMatrixT<double>*M, std::vector<double>& spec);
            TMatrixT<double> CalcCovarianceMatrix(TMatrixT<double>*M, std::vector<double>& spec, std::vector<double> &mcerr);
            TMatrixT<double> CalcCovarianceMatrix(TMatrixT<double>*M, std::vector<double>& spec,bool);
            TMatrixT<double> CalcCovarianceMatrixCNP(TMatrixT<double> M, std::vector<double>& spec, std::vector<double>& spec_collapse, const std::vector<double>& datavec );
            TMatrixT<double> CalcCovarianceMatrixCNP(TMatrixT<double>* M, std::vector<double>& spec, const std::vector<float>& datavec );
            TMatrixT<double> CalcCovarianceMatrixCNP(TMatrixT<double>* M, std::vector<double>& spec, std::vector<double>& spec_collapse, std::vector<double>& spec_mcerr, const std::vector<float>& datavec );
            */

    TMatrixT<double> CalcCovarianceMatrix(TMatrixT<double>*M, TVectorT<double>& spec,bool add_stats);
    TMatrixT<double> CalcCovarianceMatrix(TMatrixT<double>*M, std::vector<double>& spec, bool add_stats);
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
            int FillCovarianceMatrix(TMatrixT<double>*);
            int FillCorrelationMatrix(TMatrixT<double>*);
            int FillFractionalMatrix(TMatrixT<double>*);

        void FakeFillMatrix(TMatrixT <double>&  M);
        void FillStatsMatrix(TMatrixT <double>&  M, std::vector<double> diag);
        TMatrixT<double> FillSystMatrix(const TMatrixT<double>& M, const std::vector<double>& spec, const std::vector<double>& spec_err);
        TMatrixT<double> FillSystMatrix(const TMatrixT<double>& M, const std::vector<double>& spec, const std::vector<double>& spec_err, bool);


	        double CalcChi_statonlyCNP(std::vector<double>&, std::vector<double>&);
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

            float PoissonLogLiklihood(float * h0_corein, float *collapsed);
            float CalcChi_CNP(float * pred, float* data);
            float CalcChi_Pearson(float * pred, float* data);
            double CalcChi(TMatrixT<double> M, std::vector<double>& spec, std::vector<double>& data);
	        double CalcChi(TMatrixT<double> M, std::vector<double>& spec, std::vector<double>& data, bool);
        	double CalcShapeChi(TMatrixT<double> M, std::vector<double>& data, std::vector<double>& constrained_bkgd, std::vector<double>& mc, bool is_background);
        
    
            std::vector<std::vector<double >> TMatrixDToVector(TMatrixT <double> McI);

            double setTolerance(double ep){m_tolerance = ep;};

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



            std::vector<CLSresult> Mike_NP(SBNspec *specin, SBNchi &chi_h0, SBNchi & chi_h1, int num_MC, int which_sample,int id);
            TH1D SamplePoisson_NP(SBNspec *specin, SBNchi &chi_h0, SBNchi & chi_h1, int num_MC, std::vector<double> *chival,int which_sample);
            TH1D SamplePoisson_NP(SBNspec *specin, SBNchi &chi_h0, SBNchi & chi_h1, int num_MC, double,int which_sample);

            int SetFracPlotBounds(double cmin,double cmax){
                     m_cmin=cmin;
                     m_cmax=cmax;
                     return 0;
            }

            double max_sample_chi_val;

            int CollapseVectorStandAlone(std::vector<double> * full_vector, std::vector<double> *collapsed_vector);

            int CollapseVectorStandAlone(double* full_vector, double* collapsed_vector);
            int CollapseVectorStandAlone(float* full_vector, float* collapsed_vector);



            int SingleValueDecomposition(double ** matrix, double ** U, double**V, double *single_values );

            bool pseudo_from_collapsed;
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


    struct SBNconstraint{

        public:
        TMatrixD con_inv_B;
        TMatrixD con_inv_C;

        TMatrixD con_B;
        TMatrixD con_C;

        TMatrixD con_input_inv_Cov;

        int n_bins;
        int n_constrain_below;

        float *con_pred;
        float *con_data;
        
        std::vector<float> con_nfit;

        SBNconstraint(TMatrixD &inv_Cov, float *pred, float *data, int bins, int constrain_below){
            con_pred = pred;
            con_data = data;
            n_bins =  bins;
            n_constrain_below = constrain_below;
            con_input_inv_Cov.ResizeTo(n_bins,n_bins);

            con_inv_C.ResizeTo(n_bins,n_bins); con_inv_C.Zero();
            con_inv_B.ResizeTo(n_bins,n_bins); con_inv_B.Zero();
            con_C.ResizeTo(n_bins,n_bins); con_C.Zero();
            con_B.ResizeTo(n_bins,n_bins); con_B.Zero();

            TVectorD v_pred(n_bins);

            for(int i=0; i< n_bins; i++){
                v_pred[i] = pred[i];

                for(int j=0; j< n_bins; j++){ 
                    con_input_inv_Cov(i,j)=inv_Cov(i,j);
                    if(i <= n_constrain_below || j <=n_constrain_below){
                        con_inv_B(i,j) = inv_Cov(i,j);    
                        con_inv_C(i,j) = inv_Cov(i,j);    
                    }else{
                        con_inv_B(i,j) = inv_Cov(i,j);
                        con_inv_C(i,j) = inv_Cov(i,j);
                        if(i==j){
                            con_inv_B(i,j)+=1.0/data[i];
                            con_inv_C(i,j)+=1.0/pred[i];
                        }
                    }

                }
            }//Filled inverse matricies.

            //Invert B
            TDecompSVD svd(con_inv_B);

            if (!svd.Decompose()){
                std::cout<<"SBNconstraint SVD Decomposition failed, matrix not symettric?, has nans?" << std::endl;
                std::cout<<"ERROR: The matrix to invert failed a SVD decomp!"<<std::endl;

                for(int i=0; i< con_inv_B.GetNrows(); i++){
                    for(int j=0; j< con_inv_B.GetNrows(); j++){
                        std::cout<<i<<" "<<j<<" "<<con_inv_B(i,j)<<std::endl;
                    }
                }

                exit(EXIT_FAILURE);

        } else {
            con_B = svd.Invert();
        }
        if( !con_B.IsValid()){
            std::cout<<"SBNconstraint ERROR: The inverted matrix isnt valid! Something went wrong.."<<std::endl;
            exit(EXIT_FAILURE);
        }

        // Multiply out to find con_nfit;
    
        TVectorD v_fit(n_bins);
        con_nfit.resize(n_bins);
        v_fit = con_B*con_inv_C*v_pred;
        for(int i=0; i< n_bins;i++){
            con_nfit[i]=v_fit[i];
        }


        }//end constructor

        int print(){
   
            std::cout<<"### We have "<<n_bins<<" and constraint is applied below bin "<<n_constrain_below<<std::endl;
            std::cout<<"### MC Prediction is : "<<std::endl;
            for(int i=0; i< n_bins; i++)std::cout<<con_pred[i]<<" ";
            std::cout<<std::endl;
            std::cout<<"### Data is : "<<std::endl;
            for(int i=0; i< n_bins; i++)std::cout<<con_data[i]<<" ";
            std::cout<<std::endl;
            std::cout<<"### input inverse matrix is "<<std::endl;
            con_input_inv_Cov.Print();
            std::cout<<"### Calculated Inverse B is "<<std::endl;
            con_inv_B.Print();
            std::cout<<"### Calculated Inverse C is "<<std::endl;
            con_inv_C.Print();
            std::cout<<"### Calculated BC is "<<std::endl;
            con_B.Print();
            std::cout<<"### Calculated N Fit is "<<std::endl;
            for(int i=0; i< n_bins; i++)std::cout<<con_nfit[i]<<" ";
            std::cout<<std::endl;
            return 0;
        }

    };//end SBNconstraint class



};
#endif
