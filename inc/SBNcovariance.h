#ifndef SBNCOVARIANCE_H_
#define SBNCOVARIANCE_H_

#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <map>
#include <ctime>

#include "SBNspec.h"
#include "SBNconfig.h"
#include "SBNchi.h"

#include "TMatrixDEigen.h"
#include "TMatrixDSymEigen.h"
#include "TMatrixD.h"
#include "TMatrixT.h"
#include "TDirectory.h"
#include "TVectorD.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include "TNtuple.h"
#include "TLine.h"

#include "TText.h"
#include "TROOT.h"
#include "TRint.h"
#include "TGraph.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "THnSparse.h"
#include "TStopwatch.h"
#include "TTreeFormula.h"
#include "TCanvas.h"

#include "params.h"
#include <sys/stat.h> 


namespace sbn{

  int analyzeCovariance(std::string xml, std::string signal_file, std::string tag, std::vector<std::string> covar_files, std::vector<std::string> covar_names);

  class SBNcovariance : public SBNconfig {

    bool is_small_negative_eigenvalue;

  public:
		
    SBNspec spec_central_value;	

    std::vector<TH2D> h_spec2d_CV;
    std::vector<TH2D> h_spec2d_sig;
	
    SBNcovariance(std::string xmlname);
    SBNcovariance(std::string xmlname, bool);
    SBNcovariance(std::string xmlname, std::string);


    void ShapeOnlyProcessing();
    int FormCovarianceMatrix(std::string tag);
    void WriteOutVariation(std::string signal_tag) ; //write out variation spectrum in two histograms: signal/background
    int PrintMatricies(std::string tag);
    int PrintMatricies(std::string tag,bool);
    std::string output_tag;
    void GrabSubMatrix(std::string filename, std::string matrix_name, const std::vector<std::string>& channel_list);

    int plot_one(TMatrixD matrix, std::string tag, TFile *fin,bool,bool);
    int plot_one(TMatrixD matrix, std::string tag, TFile *fin,bool,bool,bool);
    int qualityTesting();
    virtual bool EventSelection(int file);
    virtual int FillHistograms(int file, int uni, double wei);
    
    std::vector<SBNspec> multi_sbnspec;
    std::vector<std::vector<double>> multi_vecspec;

    TMatrixD full_covariance;
    TMatrixD frac_covariance;
    TMatrixD full_correlation;
    TMatrixD full_collapse;
    TMatrixD full_constrain;

    std::vector<TMatrixD> vec_full_covariance;
    std::vector<TMatrixD> vec_frac_covariance;
    std::vector<TMatrixD> vec_full_correlation;
   std::vector<TMatrixD> vec_full_collapse;
      std::vector<TMatrixD> vec_full_constrain;

    TMatrixD input_frac_covariance;

    //for plotting, hsa been superseeded by PrintMatricies a bit
    TH2D * hist_frac_cov;
    TH2D * hist_full_cor;
    TH2D * hist_full_cov;

    //Some checks on montecarlos
    double tolerence_positivesemi;
    int universes_used;
    double abnormally_large_weight;

    //Multisim input variables
    std::vector<std::string> variations;
    std::vector<std::vector<double> > vars;

    std::vector<int> num_universes_per_variation;
    std::map<int, std::string> map_universe_to_var;
    std::vector<int> vec_universe_to_var;
    std::map<std::string, int> map_var_to_num_universe;

    std::map<std::string, int> map_var_to_matrix;
    std::vector<int> vec_var_to_matrix;

    int PrintVariations(std::string tag);
    int PrintVariations_2D(std::string tag);
    std::string bnbcorrection_str;

    int num_files;

    std::vector<TFile*> files;
    std::vector<TTree*> trees;

    std::vector<int> nentries;
    std::vector<TBranch*>* branch_weight;
    std::vector<std::map<std::string, std::vector<eweight_type> >* > f_weights;

    std::vector<std::vector<int> > vars_i;
    std::vector<std::vector<double> > vars_d;

    std::vector<std::string> variations_to_use;
    std::map<std::string,bool> m_variations_to_use;

    // In testing 2D fits with 4D covariance..
    SBNspec template_spec;	
    SBNspec spec_central_value2;	
    THnSparseD *frac_covariance_4d;
    TMatrixD full_covariance2D;
    TMatrixD frac_covariance2D;
    TMatrixD full_correlation2D;
    std::vector<std::vector<double> > multi_vecspec2D;
    std::vector<double> vecspec2DCV;
    TStopwatch watch;
    
    
    void ProcessEvent(const std::map<std::string, std::vector<eweight_type> >& thisfWeight,   size_t fileid,     int entryid);
   
    std::ofstream variation_scale;
    std::ofstream sorted_variation_scale;
    std::ofstream correlation_scale;
    std::ofstream correlation_constraint;
    std::ofstream ratio_con;
    
    std::vector<std::string> buildWeightMaps();
    std::vector<std::vector<TTreeFormula*>> m_variation_weight_formulas;
    std::vector<int> m_variation_modes;

    int DoConstraint(int which_signal, int which_constraint);
    std::vector<double> DoConstraint(int which_signal, int which_constraint, std::string tag);
    std::vector<double> DoConstraint(int which_signal, int which_constraint, std::string tag,int which_var);
    std::vector<double> DoConstraint_test(int which_signal, int which_constraint, std::string tag);

    int make_tables(std::string tag);
    std::ofstream constraint_table;
    
     bool check_inversion(TMatrixD matrix, TMatrixD inverse, int num_bins){
     double tol=1e-6;
     TMatrixT<double>* m_check1=new TMatrixT<double>(num_bins,num_bins);
     m_check1->Mult(matrix,inverse);
     std::cout<<"made check 1"<<std::endl;
     TMatrixT<double>* m_check2=new TMatrixT<double>(num_bins,num_bins);
     m_check1->MultT(matrix,inverse);
     m_check2->MultT(inverse,matrix);
     m_check1->Print();
     m_check2->Print();
     std::cout<<"made unity matrices"<<std::endl;
     std::cout<<"made unity matrices"<<std::endl;
     for(int i=0; i<num_bins;i++)
       {
	 for(int j=0; j<num_bins;j++)
	   {
	     if (i == j){
	       if ( abs((*m_check1)(i,j)-1) < tol && abs((*m_check2)(i,j)-1) < tol ) continue;
						   
	       else return false;
	     }
	     else{
	       if (abs((*m_check1)(i,j)) < tol && abs((*m_check2)(i,j)) < tol ) continue;
						   
	       else return false;
	     }
	   }
       }
     return true;
     }

   

  };


};
#endif
