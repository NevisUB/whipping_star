#ifndef SBNCLS_H_
#define SBNCLS_H_

#include <cmath>
#include <vector>
#include <iostream>
#include "SBNspec.h"
#include "SBNchi.h"
#include "SBNconfig.h"

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

#include "ngrid.h"

#include "bayes.h"

#include <gsl/gsl_randist.h>


namespace sbn{


class SBNcls{

	public:

	SBNcls(SBNspec *inh0, SBNspec * inh1, TMatrixD matin) : h0(inh0), h1(inh1), covariance_matrix(matin), chi_h0(*inh0, matin),chi_h1(*inh1,matin){
		which_sample = 0; //default Poisson
        which_mode = 1; //default delta chi
        use_CNP=false;
        maxchival = 210;
		rangen= new TRandom3(0);
        draw_pseudo_from_collapsed = false;
        m_tolerance = 1e-12;
        sample_with_gaussian = false;
colH0 = kRed-7;
        colH1 = kBlue-4;
        has_data_overlayed = false;
        data_result_chis.resize(3);
        data_result_pvals.resize(3);
	}
	SBNcls(SBNspec *inh0, SBNspec * inh1) : h0(inh0), h1(inh1), chi_h0(*inh0),chi_h1(*inh1){
		which_sample = 0; //default Poisson
        which_mode = 1; //default delta chi
        use_CNP = false;
        maxchival = 210;
		rangen= new TRandom3(0);
        draw_pseudo_from_collapsed = false;
        m_tolerance = 1e-12;
        sample_with_gaussian = false;
        colH0 = kRed-7;
        colH1 = kBlue-4;
        has_data_overlayed = false;
        data_result_chis.resize(3);
        data_result_pvals.resize(3);

	}



	SBNspec * h0;
	SBNspec * h1;
	
	SBNchi chi_h0;//previously just chi
	SBNchi chi_h1;

	TMatrixD covariance_matrix;
    bool m_tolerance;

	TRandom3 * rangen;

    bool use_CNP;
    int which_mode;
	int which_sample;
    double maxchival;
	bool draw_pseudo_from_collapsed;
    bool sample_with_gaussian;
    std::vector<std::string> legends;
    int colH0;
    int colH1;
        /****************** Member Functions *************/
    int ReverseColours(){
         colH1 =kRed-7;
         colH0 = kBlue-4;
             return 0;
    }
    int SetLegends(std::string in){
        std::string s = in;
        std::string delimiter = "|";

        size_t pos = 0;
        std::string token;
        while ((pos = s.find(delimiter)) != std::string::npos) {
                token = s.substr(0, pos);
                    legends.push_back(token);
                        s.erase(0, pos + delimiter.length());
        }
        legends.push_back(s);
        return 0;
    };

    int SetTolerance(double epsilon){
        m_tolerance = epsilon;
        std::cout<<"SBNcls::SetTolerance || Set Tolerance of SBNchi's to "<<epsilon<<std::endl;
        chi_h0.setTolerance(epsilon);            
        chi_h1.setTolerance(epsilon);            
    };
    int SetSampleFromCollapsed(){draw_pseudo_from_collapsed = true;};
    int CalcCLS(int,std::string);
	int SetSampleCovariance();
	int SetSamplePoisson();

    int SetGaussianSampling(){sample_with_gaussian = true;};
    double pval2sig(double p);
    double pval2sig1sided(double p);
    double pval2sig2sided(double p);
    int DrawSampleCovariance(std::string);

    int setMode(int);
    int makePlots(CLSresult &h0_result, CLSresult & h1_result, std::string tag,  int which_mode);
    int makePlots(CLSresult &h0_result, CLSresult & h1_result, std::string tag,  int which_mode,double dataval,double datapval);
    int runConstraintTest();

    int compareToRealData(SBNspec * data);
    bool has_data_overlayed;
    std::vector<double> data_result_chis;
    std::vector<double> data_result_pvals;


    int runPi0Tests();

};


};
#endif
