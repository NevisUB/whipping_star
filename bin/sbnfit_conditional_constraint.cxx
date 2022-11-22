#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <vector>
#include <unistd.h>
#include <getopt.h>
#include <cstring>
#include <cstdlib>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TString.h"
#include "TNtuple.h"
#include "TChain.h"
#include "TMath.h"
#include "TSystem.h"
#include "TMatrixT.h"
#include "TRandom.h"
#include "TError.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TLegend.h"

#include "params.h"
#include "SBNconditional.h"
#include "SBNconfig.h"
#include "SBNchi.h"
#include "SBNspec.h"
#include "SBNosc.h"
#include "SBNfit.h"
#include "SBNfit3pN.h"
#include "SBNcovariance.h"
#include "loghelper.h"
log_level_t GLOBAL_LEVEL = LOG_DEBUG;


#define no_argument 0
#define required_argument 1
#define optional_argument 2

using namespace sbn;

/*************************************************************
 *************************************************************
 *		BEGIN sbnfit_conditional_constraint.cxx
 ************************************************************
 ************************************************************/
int main(int argc, char* argv[])
{

    std::string xml = "example.xml";

    /*************************************************************
     *************************************************************
     *		Command Line Argument Reading
     ************************************************************
     ************************************************************/
    const struct option longopts[] =
    {
        {"xml", 		required_argument, 	0, 'x'},
        {"tag", 		required_argument,	0, 't'},
        {"signal", 		required_argument,	0, 's'},
        {"data", 		required_argument,	0, 'd'},
	{"num_channel",		required_argument,      0, 'n'},
	{"bestfit",             required_argument,      0, 'b'},
        {"help", 		no_argument,	        0, 'h'},
    	{"covar",		required_argument,      0, 'c'},
    	{"genie",		required_argument,      0, 'g'},
        {"flat", 		required_argument,      0,'f'},
        {"zero",		no_argument,            0,'z'},
	{"overlaydata", 	no_argument,            0, 'o'},
        {"cmin",		required_argument,      0,'k'},
        {"cmax",		required_argument,      0,'p'},
        {"public",		no_argument, 		0,  'P'},
    };

    int iarg = 0;
    opterr=1;
    int index;
    int num_channel_to_constrain = 2;   //normally we'd like to constrain 2 channels: 1g1p+1g0p
    std::string delimiter=",";          //delimiter for the comma-separated argument
    std::string fitting_subchannel="NONE";
    std::string covar_file = "Stats_Only";
    std::string genie_file = "NONE";
    bool stats_only = true;
    bool use_genie = false;
    bool overlay_data = false;
    bool is_bestfit = false;
    bool public_output = false;
    std::string signal_file;
    std::string data_file;
    std::ostringstream public_txt_tag;
    public_txt_tag.precision(2);
    std::ofstream txt_file;
    txt_file.setf(std::ios::left);

    bool bool_flat_det_sys = false;
    double flat_det_sys_percent = 0.0;

    bool remove_correlations = false;


    double best_fit_value = 0;
    double cmin = 0;
    double cmax = -9;

    //a tag to identify outputs and this specific run. defaults to EXAMPLE1
    std::string tag = "TEST";

    while(iarg != -1)
    {
        iarg = getopt_long(argc,argv, "f:x:s:d:n:t:b:c:g:k:p:zohP", longopts, &index);

        switch(iarg)
        {
              case 'f':
                bool_flat_det_sys = true;
                flat_det_sys_percent = (double)strtod(optarg,NULL);
                break;
              case 'p':
                cmax = (double)strtod(optarg,NULL);
                break;
              case 'k':
                cmin = (double)strtod(optarg,NULL);
                break;
            case 'x':
                xml = optarg;
                break;
            case 'z':
                remove_correlations = true; 
                break;
            case 't':
                tag = optarg;
                break;
     	    case 'c':
    	        covar_file=optarg;
                stats_only = false;
	          break;
	    case 'g':
	    {	use_genie=true;
		std::string comma_separated_argument = optarg; //assume first is the subchannel name, second is genie file
		fitting_subchannel = comma_separated_argument.substr(0, comma_separated_argument.find(delimiter));
		genie_file = comma_separated_argument.substr(comma_separated_argument.find(delimiter)+delimiter.length());
		//genie_file = optarg;
		break;
	    }		
            case 's':
                signal_file = optarg;
                break;
	    case 'd':
		data_file = optarg;
		break;
	    case 'n':
		num_channel_to_constrain = atoi(optarg);
		break; 
 	    case 'b':
	    {   is_bestfit = true;
		std::string comma_separated_argument = optarg; //assume first is the subchannel name, followed by BF value
		fitting_subchannel = comma_separated_argument.substr(0, comma_separated_argument.find(delimiter));
		best_fit_value = std::stod(comma_separated_argument.substr(comma_separated_argument.find(delimiter)+delimiter.length()));
		break;
            }
	    case 'o':
		overlay_data = true;
		break;
	    case 'P':
		public_output = true;
		break;
            case '?':
            case 'h':
                std::cout<<"---------------------------------------------------"<<std::endl;
                std::cout<<"sbnfit_conditional_constraint allows for the plotting of constrained channels and covariance matricies from input root files containing reconstructed variables and covariance matricies; constrained channels should be the first few channels in the xml provided"<<std::endl;
                std::cout<<"---------------------------------------------------"<<std::endl;
                std::cout<<"--- Required arguments: ---"<<std::endl;
                std::cout<<"\t-x\t--xml\t\tInput configuration .xml file for SBNconfig"<<std::endl;
                std::cout<<"\t-s\t--signal\t\tInput signal SBNspec.root file"<<std::endl;
		std::cout<<"\t-d\t--data\t\tInput data root file" << std::endl;
                std::cout<<"\t-c\t--covar\t\tInput Systematic Fractional Covariance"<<std::endl;
                std::cout<<"--- Optional arguments: ---"<<std::endl;
		std::cout<<"\t-n\t--num_channel\t\tNumber of channels to constrain [Default to 2] [Constrained channels are assumed to be in the front of xml]"<<std::endl;
                std::cout<<"\t-t\t--tag\t\tA unique tag to identify the outputs [Default to TEST]"<<std::endl;
		std::cout<<"\t-b\t--bestfit\t\tIf we are constraining the bestfit spectrum. [Default to false] [Input should be: subchannel,BFvalue]"<< std::endl;
                std::cout<<"\t-g\t--genie\t\tInput GENIE Systematic Fractional Covariance, remove Genie correlation between fitting_subchannel and other components. [Format should be: subchannel,GENIEfile]"<<std::endl;
                std::cout<<"\t-f\t--flat\t\tAdd a flat percent systematic to fractional covariance matrix (all channels) (default false, pass in percent, i.e 5.0 for 5\% experimental)"<<std::endl;
                std::cout<<"\t-z\t--zero\t\tZero out all off diagonal elements of the systematics covariance matrix (default false, experimental!)"<<std::endl;
		std::cout<<"\t-o\t--overlaydata\t\tOverlay data points on the constrained plots(default to false)"<<std::endl;
		std::cout<<"\t-P\t--public\t\t Controls if unconstrained, constrained prediction and data will be write out to txt files [Public data release related]" << std::endl;
                std::cout<<"\t--cmax\t max for fractional covariance plot" << std::endl;
                std::cout<<"\t--cmin\t min for fractional covariance plot" << std::endl;
                std::cout<<"\t-h\t--help\t\tThis help menu."<<std::endl;
                std::cout<<"---------------------------------------------------"<<std::endl;

                return 0;
        }
    }
    /*************************************************************
     *************************************************************
     *			Main Program Flow
     ************************************************************
     ************************************************************/
    time_t start_time = time(0);

    std::cout<<"Begining Covariance Plotting for tag: "<<tag<<std::endl;
    std::cout<<"Number of constrained channels: " << num_channel_to_constrain << std::endl; 
    std::cout<<"Loading SBNspec file : "<<signal_file<<" with xml "<<xml<<std::endl;
    SBNspec cv(signal_file,xml, false);
    SBNspec data(data_file,xml, false);
    cv.CollapseVector();
    cv.CalcErrorVector();
    data.CollapseVector();

    //configuration
    //case1: configuration for 2 channels to be constrained, ie. 1g1p and 1g0p
    //int Nsubchannel = 9; 
    //int start_pt_1g0p=cv.hist[0].GetXaxis()->GetNbins();

    //case2: configuration for 1 channel to be constrained, ie 1g1p for near sideband
    //in this case, you should trust result of `1g0p` printed out in the log
    //int Nsubchannel = 0;
    //int start_pt_1g0p=0;
    //int start_pt = start_pt_1g0p + cv.hist[Nsubchannel].GetXaxis()->GetNbins();


    // setup the bins and starting point for the constrained channels
    std::vector<int> channel_bin_index{0};   //starting index for constraind channel in the collapsed covar matrix
    std::vector<int> channel_hist_index{0};
    std::vector<std::string> channel_string; //name of the channel, gonna be title of the plot 
    std::vector<std::string> channel_unit;   //x axis name of the plot
    for(int i=0; i != num_channel_to_constrain; ++i){
	channel_string.push_back(cv.channel_names[i]);
	channel_unit.push_back(cv.channel_units[i]);
	channel_hist_index.push_back(channel_hist_index.back() + cv.num_subchannels[i]);
	channel_bin_index.push_back(channel_bin_index.back() + cv.num_bins[i]);
    }


    //starting point that separates constrained channel and the one constraining
    int start_pt = channel_bin_index.back();

    std::vector<SBNspec> vec_spec(1, cv);
    //get the BF values
    if(is_bestfit){
	 std::cout << "Constraining spectrum at Best-Fit value: " << best_fit_value << " for " << fitting_subchannel << std::endl;
	 vec_spec[0].Scale(fitting_subchannel, best_fit_value);
	 public_txt_tag << "_" <<  fitting_subchannel << "_" << std::fixed << best_fit_value;
    }else{ 
	std::cout << "Constraining CV spectrum " << std::endl;
	public_txt_tag << "_CV";
    }

    std::cout<<"Loading fractional covariance matrix from "<<covar_file<<std::endl;

    TFile * fsys;
    TMatrixT<double> * cov;
    TMatrixT<double>* genie_cov;
    
    if(!stats_only){
        fsys = new TFile(covar_file.c_str(),"read");
        cov = (TMatrixT<double>*)fsys->Get("frac_covariance");
	fsys->Close();
   
        TMatrixT<double> frac_flat_matrix(cv.num_bins_total, cv.num_bins_total);
        if(bool_flat_det_sys){
            std::cout << "RUNNING with flat systematics: " << flat_det_sys_percent << "%!" << std::endl;
            frac_flat_matrix.ResizeTo(cv.num_bins_total,cv.num_bins_total);
            frac_flat_matrix.Zero();
            for(int i=0 ; i< cv.num_bins_total; i++){
                   frac_flat_matrix(i,i)=flat_det_sys_percent*flat_det_sys_percent/10000.;
            }
            (*cov) = (*cov)+(frac_flat_matrix);
        }


        if(remove_correlations){
            std::cout<<"WARNING! We are running in   `Remove All Off Diagional Covariances/Correlations Mode` make sure this is what you want. "<<std::endl;
            for(int i=0; i<cv.num_bins_total;i++){ 
                for(int j=0; j<cv.num_bins_total;j++){ 
                    if(i==j)continue;
                    (*cov)(i,j) =0.0;
                }
            }
        }

       int num_nans = 0;
       for(int i=0; i<cv.num_bins_total;i++){ 
                for(int j=0; j<cv.num_bins_total;j++){ 
                    double val = (*cov)(i,j);
                    if( isinf(val) || isnan(val) || val!=val){
                            (*cov)(i,j) = 0.0;
                            num_nans++;
                    }
                }
            }
        std::cout<<"We have removed "<<num_nans<<" / "<<pow(cv.num_bins_total,2)<<" nans in the fractional covariance matrix"<<std::endl;
    }

    //grab genie covariance matrix
    if(use_genie){
	std::cout << "Running with GENIE correlation between " << fitting_subchannel << " and other subchannels removed !" << std::endl;
	std::cout << "GENIE file used: " << genie_file << std::endl;
	fsys = new TFile(genie_file.c_str(), "read");
	genie_cov  = (TMatrixT<double>*)fsys->Get("individualDir/All_UBGenie_frac_covariance");
        *genie_cov += *((TMatrixT<double>*)fsys->Get("individualDir/AxFFCCQEshape_UBGenie_frac_covariance"));
        *genie_cov += *((TMatrixT<double>*)fsys->Get("individualDir/DecayAngMEC_UBGenie_frac_covariance"));
        *genie_cov += *((TMatrixT<double>*)fsys->Get("individualDir/NormCCCOH_UBGenie_frac_covariance"));
        *genie_cov += *((TMatrixT<double>*)fsys->Get("individualDir/NormNCCOH_UBGenie_frac_covariance"));
        *genie_cov += *((TMatrixT<double>*)fsys->Get("individualDir/RPA_CCQE_UBGenie_frac_covariance"));
        *genie_cov += *((TMatrixT<double>*)fsys->Get("individualDir/Theta_Delta2Npi_UBGenie_frac_covariance"));
        *genie_cov += *((TMatrixT<double>*)fsys->Get("individualDir/VecFFCCQEshape_UBGenie_frac_covariance"));
        *genie_cov += *((TMatrixT<double>*)fsys->Get("individualDir/XSecShape_CCMEC_UBGenie_frac_covariance"));
	fsys->Close();

	//remove nan's
	for(int i=0; i<cv.num_bins_total;i++){
             for(int j=0; j<cv.num_bins_total;j++){
                    double val = (*genie_cov)(i,j);
                    if( isinf(val) || isnan(val) || val!=val){
                            (*genie_cov)(i,j) = 0.0;
                    }       
             }   
        }  


	*cov -= *genie_cov;

	 //remove correlation in genie
	SBNspec temp_cv(cv);
        temp_cv.Keep(fitting_subchannel, 1.0);
        temp_cv.CalcFullVector();
        std::vector<double> temp_full = temp_cv.full_vector;
        temp_cv=cv;
        temp_cv.Scale(fitting_subchannel, 0.0);
        std::vector<double> temp_others = temp_cv.full_vector;

        for(int i=0; i<temp_full.size(); i++){
           for(int j=0;j<temp_full.size(); j++){
                if(temp_full[i]*temp_full[j]==0 && temp_others[i]*temp_others[j]==0) (*genie_cov)(i,j)=0.0;
           }
        }

	*cov += *genie_cov;
    }


    TFile *fout = new TFile(("Constraint_"+tag+"_output.root").c_str(),"recreate");


    int which_spec_index=0;
    //loop over the spectra vector, and calculating the constrained chi2
    for(SBNspec& fspec:vec_spec){
	    fspec.CollapseVector();
	    fspec.CalcErrorVector();
	    std::cout << "\n\n\nOn spec:  " << ( is_bestfit ? "BF" : "CV" ) << std::endl;;
	    SBNchi SigChi(fspec, *cov);
	    if(stats_only) SigChi.is_stat_only=true;
	//    SigChi.SetFracPlotBounds(cmin,cmax);
	//    SigChi.PrintMatricies(tag);
	//    sig.WriteOut(tag);


	    double chi_constrain=0;
	    double chi_nueoriginal=0;
	    double chi_total=0;

	    TMatrixT<double> collapsed_covar(cv.num_bins_total_compressed, cv.num_bins_total_compressed);
	    SigChi.FillCollapsedCovarianceMatrix(&collapsed_covar);   //gq: collapsed_covar now has both nue and numu stats
	    //everything from here down is collapsed.

	    // Remove the stats for nue
	    for(int i=0; i< cv.num_bins_total_compressed;i++){
	         if(i < start_pt) collapsed_covar(i,i) -= fspec.collapsed_vector.at(i);  //use the MC stats for constrained matrix
	    }
   

	    // ****************************calculate chi constrain***************************************
	    //constrained covar matrix
	    std::vector<TMatrixT<double>> v_mat =  sbn::splitCovariance(collapsed_covar, start_pt);
	    TMatrixT<double> constrained_mat = sbn::getConstrainedCovariance(v_mat);
	    //constrained nue prediction
	    std::vector<double> constrained_pred = sbn::getConstrainedPrediction(v_mat, fspec.collapsed_vector, data.collapsed_vector, start_pt);

	    TMatrixT<double> full_constrained_mat = constrained_mat;
	    for(int i=0; i<start_pt; i++){
		full_constrained_mat(i,i) +=  ( data.collapsed_vector.at(i) >0.001 ? 3.0/(1.0/data.collapsed_vector.at(i) +  2.0/constrained_pred[i])  : constrained_pred[i]/2.0 );
	    }

	    TMatrixT<double> constrained_invert = SigChi.InvertMatrix(full_constrained_mat);
	    for(int i=0; i< start_pt; i++){
		for(int j=0; j< start_pt; j++){
			chi_constrain += constrained_invert(i,j)*(data.collapsed_vector.at(i)-constrained_pred[i])*(data.collapsed_vector.at(j)-constrained_pred[j]);
		}
	    }
	   
	    // *************************calculate chi nue original only**************************************
	    TMatrixT<double> original_mat = v_mat[0];
	    TMatrixT<double> original_full_mat = original_mat;
	    for(int i=0; i<start_pt; i++){
                original_full_mat(i,i) +=  ( data.collapsed_vector.at(i) >0.001 ? 3.0/(1.0/data.collapsed_vector.at(i) +  2.0/fspec.collapsed_vector.at(i))  : fspec.collapsed_vector.at(i)/2.0 );
            }
	    TMatrixT<double> original_invert = SigChi.InvertMatrix(original_full_mat);
	   for(int i=0; i< start_pt; i++){
                for(int j=0; j< start_pt; j++){
                        chi_nueoriginal += original_invert(i,j)*(data.collapsed_vector.at(i)-fspec.collapsed_vector.at(i))*(data.collapsed_vector.at(j)-fspec.collapsed_vector.at(j));
                }
            }
	    // **************************calculate original nue+numu chi2 ********************************
	    TMatrixT<double> total_collapsed_mat = SigChi.FillSystMatrix(*cov, fspec.full_vector, fspec.full_err_vector, true);
	    TMatrixT<double> total_mat = SigChi.AddStatMatrixCNP(&total_collapsed_mat, fspec.collapsed_vector, data.collapsed_vector);
	    TMatrixT<double> total_invert = SigChi.InvertMatrix(total_mat);
	    for(int i=0; i<cv.num_bins_total_compressed; i++){
		for(int j=0; j< cv.num_bins_total_compressed; j++){
		    chi_total+=total_invert(i,j)*(data.collapsed_vector.at(i)-fspec.collapsed_vector.at(i))*(data.collapsed_vector.at(j)-fspec.collapsed_vector.at(j));
		}
	    }
     
	     //print out the comparison
	     //SigChi.DrawComparisonIndividual(fspec, data, total_collapsed_mat, "Spec"+std::to_string(which_spec_index)+"_VS_Data", false);

	    // **************************done calculating chi2 *******************************************

	    std::cout << "\n=============== Overall Summary ======================" << std::endl;
	    std::cout << "Numu-nue side-by-side chi2 value is " << chi_total << std::endl;
	    std::cout << "Nue only original  chi2 value is " << chi_nueoriginal << std::endl;
	    std::cout << "Nue constrained  chi2 value is " << chi_constrain << std::endl;
	    std::cout << "=============== Overall Summary ======================\n" << std::endl;


            // *********************** loop over each constrained channel *********************************
	    for(int i = 0; i != num_channel_to_constrain; ++i){
		TMatrixT<double> sub_unconstrain = collapsed_covar.GetSub(channel_bin_index[i], channel_bin_index[i+1] -1, channel_bin_index[i], channel_bin_index[i+1] -1);
		TMatrixT<double> sub_constrain = constrained_mat.GetSub(channel_bin_index[i], channel_bin_index[i+1] -1, channel_bin_index[i], channel_bin_index[i+1] -1);


		double intrinsic_error = std::inner_product(fspec.collapsed_err_vector.begin() +channel_bin_index[i], fspec.collapsed_err_vector.begin() +channel_bin_index[i+1], fspec.collapsed_err_vector.begin() + channel_bin_index[i], 0.0);
		double overall_sys_unconstrain = sub_unconstrain.Sum();
		double overall_sys_constrain = sub_constrain.Sum();

		std::cout << "====== Channel: " << channel_string[i] << ", MC intrinsic error: " << sqrt(intrinsic_error) << ", overall unconstrained error: " << sqrt(overall_sys_unconstrain) << ", overall constrained error: " << sqrt(overall_sys_constrain) << "\n\t\t unconstrained events: " << std::accumulate(fspec.collapsed_vector.begin()  +channel_bin_index[i], fspec.collapsed_vector.begin() +channel_bin_index[i+1], 0.0) << ", constrained events: " << std::accumulate(constrained_pred.begin() + channel_bin_index[i], constrained_pred.begin()+channel_bin_index[i+1], 0.0) << std::endl;


             	//calculate chi2 of individual histogram
		double sub_chi_constrain = 0, sub_chi_unconstrain = 0;
		TMatrixT<double> full_sub_unconstrain = sub_unconstrain;
		TMatrixT<double> full_sub_constrain = sub_constrain;

             	for(int j= channel_bin_index[i], row=0; j!=channel_bin_index[i+1]; ++j, ++row){
                     full_sub_constrain(row, row) +=  ( data.collapsed_vector.at(j) >0.001 ? 3.0/(1.0/data.collapsed_vector.at(j) +  2.0/constrained_pred[j])  : constrained_pred[j]/2.0 );
                     full_sub_unconstrain(row, row) +=  ( data.collapsed_vector.at(j) >0.001 ? 3.0/(1.0/data.collapsed_vector.at(j) +  2.0/fspec.collapsed_vector.at(j))  : fspec.collapsed_vector.at(j)/2.0 );
                }

                TMatrixT<double> sub_constrained_invert = SigChi.InvertMatrix(full_sub_constrain);
		TMatrixT<double> sub_unconstrain_invert = SigChi.InvertMatrix(full_sub_unconstrain);
                for(int j=channel_bin_index[i], row=0; j!=channel_bin_index[i+1]; ++j, ++row){
                    for(int k=channel_bin_index[i], col=0; k!=channel_bin_index[i+1]; ++k, ++col){
                        sub_chi_constrain += sub_constrained_invert(row, col)*(data.collapsed_vector.at(j)-constrained_pred[j])*(data.collapsed_vector.at(k)-constrained_pred[k]);
			sub_chi_unconstrain += sub_unconstrain_invert(row, col)*(data.collapsed_vector.at(j) - fspec.collapsed_vector.at(j))*(data.collapsed_vector.at(k) - fspec.collapsed_vector.at(k));
                    }
                } 

   	       // ***************************start drawing histograms******************************************
 	       TH1D* h_original=(TH1D*)(cv.hist[channel_hist_index[i]]).Clone(); h_original->Reset();
	       TH1D* h_constrain=(TH1D*)h_original->Clone();  h_constrain->Reset();
	       TH1D* h_data=(TH1D*)h_original->Clone();  h_data->Reset();

	       for(int j=channel_bin_index[i], bin = 1; j< channel_bin_index[i+1]; ++j, ++bin){
		   h_original->SetBinContent(bin, fspec.collapsed_vector.at(j));
		   h_constrain->SetBinContent(bin, constrained_pred[j]);
		   h_data->SetBinContent(bin, data.collapsed_vector.at(j));

		   //reset error bars
		   h_original->SetBinError(bin, sqrt(sub_unconstrain(bin-1, bin-1)));
		   h_constrain->SetBinError(bin, sqrt(sub_constrain(bin-1, bin-1)));

		   h_data->SetBinErrorOption(TH1::kPoisson);  //use the Poisson error for data point
		   //h_data->SetBinError(bin, sqrt(h_data->GetBinContent(bin)));
	       }


	 	// write out txt file for HEPdata public release
	       std::map<std::string, double> norm_bin_width = {{"1g1p", 0.1}, {"1g0p",0.05}};
  	       bool do_norm_bin = false;
	       if(public_output){

		   std::string filename = tag + "_ConditionalConstraint_" + channel_string[i] + public_txt_tag.str() + ".txt";
		   txt_file.open(filename, std::ofstream::out | std::ofstream::trunc);

		   //print out some helper info
		   if(do_norm_bin) txt_file << "Event Counts Per " << norm_bin_width[channel_string[i]]  << " GeV\n" << std::endl;
		   else txt_file << "Event Counts\n" << std::endl; 
		   txt_file << std::setw(15) << "LowBinEdge" << std::setw(15) << "UpBinEdge" << std::setw(15) << "Unconstrained" << std::setw(20) << "UnconstrainedError" << std::setw(15) << "Constrained" << std::setw(20) << "ConstrainedError" << std::setw(15) << "Data" << std::setw(20) << "PoissonErrorUp" << std::setw(20) << "PoissonErrorLow" << std::endl; 

		   for(int bin = 1; bin <= h_original->GetNbinsX(); ++bin){
			if(do_norm_bin){
			    double bin_scale = h_original->GetBinWidth(bin)/norm_bin_width[channel_string[i]];
			    txt_file << std::setw(15) << h_original->GetBinLowEdge(bin) << std::setw(15) << h_original->GetBinLowEdge(bin) + h_original->GetBinWidth(bin)  << std::setw(15) << h_original->GetBinContent(bin)/bin_scale  << std::setw(20) << h_original->GetBinError(bin)/bin_scale  << std::setw(15) << h_constrain->GetBinContent(bin)/bin_scale << std::setw(20) << h_constrain->GetBinError(bin)/bin_scale  << std::setw(15) << h_data->GetBinContent(bin)/bin_scale  << std::setw(20) << h_data->GetBinErrorUp(bin)/bin_scale << std::setw(20) << h_data->GetBinErrorLow(bin)/bin_scale << std::endl;
			}else{
			    txt_file << std::setw(15) << h_original->GetBinLowEdge(bin) << std::setw(15) << h_original->GetBinLowEdge(bin) + h_original->GetBinWidth(bin)  << std::setw(15) << h_original->GetBinContent(bin)  << std::setw(20) << h_original->GetBinError(bin)  << std::setw(15) << h_constrain->GetBinContent(bin) << std::setw(20) << h_constrain->GetBinError(bin)  << std::setw(15) << h_data->GetBinContent(bin)  << std::setw(20) << h_data->GetBinErrorUp(bin) << std::setw(20) << h_data->GetBinErrorLow(bin) << std::endl;
			}
		   }

		   txt_file.close();
  	       }


	       TCanvas* c=new TCanvas(Form("c_%d_%d",which_spec_index, i), "c");
	       gStyle->SetOptStat(0);
	       //gStyle->SetErrorX();
	       TLegend* leg=new TLegend(0.45,0.7,0.89,0.89);
	       leg->SetBorderSize(0);
	       leg->SetTextSize(0.032);
	       h_original->SetLineColor(kMagenta+3);
	       h_original->SetLineWidth(2);
	       h_original->SetLineStyle(kDashed);
	       TH1D* h_original_copy= (TH1D*)h_original->Clone();
	       h_original->SetFillColorAlpha(kMagenta -10, 0.8);
	       //h_original->SetFillColor(kMagenta -10);
	       h_original->SetTitle(Form("%s;%s; Events", channel_string[i].c_str(), channel_unit[i].c_str()));
	       //h_1g1p_original->SetFillStyle(4050);  //only useful for TPad
               h_original->SetMaximum(std::max(std::initializer_list<double>{h_original->GetMaximum(), h_data->GetMaximum(), h_constrain->GetMaximum()})*2.5);
	       h_original->SetMinimum(0);

	       h_constrain->SetLineColor(kGreen+3);
	       h_constrain->SetLineWidth(2);
	       TH1D* h_constrain_copy = (TH1D*)h_constrain->Clone();
	       h_constrain->SetFillColorAlpha(kGreen-7, 0.9);
	       //h_constrain->SetFillColor(kGreen -10);
	       gStyle->SetHatchesLineWidth(3);
	       h_constrain->SetFillStyle(3345);

	       h_original->Draw("E2");
	       h_original_copy->Draw("HIST SAME");
	       h_constrain->Draw("E2same");
	       h_constrain_copy->Draw("HIST SAME");
	       leg->AddEntry(h_original, Form("%s MC prediction - %.1f", channel_string[i].c_str(), h_original->Integral()), "LF");
	       leg->AddEntry(h_constrain, Form("%s Constrained prediction - %.1f",channel_string[i].c_str(), h_constrain->Integral()), "LF");


	       //legend to print chi^2 and pvalues
	       TLegend* chi_leg = new TLegend(0.11, 0.76, 0.35, 0.89);
	       //chi_leg->SetNColumns(2);
	       chi_leg->SetFillStyle(0);
               chi_leg->SetBorderSize(0);
	       chi_leg->SetTextSize(0.03);
	       chi_leg->AddEntry(h_original, Form("#chi^{2}_{CNP}/ndof: %.2f/%d, P^{Wilks}_{val}: %.3f", sub_chi_unconstrain, channel_bin_index[i+1] - channel_bin_index[i], TMath::Prob(sub_chi_unconstrain, channel_bin_index[i+1] - channel_bin_index[i]) ));
	       chi_leg->AddEntry(h_constrain, Form("#chi^{2}_{CNP}/ndof: %.2f/%d, P^{Wilks}_{val}: %.3f", sub_chi_constrain, channel_bin_index[i+1] - channel_bin_index[i],  TMath::Prob(sub_chi_constrain, channel_bin_index[i+1] - channel_bin_index[i]) ));

	       if(overlay_data){
		   h_data->SetMarkerStyle(20);
		   h_data->SetMarkerColor(kBlack);
		   h_data->SetMarkerSize(1.2);
		   h_data->SetLineWidth(2);
		   h_data->SetLineColor(kBlack);
		   h_data->Draw("E1X0same");
		   leg->AddEntry(h_data, Form("%s Data - %d",channel_string[i].c_str(), (int)h_data->Integral()), "ep");
	       }
	       leg->Draw();
	       c->Update();
	       fout->cd();
	       if(is_bestfit){
		   c->SaveAs((tag+"_BF_"+channel_string[i]+"_constrain_comparison.pdf").c_str(), "pdf");
	           c->Write(Form("BF_%s_comparison",channel_string[i].c_str()));
	       }else{
	           c->SaveAs((tag+"_CV_"+channel_string[i]+"_constrain_comparison.pdf").c_str(), "pdf");
		   c->Write(Form("CV_%s_comparison",channel_string[i].c_str()));
	       }

	       if(overlay_data){
	       	   chi_leg->Draw();
		   c->Update();
		   if(is_bestfit){
                   	c->SaveAs((tag+"_Wchi_BF_"+channel_string[i]+"_constrain_comparison.pdf").c_str(), "pdf");
                   	c->Write(Form("BF_Wchi_%s_comparison",channel_string[i].c_str()));
               	   }else{
                   	c->SaveAs((tag+"_Wchi_CV_"+channel_string[i]+"_constrain_comparison.pdf").c_str(), "pdf");
                   	c->Write(Form("CV_Wchi_%s_comparison",channel_string[i].c_str()));
                   }	
	       }


   	   }  //end of constrained channel loop

	   ++which_spec_index;
    }
/*    TFile *fout = new TFile(("Constraint_"+tag+"_output.root").c_str(),"recreate");
    fout->cd();
    constrained_mat.Write("full_covariance");
    fout->Close();

    TH1D h1g1p  = sig.hist[0];  
    TH1D * a1g1p  =(TH1D*)h1g1p.Clone("after");
    TH1D h1g0p  = sig.hist[12];  
    TH1D * a1g0p  =(TH1D*)h1g0p.Clone("after2");

    h1g1p.Clear();
    h1g0p.Clear();

    a1g0p->Clear();
    a1g1p->Clear();


    int i = 0;
    for(i=0; i< h1g1p.GetNbinsX();i++){
            h1g1p.SetBinContent(i+1,0);
            h1g1p.SetBinError(i+1,100*sqrt(collapsed_covar(i,i))/sig.collapsed_vector.at(i));
    }
     for(int j=0; j< h1g0p.GetNbinsX();j++){
            h1g0p.SetBinContent(j+1,0);
            h1g0p.SetBinError(j+1,100*sqrt(collapsed_covar(i,i))/sig.collapsed_vector.at(i));
            i++;
    }
 
    i = 0;
    for(i=0; i< a1g1p->GetNbinsX();i++){
            a1g1p->SetBinContent(i+1,0);
            a1g1p->SetBinError(i+1,100*sqrt(constrained_mat(i,i))/sig.collapsed_vector.at(i));
    }
     for(int j=0; j< a1g0p->GetNbinsX();j++){
            a1g0p->SetBinContent(j+1,0);
            a1g0p->SetBinError(j+1,100*sqrt(constrained_mat(i,i))/sig.collapsed_vector.at(i));
            i++;
    }
    

    TCanvas *c = new TCanvas("c","c",1400,600);
    c->Divide(2,1);
    TPad*p1=(TPad*)c->cd(1);
   
    TLegend l1g1p(0.1,0.8,0.89,0.89);
    l1g1p.SetNColumns(2);
    l1g1p.SetLineColor(kWhite);
    l1g1p.SetLineWidth(0);
    l1g1p.SetFillStyle(0);

    h1g1p.Draw("hist");
    h1g1p.SetFillColor(kRed-7);
    h1g1p.SetLineColor(kRed-7);
    h1g1p.Draw("E2 same");
    h1g1p.SetTitle("1#gamma1p Final Selection");
    h1g1p.GetXaxis()->SetTitle("Reconstructed Shower Energy [GeV]");
    h1g1p.GetYaxis()->SetTitle("Fractional Error [\%]");
    h1g1p.SetFillStyle(3333);
    h1g1p.SetMaximum(60);
    h1g1p.SetMinimum(0);
    
    l1g1p.AddEntry(&h1g1p,"Flux & XS Systematics","F");

    a1g1p->SetFillColor(kRed-7);
    a1g1p->SetLineColor(kRed-7);
    a1g1p->Draw("E2 same");

    l1g1p.AddEntry(a1g1p,"Constrained Systematics","F");
    l1g1p.Draw();
    p1->RedrawAxis();


    TPad*p2 = (TPad*)c->cd(2);
    
    TLegend l1g0p(0.1,0.8,0.89,0.89);
    l1g0p.SetNColumns(2);
    l1g0p.SetLineColor(kWhite);
    l1g0p.SetLineWidth(0);
    l1g0p.SetFillStyle(0);


    h1g0p.Draw("hist");
    h1g0p.SetFillColor(kBlue-7);
    h1g0p.SetLineColor(kBlue-7);
    h1g0p.Draw("E2 same");
    h1g0p.SetTitle("1#gamma0p Final Selection");
    h1g0p.GetXaxis()->SetTitle("Reconstructed Shower Energy [GeV]");
    h1g0p.GetYaxis()->SetTitle("Fractional Error [\%]");
    h1g0p.SetFillStyle(3333);
    h1g0p.SetMaximum(60);
    h1g0p.SetMinimum(0);
    
    a1g0p->SetFillColor(kBlue-7);
    a1g0p->SetLineColor(kBlue-7);
    a1g0p->Draw("E2 same");
    
    l1g0p.AddEntry(&h1g0p,"Flux & XS Systematics","F");
    l1g0p.AddEntry(a1g0p,"Constrained Systematics","F");
    l1g0p.Draw();
    
    p2->RedrawAxis();
    c->SaveAs(("Constraint_"+tag+".pdf").c_str(),"pdf");
*/
    
    fout->Close();
    std::cout << "Total wall time: " << difftime(time(0), start_time)/60.0 << " Minutes.\n";
    return 0;
}
