#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <vector>
#include <unistd.h>
#include <getopt.h>
#include <float.h>
#include <cstring>
#include <gsl/gsl_math.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TH3.h"
#include "TMultiGraph.h"
#include "TString.h"
#include "TNtuple.h"
#include "TChain.h"
#include "TMath.h"
#include "TSystem.h"
#include "TMatrixT.h"
#include "TObjArray.h"
#include "TList.h"
#include "TRandom.h"
#include "TError.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TROOT.h"
#include "TMarker.h"

#include "params.h"
#include "SBNconfig.h"
#include "SBNchi.h"
#include "SBNspec.h"
#include "SBNosc.h"
#include "SBNfit.h"
#include "SBNfit3pN.h"
#include "SBNcovariance.h"
#include "SBNfeld.h"
#include "loghelper.h"
log_level_t GLOBAL_LEVEL = LOG_DEBUG;



#define no_argument 0
#define required_argument 1
#define optional_argument 2

using namespace sbn;

/*************************************************************
 *************************************************************
 *		BEGIN sbnfit_make_covariance.cxx
 ************************************************************
 ************************************************************/
int main(int argc, char* argv[])
{

    std::string xml = "oscillate_example.xml";

    /*************************************************************
     *************************************************************
     *		Command Line Argument Reading
     ************************************************************
     ************************************************************/
    const struct option longopts[] =
    {
        {"xml", 		required_argument, 	0, 'x'},
        {"printall", 		no_argument, 		0, 'p'},
        {"stat", 		no_argument, 		0, 's'},
        {"montecarlo", 		required_argument,	0, 'm'},
        {"data", 		required_argument,	0, 'd'},
	{"interpolation",       required_argument,      0, 'i'},
	{"covarmatrix",         required_argument,      0, 'c'},
	{"flat",        required_argument, 0 ,'f'},
        {"randomseed",        required_argument, 0 ,'r'},
        {"help", 		no_argument,	0, 'h'},
        {0,			    no_argument, 		0,  0},
    };

    int iarg = 0;
    opterr=1;
    int index;

    //a tag to identify outputs and this specific run. defaults to EXAMPLE1
    std::string tag;
    bool bool_stat_only = false;
    int interpolation_number = -99;  //number of points for chi2 value interpolation
    double random_number_seed = -1;

    bool input_data = false;
    std::string data_filename;  //root file containing data/mc spectra
    std::string mc_filename;
    std::string covmatrix_file;  //root file containing covariance matrix

    bool bool_flat_det_sys = false;
    double flat_det_sys_percent = 0.0;

    while(iarg != -1)
    {
        iarg = getopt_long(argc,argv, "d:x:m:r:p:i:f:sh", longopts, &index);

        switch(iarg)
        {
            case 'x':
                xml = optarg;
                break;
            case 'd':
                input_data = true;
                data_filename = optarg;
                break;
	    case 'm':
		mc_filename = optarg;
		break;
	    case 'i':
		interpolation_number = (int)strtod(optarg, NULL);
		break;
	    case 'c':
		covmatrix_file = optarg;
		break;
            case 'f':
                bool_flat_det_sys = true;
                flat_det_sys_percent = (double)strtod(optarg,NULL);
                break;
            case 'r':
                random_number_seed = (double)strtod(optarg,NULL);
                std::cout<<"Reading in random seed argument: "<<random_number_seed<<std::endl;
                break;
            case 's':
                bool_stat_only = true;
                break;
            case '?':
            case 'h':
                std::cout<<"---------------------------------------------------"<<std::endl;
                std::cout<<"sbnfit_fraction_fit is a work in progress."<<std::endl;
                std::cout<<"---------------------------------------------------"<<std::endl;
                std::cout<<"--- Required arguments: ---"<<std::endl;
                std::cout<<"\t-x\t--xml\t\tInput configuration .xml file for SBNconfig"<<std::endl;
                std::cout<<"\t-d\t--data\t\tInput observed data for a global scan"<<std::endl;
                std::cout<<"\t-m\t--montecarlo\t\tInput monte carlo for global scan"<<std::endl;
		std::cout<<"\t-i\t--interpolation\t\tInput number of points for interpolation"<< std::endl;
		std::cout<<"\t-c\t--covariance matrix\t\tInput syst covariance matrix for global scan"<< std::endl;
                std::cout<<"\t-f\t--flat\t\t Input flat systematic fractional covariance matrix"<<std::endl;
                std::cout<<"--- Optional arguments: ---"<<std::endl;
                std::cout<<"\t-s\t--stat\t\tStat only runs"<<std::endl;
                std::cout<<"\t-r\t--randomseed\t\tRandomNumber Seed (default from machine)"<<std::endl; 
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

    std::cout<<"Begining Fraction Fit"<<std::endl;

    NGrid mygrid;

    //now only available for 2 subchannels only
    mygrid.AddDimension("NCPi0Coh", 0, 5, 0.1);   //0.1 FULL
    mygrid.AddDimension("NCPi0NotCoh", 0., 3, 0.1);   //0.1 FULL
    mygrid.AddDimension("NCDeltaRadOverlaySM", 0., 4, 0.2);   //0.1 FULL


    std::cout << "Fraction Fit|| "<< "\tStart initializing MC and data spectrum" << std::endl;

    //initialize the MC spectrum
    SBNspec mc_spec(mc_filename, xml);

    //initlaize the data spectrum
    SBNspec data_spec(data_filename, xml);
    data_spec.CollapseVector();  //collapse full vector
    
    //compare plots before the fit
    tag = "before_fit";
    TMatrixT<double> tt(data_spec.num_bins_total_compressed,data_spec.num_bins_total_compressed);
    tt.Zero();
    mc_spec.CompareSBNspecs(tt,&data_spec, tag);
    tag.clear();


    std::cout << "Fraction Fit||" <<  "\tInitialize fractional systematric covariance matrix" << std::endl;
    //initialize covariance matrix
    TMatrixT<double> frac_syst_matrix(mc_spec.num_bins_total, mc_spec.num_bins_total);


    if(bool_stat_only){
    	frac_syst_matrix.Zero();
        std::cout<<"RUNNING Stat Only!"<<std::endl;
    }else if(bool_flat_det_sys){
	std::cout << "RUNNING with flat systematics: " << flat_det_sys_percent << "%!" << std::endl;
	frac_syst_matrix.Zero();
	//set up flat fractional syst covariance
	for(int i=0 ; i< mc_spec.num_bins_total; i++)
		for(int j=0;j<mc_spec.num_bins_total; j++)
			frac_syst_matrix(i,j)=flat_det_sys_percent*flat_det_sys_percent/10000.;
    }
    else{
	std::cout<< "RUNNING with Systematics!" << std::endl;
	frac_syst_matrix.Zero();

	std::cout << "Open covariance matrix root file: " << covmatrix_file << std::endl;
	TFile* f_covar = new TFile(covmatrix_file.c_str(), "read");
	TMatrixT<double>* m_temp;
	m_temp = (TMatrixD*)f_covar->Get("frac_matrix_name");
	frac_syst_matrix = *m_temp;
	f_covar->Close();	
    }



   std::cout<< "Fraction Fit||" << "\tGrab info of the grid"<< std::endl;
   //check grid size
   std::vector<std::vector<double>> grid = mygrid.GetGrid();
   if(grid.size() != mygrid.f_num_total_points){
	std::cout <<  "the number of points don't match: something wrong with the grid setup!!" << std::endl;
	return 1;
   }


   //collect the name of dimensions: subchannels you want to vary; and the range
   std::vector<std::string> dimension_name;
   const double range_x_low = mygrid.f_dimensions.at(0).f_min;
   const double range_x_up = mygrid.f_dimensions.at(0).f_max;
   const double range_y_low = mygrid.f_dimensions.at(1).f_min;
   const double range_y_up = mygrid.f_dimensions.at(1).f_max;
   const double range_z_low = mygrid.f_dimensions.at(2).f_min;
   const double range_z_up = mygrid.f_dimensions.at(2).f_max;
   int nbin_x = mygrid.f_dimensions.at(0).f_N;  //number of point in x axis
   int nbin_y = mygrid.f_dimensions.at(1).f_N;  
   int nbin_z = mygrid.f_dimensions.at(2).f_N;  
 
   dimension_name.clear(); 
   for(int i=0; i< mygrid.f_num_dimensions ; i++){
	dimension_name.push_back(mygrid.f_dimensions.at(i).f_name);
   }



   //*********************loop over grid points, calc chi square********************************
   TFile* f_output = new TFile("chi_contour.root", "recreate");
   TH3D* h_chi2_raw = new TH3D("h_chi2_raw", Form("h_chi2_raw;%s;%s;%s", dimension_name[0].c_str(), dimension_name[1].c_str(), dimension_name[2].c_str()), nbin_x, range_x_low,range_x_up, nbin_y, range_y_low, range_y_up, nbin_z, range_z_low, range_z_up);
   TH3D* h_chi2_delta = new TH3D("h_chi2_delta", Form("h_chi2_raw;%s;%s;%s", dimension_name[0].c_str(), dimension_name[1].c_str(), dimension_name[2].c_str()), nbin_x, range_x_low,range_x_up, nbin_y, range_y_low, range_y_up, nbin_z, range_z_low, range_z_up);
   TH3D* h_chi2_inter=NULL;  


   std::cout << "Fraction Fit||"<<  "\tStart GLOBAL SCAN" <<std::endl;
   std::vector<double> chi;  //vector to save chi square values.
   chi.reserve(grid.size());  //reserve the memory

   for(int i =0; i< grid.size() ;i++){

	std::cout << "on Point " << i << std::endl;
	//set a temperary SBNspec, assume there is already a MC CV root file with corresponding sequence defined in xml	
	SBNspec spec_temp(mc_filename, xml, false);
	
	//access grid point
	std::vector<double> point = grid[i];

	//check dimension point
	if(point.size() != dimension_name.size()){
		std::cout << "dimension doesn't match: something wrong with grid setup!" << std::endl;
		return 1;
	}

	//scale chosen subchannels
	for(int j=0; j< point.size(); j++ ){
		spec_temp.Scale(dimension_name[j], point[j]);
	}

	spec_temp.CollapseVector();

	//initialize a SBNchi
	SBNchi chi_temp(xml);
	//inverted collapsed covariance matrix	
	TMatrixT<double> collapsed_temp(spec_temp.num_bins_total_compressed, spec_temp.num_bins_total_compressed);
	TMatrixT<double> invert_collapsed_temp(spec_temp.num_bins_total_compressed, spec_temp.num_bins_total_compressed);
	collapsed_temp = chi_temp.CalcCovarianceMatrixCNP(frac_syst_matrix, spec_temp.full_vector, spec_temp.collapsed_vector, data_spec.collapsed_vector);
	invert_collapsed_temp = chi_temp.InvertMatrix(collapsed_temp);
	//calculate chi2
	double chi_value = chi_temp.CalcChi(invert_collapsed_temp, spec_temp.collapsed_vector, data_spec.collapsed_vector);
	chi.push_back(chi_value);
	//chi.push_back(chi_temp.CalcChi(data_spec.collapsed_vector));
	h_chi2_raw->Fill(point[0], point[1], point[2],chi.back());

	//print out chi square value
	//std::cout << "CHI2: " << i << " " << chi.back();
	//for(int j=0; j<point.size() ; j++) std::cout << " " << point[j];
	//std::cout << std::endl;

   }
   
   f_output->cd();
   h_chi2_raw->Write();   

   double chi_min=DBL_MAX; // minimum of chi2.
   int chi_min_index = -99;   // index of the min chi value
   chi_min_index = std::min_element(chi.begin(), chi.end()) - chi.begin();
   chi_min = *std::min_element(chi.begin(), chi.end());

   //compare spectrum between data and the best fit point
   SBNspec spec_temp(mc_filename, xml, false);
   std::vector<double> best_point = grid[chi_min_index];
   for(int i=0; i< best_point.size(); i++ ){
                spec_temp.Scale(dimension_name[i], best_point[i]);
   }
   tag = "best_fit";
   spec_temp.CompareSBNspecs(&data_spec, tag);
   std::cout << "Fraction Fit||" <<  "\tBest Fit Point: (" << best_point[0] << ", " << best_point[1] << ", " << best_point[2] << ")"<<std::endl;


   std::cout << "Fraction Fit||\tPrint out delta chi square for grid points now" << std::endl;
   //delta chi2;
   for(int i=0; i< grid.size(); i++){
	 chi[i] -= chi_min;
	 std::vector<double> point = grid[i];
	 h_chi2_delta->Fill(point[0], point[1], point[2], chi[i]);	 
	 //std::cout << "DeltaChi: " << i << " " << chi[i];
	 //for(int j =0;j <point.size(); j++) std::cout << " " << point[j];
	 //std::cout << std::endl;
   }
   h_chi2_delta->SetMinimum(-1);  // to show the zero z values
   f_output->cd();
   h_chi2_delta->Write();
   std::cout << "Fraction Fit||" <<  "\t End of Global Scan--Chi2 calculation" << std::endl;

   std::cout << "Fraction Fit||" <<  "\t Start projecting on 2D plots" << std::endl;

   //histograms saving delta chi values with one parameter being marginalized
   TH2D* h_mchi2_xy = new TH2D("h_mchi2_xy", Form("h_mchi2_xy; %s;%s",dimension_name[0].c_str(), dimension_name[1].c_str()),nbin_x, range_x_low,range_x_up, nbin_y, range_y_low, range_y_up);
   TH2D* h_mchi2_xz = new TH2D("h_mchi2_xz", Form("h_mchi2_xz; %s;%s",dimension_name[0].c_str(), dimension_name[2].c_str()),nbin_x, range_x_low,range_x_up, nbin_z, range_z_low, range_z_up);
   TH2D* h_mchi2_yz = new TH2D("h_mchi2_yz", Form("h_mchi2_yz; %s;%s",dimension_name[1].c_str(), dimension_name[2].c_str()),nbin_y, range_y_low,range_y_up, nbin_z, range_z_low, range_z_up);

   //histograms saving delta chi values with respect to global minimum
   TH2D* h_gchi2_xy = new TH2D("h_gchi2_xy", Form("h_gchi2_xy; %s;%s",dimension_name[0].c_str(), dimension_name[1].c_str()),nbin_x, range_x_low,range_x_up, nbin_y, range_y_low, range_y_up);
   TH2D* h_gchi2_xz = new TH2D("h_gchi2_xz", Form("h_gchi2_xz; %s;%s",dimension_name[0].c_str(), dimension_name[2].c_str()),nbin_x, range_x_low,range_x_up, nbin_z, range_z_low, range_z_up);
   TH2D* h_gchi2_yz = new TH2D("h_gchi2_yz", Form("h_gchi2_yz; %s;%s",dimension_name[1].c_str(), dimension_name[2].c_str()),nbin_y, range_y_low,range_y_up, nbin_z, range_z_low, range_z_up);


   double max_bin=*std::max_element(chi.begin(), chi.end());;
   //set the bin content to DBL_MAX
   for(int ix=1;ix <= nbin_x; ix++){
	for(int iy=1; iy <= nbin_y; iy++) h_mchi2_xy->SetBinContent(ix, iy, max_bin);
	for(int iz=1; iz <= nbin_z; iz++) h_mchi2_xz->SetBinContent(ix, iz, max_bin);
   }
   for(int iy=1; iy <= nbin_y; iy++){
	for(int iz=1; iz <= nbin_z; iz++) h_mchi2_yz->SetBinContent(iy, iz, max_bin);
   }


   
   //marginalize one parameter
   for(int ix=0; ix < nbin_x; ix++){
	for(int iy=0; iy < nbin_y; iy++){
	   for(int iz=0 ; iz< nbin_z; iz++){
		int ip = ix*nbin_y*nbin_z + iy*nbin_z + iz; // index of grid point
		std::vector<double> point = grid[ip]; 


		//marginalized minimum
		//conditional operator, saver the smaller chi.
		if(chi[ip]< h_mchi2_xy->GetBinContent(ix+1, iy+1)){
			 h_mchi2_xy->SetBinContent(ix+1, iy+1, chi[ip]);
			std::cout << "chi2 value: " << chi[ip] << std::endl;
		}
		if(chi[ip]< h_mchi2_xz->GetBinContent(ix+1, iz+1)) h_mchi2_xz->SetBinContent(ix+1, iz+1, chi[ip]);
		if(chi[ip]< h_mchi2_yz->GetBinContent(iy+1, iz+1)) h_mchi2_yz->SetBinContent(iy+1, iz+1, chi[ip]);


		//global minimum
		if(point[2] == best_point[2]) h_gchi2_xy->Fill(point[0], point[1], chi[ip]);
		if(point[1] == best_point[1]) h_gchi2_xz->Fill(point[0], point[2], chi[ip]);
		if(point[0] == best_point[0]) h_gchi2_yz->Fill(point[1], point[2], chi[ip]);
	   }
	}
   } 

  f_output->cd();
  h_gchi2_xy->Write();
  h_gchi2_yz->Write();
  h_gchi2_xz->Write();
  h_mchi2_xy->Write();
  h_mchi2_yz->Write();
  h_mchi2_xz->Write();

   /*std::cout << "DeltaChi:" <<std::endl;
   for(int i =0; i< chi.size(); i++) std::cout << " " << chi[i];
   std::cout<< std::endl;
   */


   //****************** START INTERPOLATION*********************************

/*   if(interpolation_number != -99){

	//the upper range of x, y axis for interpolation
	double range_x_inter = mygrid.f_dimensions.at(0).f_max - mygrid.f_dimensions.at(0).f_step;
	double range_y_inter = mygrid.f_dimensions.at(1).f_max - mygrid.f_dimensions.at(1).f_step;
	std::cout << "range x interpolation " << range_x_inter << " Y: "<< range_y_inter <<std::endl;

	//save interpolated 2D plot
   	h_chi2_inter = new TH2D("h_chi2_interpolation", "h_chi2_interpolation", interpolation_number, range_x_low,range_x_inter, interpolation_number, range_y_low, range_y_inter);
	h_chi2_inter->SetMinimum(-1);
	//h_chi2_inter->SetAxisRange(range_x_low, range_x_inter, "X");
	//h_chi2_inter->SetAxisRange(range_y_low, range_y_inter, "Y");


	std::cout << "Fraction Fit||" << "\tStart interpolation with number "<< interpolation_number <<std::endl;
   	const gsl_interp2d_type *T = gsl_interp2d_bicubic;
	std::vector<double> grid_x;   //to save grid point values.
	std::vector<double> grid_y;
	grid_x = mygrid.f_dimensions.at(0).f_points;  //initialize grid point
	grid_y = mygrid.f_dimensions.at(1).f_points;   //initialize grid point
	//double x_step = fabs(range_x_low - range_x_inter)/(double)interpolation_number;
	//double y_step = fabs(range_y_low - range_y_inter)/(double)interpolation_number;
	double x_step = fabs(range_x_low - range_x_inter)/(interpolation_number - 1.);
	double y_step = fabs(range_y_low - range_y_inter)/(interpolation_number - 1.);
	std::cout << "Step size X:" << x_step << " Step size Y:" << y_step <<std::endl;

  	gsl_spline2d *spline = gsl_spline2d_alloc(T, nbin_y, nbin_x);  //due to the way grid points are sequenced in our code, we change the y and x position in gsl interpolation

  	gsl_interp_accel *xacc = gsl_interp_accel_alloc();
  	gsl_interp_accel *yacc = gsl_interp_accel_alloc();

	// initialize interpolation 
	gsl_spline2d_init(spline,  &grid_y[0], &grid_x[0], &chi[0] , nbin_y, nbin_x);


	// interpolate N values in x and y and print out grid for plotting 
	for (int i = 0; i < interpolation_number; i++){

	      double xi;

	      if(i == (interpolation_number -1)) xi = range_x_low + (i-0.5)*x_step;  //to fill the last bin
	      else xi = range_x_low + i*x_step;
	      //if(xi > range_x_inter) continue;  //skip points that are out of bounds

	      for (int j = 0; j < interpolation_number; j++)
		{

		  double yj;
	          if(j == (interpolation_number -1)) yj = range_y_low + (j-0.5)*y_step;
	          else yj = range_y_low + j*y_step;
		  //if(yj > range_y_inter) continue;   //skip points that are out of bounds

		  double zij = gsl_spline2d_eval(spline, yj, xi, yacc, xacc);  //again, y and x values are exchanged here

		  h_chi2_inter->Fill(xi, yj, zij);

		  std::cout << "Interpolated value: " << zij << " ||Grid Point:(" << xi << ", " << yj << ")" << std::endl; 
		}
	}

	  f_output->cd();
	  h_chi2_inter->Write();

	  gsl_spline2d_free(spline);
	  gsl_interp_accel_free(xacc);
	  gsl_interp_accel_free(yacc);
  }
  else{ 
	h_chi2_inter = (TH3D*)h_chi2_delta->Clone();
  }

   //*****************END OF INTERPOLATION*************************************************






  //****************draw contour, get range for fraction fit****************************

   std::cout << "Fraction Fit||" <<  "\tDraw contours"<<std::endl;
   // empty TH2D to draw x & y axis on canvas
   TH2D* hr = NULL;
   if(interpolation_number != -99)
	hr = new TH2D("hr", Form("contour;%s;%s", dimension_name[0].c_str(), dimension_name[1].c_str()), interpolation_number, range_x_low,range_x_up, interpolation_number, range_y_low, range_y_up);
   else hr = new TH2D("hr", Form("contour;%s;%s", dimension_name[0].c_str(), dimension_name[1].c_str()), nbin_x, range_x_low,range_x_up, nbin_y, range_y_low, range_y_up);
   hr->SetStats(kFALSE);

   TGraph* c_graph=NULL;
   TGraph* temp_graph = NULL;
   TList* con_list=NULL;
   TMarker* marker = new TMarker(best_point[0], best_point[1], 29);  // a marker at best fit point
   marker->SetMarkerSize(2);

   std::vector<double> contour{2.30, 4.61,6.18,11.83}; // chi2 value for 2 dof with 1 sigma, 90%, 2 sigma and 3 sigma confidence level
   std::vector<std::string> CL_string{"1sigma", "90", "2sigma", "3sigma"};
   TCanvas* c = new TCanvas("c", "c");
   h_chi2_inter->SetContour((int)contour.size(), &contour[0]);
   h_chi2_inter->Draw("CONT Z LIST");//"LIST" generates a list of TGraph for each contour
   c->Update();
   f_output->cd();
   c->Write("contour");

   //grab contour object
   TObjArray *conts = (TObjArray*)gROOT->GetListOfSpecials()->FindObject("contours");
   TCanvas* v_canvas[conts->GetSize()];
	
   //hr->Draw();

   //loop over contours
   for(int i=0; i<conts->GetSize(); i++){
	con_list = (TList*)conts->At(i);   //a list of TGraph for i'th contour.
	std::cout << "Contour Z=" << contour[i] << " has "<< con_list->GetSize() << " enclosed region(s)! " << std::endl;

	double x_min = DBL_MAX;
	double x_max = DBL_MIN;
	double y_min = DBL_MAX;
	double y_max = DBL_MIN;

	v_canvas[i] = new TCanvas(Form("c_%s", CL_string[i].c_str()), Form("c_%s", CL_string[i].c_str()));
	hr->Draw();
	marker->Draw();
	c_graph = (TGraph*)con_list->First();  //grab the TGraph
	

	for(int j=0; j< con_list->GetSize() ; j++){
		temp_graph= (TGraph*)c_graph->Clone();
                temp_graph->Draw("C");
		//std::cout << "TGraph "<< j << ": number of points " << copy_graph->GetN() << std::endl;
			
		// grab x ,y coordinate
		double x,y;
		for(int k =0; k< temp_graph->GetN(); k++){
			temp_graph->GetPoint(k, x, y);
			if(x < x_min) x_min = x;
			if(x > x_max) x_max = x;

			if(y < y_min) y_min =y;
			if(y > y_max) y_max = y;
		}
	
		//move to next graph
		c_graph=(TGraph*)con_list->After(c_graph);
	}

	v_canvas[i]->Update();
	f_output->cd();
	v_canvas[i]->Write();

	std::cout << "Contour " << CL_string[i] << ": " << contour[i] << std::endl;
	std::cout << "range for x :" << x_min << "~" << x_max << std::endl;
	std::cout << "range for y : " << y_min << "~" << y_max << std::endl;
   } // contour loop
*/
    f_output->Close(); 

    std::cout << "Fraction Fit||" << "\tFinished" <<std::endl;
    std::cout << "Total wall time: " << difftime(time(0), start_time)/60.0 << " Minutes.\n";
    return 0;

}

