#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <vector>
#include <unistd.h>
#include <getopt.h>
#include <cstring>
#include <algorithm>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TString.h"
#include "TNtuple.h"
#include "TLine.h"
#include "TChain.h"
#include "TMath.h"
#include "TSystem.h"
#include "TMatrixT.h"
#include "TRandom.h"
#include "TError.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TMultiGraph.h"

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

        

static int compare (const void * a, const void * b)
{// Comparator function for getting the critical value
      if (*(double*)a > *(double*)b) return 1;
        else if (*(double*)a < *(double*)b) return -1;
          else return 0;  
}

double getCritValue(TTree *t, std::string name, double pval){
    //Gets the critical value of DeltaChi^2 such than pval % lies below. 
    //Sort the array of all psuedo-universes and then the required value can be found at ceil(Nentries*pval)

    int Nentries = t->GetEntries(); 
    const unsigned vals = t->Draw(name.c_str(), "1");
    if (vals) {
       double* x = t->GetV1();
       qsort(x, t->GetEntries(), sizeof(double), compare);
   
       int pos = ceil(Nentries*pval);  
       //std::cout<<"getCritValue at pos "<<pos<<" val "<<x[pos]<<" (check: )"<< t->GetEntries(("delta_chi2<"+std::to_string(x[pos])).c_str())/(double)t->GetEntries()<<std::endl;

       return x[pos];
    }
    
    return -99;
}



double Median(const TH1D * h1) { 

    int n = h1->GetXaxis()->GetNbins();  
    std::vector<double>  x(n);
    h1->GetXaxis()->GetCenter( &x[0] );
    const double * y = h1->GetArray(); 
    // exclude underflow/overflows from bin content array yG
    return TMath::Median(n, &x[0], &y[1]); 
}


double quick_median(std::vector<double> &v)
{
    size_t n = v.size() / 2;
    std::nth_element(v.begin(), v.begin()+n, v.end());
    return v[n];
}

double lin_interp(double x0, double x1, double y0, double y1, double x){
    return (y0*(x1-x)+y1*(x-x0))/(x1-x0);
}
/*************************************************************
 *************************************************************
 *		BEGIN sbnfit_make_covariance.cxx
 ************************************************************
 ************************************************************/
void runHelp(){
    std::cout<<"---------------------------------------------------"<<std::endl;
    std::cout<<"Modified single subchannel scaling feldman_cousins confidence belt constructor"<<std::endl;
    std::cout<<"---------------------------------------------------"<<std::endl;
    std::cout<<"--- Required arguments: ---"<<std::endl;
    std::cout<<"\t-x\t--xml\t\tInput configuration .xml file for SBNconfig"<<std::endl;
    std::cout<<"\t-t\t--tag\t\tA unique tag to identify the inputs/outputs [Default to TEST]"<<std::endl;
    std::cout<<"\t-i\t--input\t\tInput subchannel to scale (no default, required argument)"<<std::endl;
    std::cout<<"\t-g\t--grid\t\tGrid to scan, in the form 'min max num_steps' (default '1e-4 10.0 20')"<<std::endl;
    std::cout<<"\t-m\t--mode\t\tWhat mode you want to run in. Arguments are:"<<std::endl;
    std::cout<<"\t\t\t--\t feldman : Perform the pseudo universe grid scan (run first)"<<std::endl;  
    std::cout<<"\t\t\t--\t belt: Constructs the confidence belts, must be run after'feldman'"<<std::endl;
    std::cout<<"\t\t\t--\t data: Pass in an optional datafile, that will be compared to the grid"<<std::endl;
    std::cout<<"--- Optional arguments: ---"<<std::endl;
    std::cout<<"\t-s\t--stat\t\tStatistical error only mode, will ignore any covariance matrix passed in"<<std::endl;
    std::cout<<"\t-n\t--number\t\tNumber of pseudo-experiments to simulate (default 2500)"<<std::endl; 
    std::cout<<"\t-r\t--randomseed\t\tRandomNumber Seed (default from machine)"<<std::endl; 
    std::cout<<"\t-c\t--cnp\t\tuse a Combined Newman Pearson chi2 (default false)"<<std::endl;
    std::cout<<"\t-d\t--data\t\ta data SBNspec file to input, use with mode data"<<std::endl;
    std::cout<<"\t-h\t--help\t\tThis help menu."<<std::endl;
    std::cout<<"---------------------------------------------------"<<std::endl;
    return;
}

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
        {"stat", 		no_argument, 		0, 's'},
        {"number", 		required_argument,	0,'n'},
        {"cnp", 		no_argument,	0,'c'},
        {"grid", 		required_argument,	0,'g'},
        {"tag", 		required_argument,	0, 't'},
        {"mode",        required_argument, 0 ,'m'},
        {"data",        required_argument, 0 ,'d'},
        {"input",       required_argument, 0 ,'i'},
        {"randomseed",        required_argument, 0 ,'r'},
        {"help", 		no_argument,	0, 'h'},
        {0,			    no_argument, 		0,  0},
    };

    int iarg = 0;
    opterr=1;
    int index;

    //a tag to identify outputs and this specific run. defaults to EXAMPLE1
    std::string tag = "TEST";
    std::string mode_option;
    bool bool_stat_only = false;
    int number = 2500;
    double random_number_seed = -1;
    bool use_cnp = false;

    std::string grid_string = "1e-4 8.0 33";
    std::string input_scale_subchannel = "unset";
    std::string data_file_input = "null";

    while(iarg != -1)
    {
        iarg = getopt_long(argc,argv, "x:t:m:n:r:d:p:i:g:sch", longopts, &index);

        switch(iarg)
        {
            case 'x':
                xml = optarg;
                break;
            case 'n':
                number = (int)strtod(optarg,NULL);
                break;
            case 't':
                tag = optarg;
                break;
            case 'g':
                grid_string = optarg; 
                break;
            case 'd':
                data_file_input = optarg; 
                break;

            case 'i':
                input_scale_subchannel = optarg;
                break;
            case 'm':
                mode_option = optarg;
                break;
            case 'c':
                use_cnp = true;
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
                runHelp(); 
                return 0;
        }
    }



    /*************************************************************
     *************************************************************
     *			Main Program Flow
     ************************************************************
     ************************************************************/
    time_t start_time = time(0);
    int bf_data_pt = -1;
    std::cout<<"Begining SBNfit uboone subchannel scaling Feldman Cousins confidence belt constructor for tag: "<<tag<<std::endl;

    if(input_scale_subchannel=="unset"){
        std::cout<<"Error! you must set a value for which input subchannel to scale, e.g using --input/ -i 'nu_uBooNE_1g1p_ncdelta'"<<std::endl;
        std::cout<<"Please see...."<<std::endl;
        runHelp();
        return 0;
    }

    NGrid mygrid;
    mygrid.AddDimension(input_scale_subchannel,grid_string);
   
    mygrid.Print();
    SBNfeld myfeld(mygrid,tag,xml);

    if(mode_option == "feldman"){

        std::cout<<"Begininning a full Feldman-Cousins analysis for tag : "<<tag<<std::endl;

        myfeld.SetFractionalCovarianceMatrix(tag+".SBNcovar.root","frac_covariance");
        //myfeld.SetFractionalCovarianceMatrix(Msys);
        //myfeld.SetStatOnly();
        myfeld.m_subchannel_to_scale = input_scale_subchannel;

        myfeld.SetCoreSpectrum(tag+"_CV.SBNspec.root");
        myfeld.SetBackgroundSpectrum(tag+"_CV.SBNspec.root",input_scale_subchannel,0.0);
        myfeld.GenerateScaledSpectra();

        std::cout<<"Setting random seed "<<random_number_seed<<std::endl;
        myfeld.SetRandomSeed(random_number_seed);
        myfeld.SetNumUniverses(number);

        std::cout<<"Calculating the necessary SBNchi objects"<<std::endl;
        myfeld.CalcSBNchis();
        std::cout <<"DONE calculating the necessary SBNchi objects at : " << difftime(time(0), start_time)/60.0 << " Minutes.\n";

        std::cout<<"Beginning to peform FullFeldmanCousins analysis"<<std::endl;
        myfeld.FullFeldmanCousins();

    }else if(mode_option=="data"){

        std::cout<<"Begininning a real data analysis for tag : "<<tag<<std::endl;

        myfeld.SetFractionalCovarianceMatrix(tag+".SBNcovar.root","frac_covariance");
        myfeld.m_subchannel_to_scale = input_scale_subchannel;

        if(use_cnp) myfeld.UseCNP();

        myfeld.SetCoreSpectrum(tag+"_CV.SBNspec.root");
        myfeld.SetBackgroundSpectrum(tag+"_CV.SBNspec.root",input_scale_subchannel,1.0);
        myfeld.GenerateScaledSpectra();

        std::cout<<"Setting random seed "<<random_number_seed<<std::endl;
        myfeld.SetRandomSeed(random_number_seed);
        myfeld.SetNumUniverses(number);

        std::cout<<"Calculating the necessary SBNchi objects"<<std::endl;
        myfeld.CalcSBNchis();
        std::cout <<"DONE calculating the necessary SBNchi objects at : " << difftime(time(0), start_time)/60.0 << " Minutes.\n";

        SBNspec * datain = new SBNspec(data_file_input.c_str(),xml);
        //myfeld.CompareToData(datain,{2.48,1.52,1.24},{5.68,6.64,7.36}); // 21.3
        //myfeld.CompareToData(datain,{1.84,1.08,0.88},{4.84,5.92,6.52}); //22.4
        //myfeld.CompareToData(datain,{0,0,0},{1.36,2.16,2.64}); //23.1
        //myfeld.CompareToData(datain,{1.2,0.68,0.346668},{4.08,5.12,5.6}); //24.1
        //myfeld.CompareToData(datain,{0.346668,0.0,0.0},{2.64,3.52,4.08}); //25.1

//        myfeld.CompareToData(datain,{0.320001,0.0,0.0},{2.56,3.52,4.0}); //25.1 v5
  
        
//        myfeld.CompareToData(datain,{0.0,-2.0,-2.0},{2.8,3.6,4.2}); //25.1 v6
        myfeld.CompareToData(datain); //25.1 v6


    }else if(mode_option == "belt"){

        if(bool_stat_only){
            myfeld.SetEmptyFractionalCovarianceMatrix();
            myfeld.SetStatOnly();
            std::cout<<"RUNNING Statistics uncertainty only!"<<std::endl;
        }else{
            myfeld.SetFractionalCovarianceMatrix(tag+".SBNcovar.root","frac_covariance");
        }
        myfeld.m_subchannel_to_scale = input_scale_subchannel;

        if(use_cnp) myfeld.UseCNP();

        myfeld.SetCoreSpectrum(tag+"_CV.SBNspec.root");
        myfeld.SetBackgroundSpectrum(tag+"_CV.SBNspec.root",input_scale_subchannel,1.0);
        myfeld.GenerateScaledSpectra();

        std::cout<<"Calculating the necessary SBNchi objects"<<std::endl;
        myfeld.CalcSBNchis();
        std::cout <<"DONE calculating the necessary SBNchi objects at : " << difftime(time(0), start_time)/60.0 << " Minutes.\n";

        TFile *f = new TFile("scan.root","recreate");
        f->cd();
        std::cout<<"Starting to peform a globalScan analysis"<<std::endl;
        std::vector<std::vector<double>> vec_grid = mygrid.GetGrid();

        TFile *fin = new TFile(("SBNfeld_output_"+tag+".root").c_str(),"read");

        //Some Manual Color Changing and such
        std::vector<double> plotting_pvals = {0.68, 0.90, 0.95, 0.97};
        //std::vector<double> plotting_pvals = {0.68, 0.90, 0.95};
        std::vector<std::string> plotting_strs = {"68%","90%","95%","97%"};
       // std::vector<std::string> plotting_strs = {"68%","90%","95%"};
        std::vector<int> gcols = {kGreen+3,kGreen+2,kGreen-3,kGreen-9};
        //std::vector<int> gcols = {kRed-9,kBlue-9,kGreen-9};
        
        
        std::vector<double> v_median;
        std::vector<double> v_true;
        std::vector<double> v_1sigma_p;
        std::vector<double> v_1sigma_m;
        std::vector<std::vector<double>> v_min;v_min.resize(plotting_pvals.size());
        std::vector<std::vector<double>> v_max;v_max.resize(plotting_pvals.size());

        std::cout<<"MPrinting stuff"<<std::endl;
        TH2D * f_FC = new TH2D("f_FC","f_FC",vec_grid.size(),vec_grid.front()[0],vec_grid.back()[0],vec_grid.size(),vec_grid.front()[0],vec_grid.back()[0]);

        double bfval_v[vec_grid.size()];

       //TESTING 
        std::vector<double> v_gridvals;
        for(size_t t =0; t < vec_grid.size(); t++){
            v_gridvals.push_back(vec_grid[t][0]);
        }

        for(int i=0; i< vec_grid.size(); i++){

            v_true.push_back(vec_grid[i][0]);

            //Whats the critical value? Need to get the information stored in root files after the initial calc mode.
            TTree *t =  (TTree*)fin->Get(("ttree_"+std::to_string(i)).c_str());
            TH1D * cumul = (TH1D*)fin->Get(("delta_chi2_"+std::to_string(i)+"_cumulative").c_str());
            TH1D * h_bfval = (TH1D*)fin->Get(("bf_value_"+std::to_string(i)).c_str());//Added by Ivan

            //Get the likelihood that we want to compare to. This is currently a bit slow. Note perform Iterative Fit, no longer performs an iterative fit, but an exact grid minimization
            std::vector<double> comparative_likelihood;
            TMatrixT<double> inverse_background_collapsed_covariance_matrix(myfeld.num_bins_total_compressed,myfeld.num_bins_total_compressed);   
            for(int j=0; j< vec_grid.size(); j++){
                std::vector<double> this_likelihood_vector = myfeld.PerformIterativeGridFit(myfeld.m_cv_spec_grid.at(j)->f_collapsed_vector, i, inverse_background_collapsed_covariance_matrix);
                comparative_likelihood.push_back(this_likelihood_vector[1]);
            }

            //The liklihood for a true scale_i, given asimov dataset for all other values of scalefactors 
            TGraph *g_likelihood = new TGraph(comparative_likelihood.size(),&v_gridvals[0],&comparative_likelihood[0]);

            //calculate the acceptance regions for all confidence regions expected
            for(int p =0; p< plotting_pvals.size(); ++p){
                double plotting_pval = plotting_pvals[p];
                //First lets find a critical chi^2 for this confidence level. This is the critical value that plotting_pval % of events simulated have a delta chi^2 < Critical Value
                double critical_delta_chi2 = getCritValue(t,"delta_chi2", plotting_pval);

                //initilize the limits at grid boundaries.
                double low_limit = v_true.front();
                double up_limit = vec_grid[vec_grid.size()-1][0];
                
                //Find lower limit, starting from the grid point, walk backwards until the liklihood (or deltachi) is above the critical value.
                for(int low = i; low>=0;low--){
                    double val = g_likelihood->Eval(vec_grid[low][0]); 
                    if(val >= critical_delta_chi2){
                        //If the value isn't EXACTLY the critical chi, move to grid point below. Thus we can slighly overcover, never undercover. 
                        if(val!= critical_delta_chi2){low=low-1;}
                        if(low<0)low=0; //if BF is 0, and overestimates
                        low_limit = vec_grid[low][0];
                        break;
                    }
                }

                v_min[p].push_back(low_limit);
                v_max[p].push_back(up_limit);

                std::cout<<"Grid point "<<i<<" (mu= "<<vec_grid[i][0] <<") has a critical delta chi of "<<critical_delta_chi2<<" for a pval of "<<plotting_pval<<" whose grid acceptance region is ["<<low_limit<<" , "<<up_limit<<"]"<<std::endl; 

            }//end pval loop



            std::vector<double> whatsmedian(t->GetEntries(),-9);
            double f_bfval = 0;
            double whatsmean = 0;
            t->SetBranchAddress("bf_gridvalue",&f_bfval);
            for(int k=0; k<t->GetEntries(); k++){
                t->GetEntry(k);               
                whatsmedian[k]=f_bfval;
                whatsmean+=f_bfval;
            }
            whatsmean=whatsmean/(double)t->GetEntries();


            //Get median of BF value
            h_bfval->ComputeIntegral();
            std::vector<double> pvalues = { 0.5 };
            std::vector<double> bf_val_quantiles(pvalues.size());
            h_bfval->GetQuantiles(pvalues.size(),&bf_val_quantiles[0], &pvalues[0]);
            double bfval = bf_val_quantiles[0];
            bfval_v[i] = bfval;
            bfval_v[i] = quick_median(whatsmedian);
            
            //bfval_v[i] = Median(h_bfval);
            //bfval_v[i] = whatsmean;

            delete cumul;
        }


        TCanvas *c3 = new TCanvas("h_ono_r3");
        c3->SetFillStyle(0);

        TPad *pad = new TPad("pad", "pad", 0, 0, 0.8, 1.0);
        pad->SetRightMargin(0); // Upper and lower plot are joined
        pad->Draw();             // Draw the upper pad: pad
        pad->cd();               // pad becomes the current pad

        std::vector<TGraph*> gmaxs;
        std::vector<TGraph*> gmins;
        std::vector<TGraph*> grshades;

        TLegend * l_probs = new TLegend(0.11,0.52,0.89,0.89);//69 was 29

        TMultiGraph *mg = new TMultiGraph();
        mg->SetTitle("Feldman Cousins Corrected Confidence Belt");

        for(int p=plotting_pvals.size()-1; p>=0;--p){
            pad->cd();

            gmins.push_back(new TGraph(v_true.size(),&(v_min[p])[0], &v_true[0]));
            gmaxs.push_back(new TGraph(v_true.size(),&(v_max[p])[0], &v_true[0]));

            grshades.push_back( new TGraph(2*v_true.size()));

            for (int i=0;i<v_true.size();i++) {
                int n = v_true.size();
                grshades.back()->SetPoint(i,v_min[p][i],v_true[i]);
                grshades.back()->SetPoint(n+i,v_max[p][n-i-1],v_true[n-i-1]);
            }
            grshades.back()->SetFillColor(gcols[p]);
            gmins.back()->SetLineWidth(2);
            gmins.back()->SetLineColor(kBlack);
            gmaxs.back()->SetLineWidth(2);
            gmaxs.back()->SetLineColor(kBlack);
            
            l_probs->AddEntry(grshades.back(), plotting_strs[p].c_str() ,"f");
        }
        
        //l_probs->SetHeader("#splitline{#splitline{Classical}{Confidence}}{#splitline{Level of}{Interval}}");
        l_probs->SetHeader("#splitline{Confidence}{#splitline{Level of}{Interval}}");

        for(auto &g:grshades)mg->Add(g);
        pad->cd();
        mg->Draw("ALF");

        mg->GetXaxis()->SetTitle("Measured #Delta Radiative Rate (#hat{x}_{#Delta})");
        mg->GetYaxis()->SetTitle("True #Delta Radiative Rate (x_{#Delta})");


        double mplot = v_true.back();
        double ymplot  = std::max(0.0, v_true.front());

        mg->SetMinimum(ymplot);

        mg->GetXaxis()->SetLimits(v_true.front(),mplot);      
        mg->GetHistogram()->SetMaximum(mplot);//v_true.back());          
        mg->GetHistogram()->SetMinimum(ymplot);     
 
        for(auto &g:gmins)g->Draw("l same");
        for(auto &g:gmaxs)g->Draw("l same");

        //also draw median
        double npoints = vec_grid.size();
        double maxpt = vec_grid.back()[0];
        double scale = (npoints-1)/maxpt;

        mg->GetXaxis()->SetLimits(v_true.front(),maxpt);
        mg->GetHistogram()->SetMaximum(maxpt);//v_true.back());          

        TLine lcross(ymplot,ymplot, mplot, mplot);
            lcross.SetLineStyle(9);
            lcross.SetLineWidth(1);
            lcross.SetLineColor(kBlack);
            lcross.Draw("same");

        TLine lv1(1.0,0.0,1.0,maxpt);
        lv1.SetLineStyle(2);
        lv1.SetLineWidth(1);
        lv1.SetLineColor(kBlack);
        //lv1.Draw("same");

        TLine lh3(0.0,3.1, 1.0,3.1);
        lh3.SetLineStyle(2);
        lh3.SetLineWidth(1);
        lh3.SetLineColor(kBlack);
        //lh3.Draw("same");


        double u_measured[vec_grid.size()];
        double uexp = 0;

        std::cout << "//////////HERE!!!!!!!! vec_grid.size() = " << vec_grid.size() << std::endl;
        for(int k = 0; k < npoints; k++){
            uexp = ((double) k)/scale;
            std::cout << "//////////HERE!!!!!!!! k, uexp, median = " << k << ", " << uexp << ", " << bfval_v[k] << std::endl;
            u_measured[k] = uexp;
        }

        TGraph *bf_vals_g = new TGraph(vec_grid.size(), bfval_v, &v_true[0]);
        //bf_vals_g->Draw("same *");

        pad->Update();
        pad->RedrawAxis();
        // TLine l;
        //  l.DrawLine(gPad->GetUxmin(), gPad->GetUymax(), gPad->GetUxmax(), gPad->GetUymax());
        //   l.DrawLine(gPad->GetUxmax(), gPad->GetUymin(), gPad->GetUxmax(), gPad->GetUymax());

        c3->cd();
        TPad *padl = new TPad("padl", "padl", 0.8, 0, 1, 1);
        padl->SetBottomMargin(0.2);
        padl->Draw();
        padl->cd();       // padl becomes the current pad
        l_probs->SetTextSize(0.0775);
        l_probs->Draw();
        l_probs->SetLineColor(kWhite);
        l_probs->SetLineWidth(0);



        c3->SaveAs(("FC_confidence_belt_"+tag+".pdf").c_str(),"pdf");

        std::cout<<"**************** Feldman Cousins 1D Confidence Intervals  **********************"<<std::endl;
        for(int i=0; i<v_true.size(); i++){
            std::cout<<"Grid Pt: "<<i<<", ScaleFactor: "<<vec_grid[i][0]<<std::endl;
            std::cout<<"- Median: "<<bfval_v[i]<<std::endl;
            for(int p=0; p< plotting_pvals.size();++p){

                double measured_val = vec_grid[i][0];
                //std::vector<double> reg = myfeld.getConfidenceRegion(gmins[plotting_pvals.size()-p-1],gmaxs[plotting_pvals.size()-p-1],measured_val);
                std::vector<double> reg = myfeld.getConfidenceRegion(v_min[p][i],v_max[p][i],measured_val);
                std::cout<<"-- CL: "<<plotting_pvals[p]<<"  Sigma: "<<sqrt(2)*TMath::ErfInverse(plotting_pvals[p])<<"  ConfidenceInterval: ["<<reg[0]<<" -> "<<reg[1]<<"]"<<std::endl;
            }
        }
    }else{
        std::cout<<"The mode you asked for ("<<mode_option<<") is not available, the available options are.."<<std::endl;
        runHelp();

    }
    std::cout << "Fin. Total wall time: " << difftime(time(0), start_time)/60.0 << " Minutes.\n";

    return 0;

}


