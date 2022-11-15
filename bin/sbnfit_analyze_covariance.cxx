#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <vector>
#include <unistd.h>
#include <getopt.h>
#include <cstring>

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

#include "params.h"
#include "SBNconfig.h"
#include "SBNchi.h"
#include "SBNspec.h"
#include "SBNosc.h"
#include "SBNfit.h"
#include "SBNfit3pN.h"
#include "SBNcovariance.h"

#define no_argument 0
#define required_argument 1
#define optional_argument 2

using namespace sbn;

/*************************************************************
 *************************************************************
 *		BEGIN sbnfit_analyze_covariance.cxx
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
        {"help", 		no_argument,	0, 'h'},
        {"covar",		required_argument,    0, 'c'},
        {"flat", required_argument,0,'f'},
        {"zero",no_argument,0,'z'},
        {"cmin",required_argument,0,'k'},
        {"cmax",required_argument,0,'p'},
        {0,			    no_argument, 		0,  0},
    };

    int iarg = 0;
    opterr=1;
    int index;
    std::vector<std::string> covar_files;
    bool stats_only = true;
    std::string signal_file;

    bool bool_flat_det_sys = false;
    double flat_det_sys_percent = 0.0;

    bool remove_correlations = false;

    double cmin = 0;
    double cmax = -9;

    //a tag to identify outputs and this specific run. defaults to EXAMPLE1
    std::string tag = "TEST";

    while(iarg != -1)
    {
        iarg = getopt_long(argc,argv, "f:x:s:t:c:k:p:zh", longopts, &index);

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
            case 'c':
                covar_files.push_back(optarg);
                for (int i = optind; i < argc; i++) {
                    covar_files.push_back(argv[i]);
                }
                break;

            case 't':
                tag = optarg;
                break;
            case 's':
                signal_file = optarg;
                break;

            case '?':
            case 'h':
                std::cout<<"---------------------------------------------------"<<std::endl;
                std::cout<<"sbnfit_anayze_covariance allows for the analysis of covariance matricies from input root files containing reconstructed variables and covariance matricies. "<<std::endl;
                std::cout<<"---------------------------------------------------"<<std::endl;
                std::cout<<"--- Required arguments: ---"<<std::endl;
                std::cout<<"\t-x\t--xml\t\tInput configuration .xml file for SBNconfig"<<std::endl;
                std::cout<<"\t-t\t--tag\t\tA unique tag to identify the outputs [Default to TEST]"<<std::endl;
                std::cout<<"\t-s\t--signal\t\tInput signal SBNspec.root file"<<std::endl;
                std::cout<<"\t-c\t--covar\t\tList of covariances to plot"<<std::endl;
                std::cout<<"--- Optional arguments: ---"<<std::endl;
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
    std::cout<<"Loading SBNspec file : "<<signal_file<<" with xml "<<xml<<std::endl;
    SBNspec sig(signal_file,xml);
    SBNspec sig_NoMC(signal_file,xml);
    sig.CalcFullVector();
    sig_NoMC.RemoveMCError();
    sig_NoMC.CalcFullVector();

    std::cout<<"Loading fractional covariance matrix from "<<covar_files.size()<<std::endl;

    std::cout<<"The Covar File string is length: "<<covar_files.size()<<std::endl;
    for(auto &f: covar_files) std::cout<<" "<<f<<std::endl;
    std::vector<TFile*> files;

    std::vector<TMatrixD> fmats;
    for(auto s: covar_files){
        std::cout<<"Loading Covar File String: "<<s.c_str()<<std::endl;
        files.push_back(new TFile(s.c_str(),"read") );
        TMatrixD fm = *(TMatrixD*)files.back()->Get("frac_covariance");
        std::cout<<s<<" has Dimensions: "<<fm.GetNrows()<<" "<<fm.GetNcols()<<std::endl;
        fmats.push_back(fm);
    }

    std::vector<double> summed_up(sig.num_bins_total_compressed,0.0);
    std::vector<std::vector<double>> v_sys;
    std::vector<double> v_stat;
    std::vector<double> v_mcstat;

    std::cout<<"Adding all together"<<std::endl;
    TMatrixD m_fsum = fmats[0];
    std::cout<<fmats.size()<<std::endl;


    for(int i=1; i< fmats.size(); i++){
        m_fsum = m_fsum + fmats[i];
        std::cout<<"Now adding : "<<i<<std::endl;
    }

    std::vector<SBNchi*> chis;
    SBNchi AllChi(sig, m_fsum);

    for(int i=0; i< fmats.size(); i++){
        SBNchi * chi = new SBNchi(sig_NoMC, (fmats[i]));
        chis.push_back(chi);
    }

    for(int i=0; i< fmats.size(); i++){
        std::cout<<"On Covar : "<<i<<std::endl;
        std::vector<double> tmp;
        for(int j=0; j<chis[i]->vec_matrix_collapsed.size(); j++){
               double covar = (chis[i]->vec_matrix_collapsed.at(j).at(j)-chis[i]->core_spectrum.collapsed_vector.at(j));
               double fcovar = covar/pow(chis[i]->core_spectrum.collapsed_vector.at(j),2);
               std::cout<<sqrt(fcovar)<<" ";
               summed_up[j]+=fcovar;
               tmp.push_back(sqrt(fcovar));
        }
        v_sys.push_back(tmp);
        std::cout<<std::endl;
    }

    sig.CollapseVector();
    std::cout<<"On ``Stats`` : "<<std::endl;
    for(int i=0; i< sig.collapsed_vector.size(); i++){
        double err = sqrt(sig.collapsed_vector[i]);
        double ferr = err/sig.collapsed_vector[i];
        std::cout<<ferr<<" ";
        summed_up[i]+= pow(ferr,2);
        v_stat.push_back(ferr);
    }
    std::cout<<std::endl;

    std::cout<<"On ``Intrinsic MC Stats`` : "<<std::endl;
    TMatrixD mcstats(sig.num_bins_total,sig.num_bins_total);
    TMatrixD mcstats_collapsed(sig.num_bins_total_compressed,sig.num_bins_total_compressed);
    mcstats.Zero();
    for(int i=0; i<sig.full_err_vector.size();i++){
        mcstats(i,i) = pow(sig.full_err_vector[i],2);
    }
    AllChi.CollapseModes(mcstats,mcstats_collapsed);
    for(int i=0; i<mcstats_collapsed.GetNcols();i++){
        double ferr = mcstats_collapsed(i,i)/pow(sig.collapsed_vector.at(i),2); 
        std::cout<<sqrt(ferr)<<" ";
        summed_up[i]+=ferr;
        v_mcstat.push_back(sqrt(ferr));
    }
    std::cout<<std::endl;

    std::cout<<"--------------- TOTAL Sys ---------------"<<std::endl;
    TMatrixD Allsys;
    AllChi.FillCollapsedFractionalMatrix(&Allsys);
    std::cout<<"Total Sys : "<<std::endl;
    for(int j=0; j<Allsys.GetNcols(); j++) std::cout<<sqrt(Allsys(j,j))<<" ";
    std::cout<<std::endl;

    std::cout<<"--------------- TOTAL Sys + stat ---------------"<<std::endl;
    std::cout<<"Total Sys +stat : "<<std::endl;
    for(int j=0; j<Allsys.GetNcols(); j++)std::cout<<sqrt(Allsys(j,j)+sig.collapsed_vector[j]/pow(sig.collapsed_vector[j],2))<<" ";
    std::cout<<std::endl;

    std::cout<<"--------------- Summed ---------------"<<std::endl;
    std::cout<<"SummedUp : "<<std::endl;
    for(int j=0; j<summed_up.size(); j++)std::cout<<sqrt(summed_up[j])<<" ";
    std::cout<<std::endl;

    std::cout<<"--------------- Summed ---------------"<<std::endl;
    //OK, summed_up, v_sys, v_stat, v_mcstat
    auto mapo = sig.GetCollapsedChannelIndicies();
    auto vhist = sig.GetBlankChannelHists();

    auto v_stat_hist = vhist;
    auto v_mcstat_hist = vhist;
    auto v_summed_hist = vhist;
    std::vector<std::vector<TH1D>> v_sys_hist(v_sys.size(),vhist);

    std::vector<int> cols = {kRed-7,kBlue-7,kGreen-3,kCyan};
    std::vector<std::string> nams = {"Detector Systematics","Flux Systematics","GENIE Systematics","Geant4 Systematics"};

    gStyle->SetOptStat(0);

    for(int i=0; i<vhist.size();i++){
            TCanvas *c = new TCanvas(std::to_string(i).c_str(),std::to_string(i).c_str(),1100,1000);
            c->cd();
            
            auto vec = mapo[i];
            std::cout<<"Chan "<<i<<" : from "<<vec[0]<<" to "<<vec[1]<<std::endl;

            int bincount = 1;
            v_stat_hist[i].Reset();
            v_mcstat_hist[i].Reset();
            v_summed_hist[i].Reset();
            for(int k=0; k<v_sys.size();k++){
                     v_sys_hist[k][i].Reset();
            }
            for(int j=vec[0]; j<vec[1]; j++){
                v_stat_hist[i].SetBinContent(bincount,(v_stat[j]));
                v_summed_hist[i].SetBinContent(bincount,sqrt(summed_up[j]-pow(v_stat[j],2)));
                v_mcstat_hist[i].SetBinContent(bincount,(v_mcstat[j]));
                //For individual stats
                std::cout<<"Stat "<<(v_stat[j])<<" Summed: "<<sqrt(summed_up[j])<<"  MC "<<(v_mcstat[j])<<std::endl;
                for(int k=0; k<v_sys.size();k++){
                     v_sys_hist[k][i].SetBinContent(bincount, (v_sys[k][j]));
                     std::cout<<" Sys "<<k<<" "<<(v_sys[k][j])<<std::endl;
                }
                bincount++;
            }

            TLegend *l = new TLegend(0.11,0.69,0.89,0.89);
            l->SetNColumns(2);
            l->SetLineWidth(0);
            l->SetLineColor(kWhite);
            l->SetFillStyle(0);

            v_summed_hist[i].SetLineColor(kBlack);
            v_summed_hist[i].SetLineWidth(2);
            v_summed_hist[i].Scale(100.0);
            v_summed_hist[i].Draw("hist");
            v_summed_hist[i].GetYaxis()->SetTitle("Fractional Error (%)");
//            v_summed_hist[i].GetYaxis()->SetTitleOffset(0.5);
            v_summed_hist[i].GetXaxis()->SetTitle((sig.channel_units.at(i).c_str()));
            v_summed_hist[i].SetMinimum(0);
            //v_summed_hist[i].SetMaximum(v_summed_hist[i].GetMaximum()*1.4);
            v_summed_hist[i].SetMaximum(65);
            l->AddEntry(&v_summed_hist[i],"Total Systematics","l");

            v_stat_hist[i].SetLineColor(kGray);
            v_stat_hist[i].SetLineWidth(2);
            v_stat_hist[i].Scale(100);
            v_stat_hist[i].Draw("hist same");
            v_stat_hist[i].SetLineStyle(9);
            
            v_mcstat_hist[i].SetLineWidth(2);
            v_mcstat_hist[i].SetLineColor(kMagenta);
            v_mcstat_hist[i].Scale(100);
            v_mcstat_hist[i].Draw("hist same");
            l->AddEntry(&v_mcstat_hist[i],"Intrinsic MC Stats","l");
            
            for(int k=0; k<v_sys.size();k++){
                     v_sys_hist[k][i].SetLineColor(cols[k]);
                     v_sys_hist[k][i].Scale(100);
                     v_sys_hist[k][i].Draw("hist same");
                     v_sys_hist[k][i].SetLineWidth(2);
                     l->AddEntry(&v_sys_hist[k][i],(nams[k].c_str()),"l");
            
            }
            l->AddEntry(&v_stat_hist[i],"Data Sized Stats","l");
            l->Draw();
            c->Update();
            c->SaveAs(("AnalyzeSys_"+tag+"_"+std::to_string(i)+".pdf").c_str(),"pdf");
    }



    std::cout << "Total wall time: " << difftime(time(0), start_time)/60.0 << " Minutes.\n";
    return 0;

}
