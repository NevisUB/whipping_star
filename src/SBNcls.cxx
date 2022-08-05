#include "SBNcls.h"
using namespace sbn;

int SBNcls::SetSamplePoisson(){
    which_sample = 0;
   return which_sample;
}

int SBNcls::SetSampleCovariance(){
    which_sample = 1;
   return which_sample;
}


double SBNcls::pval2sig1sided(double pval){
   return sqrt(2)*TMath::ErfInverse(1-pval);
}

double SBNcls::pval2sig2sided(double pval){
   return sqrt(2)*TMath::ErfInverse(1-pval*2.0);
}

double SBNcls::pval2sig(double pval){
   //for backward compatability, work
   return pval2sig1sided(pval);
}

int SBNcls::setMode(int input_mode){
    which_mode = input_mode;
    return which_mode;

}

int SBNcls::CalcCLS(int numMC, std::string tag){


	if(which_sample == 0){
		std::cout<<"SBNcls::CalcCLS\t|| Running in Poission sampling mode!"<<std::endl;
	}
	 else if(which_sample ==1){ 
		std::cout<<"SBNcls::CalcCLS\t|| Running in Covariance sampling mode!"<<std::endl;
	}
	
    std::vector<double> ven={0};
	std::vector<double>  vec_CLs;

    time_t start_time = time(0);
	
	double N_h1 = h1->GetTotalEvents();
	double N_h0 = h0->GetTotalEvents();

    double central_value_chi = chi_h0.CalcChi(h1);
    double central_value_chi_h1 = chi_h1.CalcChi(h1);


    std::cout<<"SBNcls::CalcCLS\t|| Central Value Chi is : "<<central_value_chi<<" and "<<central_value_chi_h1<<" should be 0"<<std::endl;
	TH1D h1_pdf;


    //mode =0: absolute chi2 value
    if(which_mode==0){
	//distribution of chi2 values (from fake experiments), given a cv spectra

	//core spec is h0, random possion sampling from h1, and calculate the chi2 value of random h1 sample and h0.
        if(which_sample == 0) h1_pdf = chi_h0.SamplePoissonVaryInput(h1, numMC, central_value_chi*50);

	//core spec is h0, sampling from h1 using cholosky decomposition, and calculate the chi2 value of random h1 and h0
        else if(which_sample==1) h1_pdf = chi_h0.SampleCovarianceVaryInput(h1, numMC, central_value_chi*50);

    //mode =1: delta chi2 value
    }else if (which_mode ==1){
        h1_pdf = chi_h1.SamplePoisson_NP(h1,chi_h0,chi_h1, numMC, central_value_chi*50,which_sample);
    }
	
	double sig1 = 0.5-(0.6827)/2.0;
	double sig2 = 0.5-(0.9545)/2.0;

	std::vector<double> prob_values = {1-sig2, 1-sig1, 0.5, sig1, sig2};
	std::vector<double> quantiles(prob_values.size());	
	h1_pdf.ComputeIntegral(); 
	h1_pdf.GetQuantiles(prob_values.size(), &quantiles[0], &prob_values[0]);
	std::cout<<" BKG Quantile: "<<quantiles[0]<<" "<<quantiles[1]<<" "<<quantiles[2]<<std::endl;

	//Now calculate the pvalues associated with those h1 variations. 
	TH1D h0_pdf;
	std::vector<double> pval = quantiles;
    pval.push_back(central_value_chi); 

    if(which_mode ==0){
    	if(which_sample == 0) h0_pdf = chi_h0.SamplePoissonVaryInput(h0, numMC, &pval);
        else if(which_sample ==1) h0_pdf = chi_h0.SampleCovarianceVaryInput(h0, numMC, &pval);
    }else if (which_mode==1){
    	h0_pdf = chi_h0.SamplePoisson_NP(h0, chi_h0,chi_h1, numMC, &pval,which_sample);
    }

	//lets do CLs
	for(int p=0; p<pval.size()-1;p++){
		vec_CLs.push_back(pval.at(p)/(1-prob_values.at(p)) );
	}

	std::cout<<"NUMBER @ scale : "<<"\t"<< pval.at(0)/(1-prob_values.at(0))<<"\tsig: "<<N_h0<<"\tbkg: "<<N_h1<<"\tsig_only "<<N_h0-N_h1<<"\t\ts/sqrt(s+b): "<<(N_h0-N_h1)/sqrt(N_h0)<<std::endl;


    std::cout << "Total wall time: " << difftime(time(0), start_time)/1.0 << " Secs.\n";


    return makePlots( h0_pdf,h1_pdf, tag,pval,which_mode);

}


int SBNcls::makePlots(TH1D & h0_pdf, TH1D & h1_pdf, std::string tag, std::vector<double> & pval, int which_mode){
    std::vector<double> ven={0};
	std::vector<double>  vec_CLs;

    
    double sig1 = 0.5-(0.6827)/2.0;
	double sig2 = 0.5-(0.9545)/2.0;

	std::vector<double> prob_values = {1-sig2, 1-sig1, 0.5, sig1, sig2};
	std::vector<double> quantiles(prob_values.size());	
	h1_pdf.GetQuantiles(prob_values.size(), &quantiles[0], &prob_values[0]);
	std::cout<<" BKG Quantile: "<<quantiles[0]<<" "<<quantiles[1]<<" "<<quantiles[2]<<std::endl;


	//lets do CLs
	for(int p=0; p<pval.size()-1;p++){
		vec_CLs.push_back(pval.at(p)/(1-prob_values.at(p)) );
	}

	TFile * fp = new TFile(("SBNfit_CLs_"+tag+".root").c_str(),"recreate");
	fp->cd();
	TCanvas *cp=new TCanvas();

	h0_pdf.SetStats(false);
	h1_pdf.SetStats(false);

    h0_pdf.Scale(1.0/h0_pdf.Integral("width"));
	h1_pdf.Scale(1.0/h1_pdf.Integral("width"));

	h0_pdf.SetLineColor(kRed-7);
	h1_pdf.SetLineColor(kBlue-4);
	h0_pdf.SetFillColor(kRed-7);
	h1_pdf.SetFillColor(kBlue-4);
	h0_pdf.SetFillStyle(3445);
	h1_pdf.SetFillStyle(3454);

	h0_pdf.Draw("hist");
		
	double maxval =std::max(  h0_pdf.GetMaximum(),h1_pdf.GetMaximum());
	double minval = 0;
	std::cout<<"SBNcls::CalcCLS() || Minimum value: "<<minval<<" Maximum value: "<<maxval<<std::endl;
	h0_pdf.SetMinimum(minval);
	h0_pdf.SetMaximum(maxval*1.35);

	double minbin = std::min(h0_pdf.GetBinLowEdge(h0_pdf.FindFirstBinAbove(0)), h1_pdf.GetBinLowEdge(h1_pdf.FindFirstBinAbove(0)));
	double maxbin = std::max(h0_pdf.GetBinLowEdge(h0_pdf.FindLastBinAbove(0))+h0_pdf.GetBinWidth(h0_pdf.FindLastBinAbove(0)), h1_pdf.GetBinLowEdge(h1_pdf.FindLastBinAbove(0))+h1_pdf.GetBinWidth(h1_pdf.FindLastBinAbove(0)));

	std::cout<<"SBNcls::CalcCLS() || Minimum Bin: "<<minbin<<" Maximum bin: "<<maxbin<<std::endl;
	h0_pdf.GetXaxis()->SetRangeUser(minbin<0 ? minbin*1.2: minbin*0.8,maxbin*0.8);

    bool draw_both = true;

	std::vector<std::string> quantile_names = {"-2#sigma","-1#sigma","Median","+1#sigma","+2#sigma"};
	std::vector<int> cols ={kYellow-7, kGreen+1, kBlack,kGreen+1, kYellow-7};	

    if(draw_both){
	for(int i=0; i< quantiles.size(); i++){	
		TLine *l = new TLine(quantiles.at(i),0.0, quantiles.at(i),maxval*1.05);
		l->SetLineColor(cols.at(i));
		l->SetLineWidth(2);
	        TLatex * qnam = new TLatex();
   	        qnam->SetTextSize(0.045);
   		qnam->SetTextAlign(12);  //align at top
		qnam->SetTextAngle(-90);
 		qnam->DrawLatex(quantiles.at(i), maxval*1.3 ,quantile_names.at(i).c_str());
		l->Draw("same");
	
		TLatex * qvals = new TLatex();
		qvals->SetTextSize(0.03);
		qvals->SetTextAlign(32);
        double sigma_val = pval2sig(pval.at(i));
        std::string whatsigma = sbnfit_to_string_prec(sigma_val,1)+"#sigma";
        if(sigma_val==0.0){
            whatsigma = "inf ";
        }

		std::string details =  ("#splitline{"+quantile_names.at(i)+"}{1-#beta(" +sbnfit_to_string_prec(1-prob_values.at(i),3) + ") #alpha("+ sbnfit_to_string_prec(pval.at(i),3) +" | "+whatsigma+ ") CL_{s}("+sbnfit_to_string_prec(vec_CLs.at(i),3)+")}");
		std::string details2 =  ("#splitline{"+quantile_names.at(i)+"}{1-#beta(" +sbnfit_to_string_prec(1-prob_values.at(i),10) + ") #alpha("+ sbnfit_to_string_prec(pval.at(i),10) +" | "+sbnfit_to_string_prec(pval2sig(pval.at(i)),1)+ "#sigma) CL_{s}("+sbnfit_to_string_prec(vec_CLs.at(i),10)+")}");
		std::cout<<details2<<std::endl;
		qvals->DrawLatexNDC(0.875, 0.2+i*0.1,details.c_str()  );
	}
    }

    /*
    TLine *lcv = new TLine(central_value_chi,minval,central_value_chi, maxval);
    lcv->SetLineColor(kBlack);
    lcv->SetLineStyle(2);
    lcv->SetLineWidth(2);
    TLatex * cvnam = new TLatex();
    cvnam->SetTextSize(0.045);
   	cvnam->SetTextAlign(12);  //align at top    
    cvnam->SetTextAngle(-90);
 	cvnam->DrawLatex(central_value_chi, maxval*1.3 ,"CV");
    lcv->Draw("same");
	*/

    std::string cv_details =  ("#splitline{CV}{#alpha("+ sbnfit_to_string_prec(pval.back(),5) +" | "+sbnfit_to_string_prec(pval2sig(pval.back()),5)+ "#sigma)}");
    std::cout<<cv_details<<std::endl;	

    //chi^2 prob bit
    std::vector<double> analytical_chi;
    std::vector<double> analytical_prob;

    double analytical_sum =0;
    for(double t=0; t< maxbin*10; t+=0.01){
        analytical_chi.push_back(t);
        analytical_prob.push_back( gsl_ran_chisq_pdf(t,h1->num_bins_total_compressed)   );
        analytical_sum +=0.01*analytical_prob.back();
    }

    for(auto &p:analytical_prob){
    //        p =p/analytical_sum;
    }

    TGraph *analytical_graph = new TGraph(analytical_chi.size(),&analytical_chi[0],&analytical_prob[0]);
    analytical_graph->SetLineColor(kRed);
    if(which_mode==0)analytical_graph->Draw("same");


	TLegend *leg = new TLegend(0.7,0.7,0.89,0.89);
	leg->SetLineWidth(0);
	leg->SetFillStyle(0);
	leg->AddEntry(&h0_pdf,"H_{0}","lf");
	if(draw_both)leg->AddEntry(&h1_pdf,"H_{1}","lf");
    if(which_mode==0)leg->AddEntry(analytical_graph,("#chi^{2} PDF "+std::to_string(h0->num_bins_total_compressed)+" dof").c_str(),"l");
	leg->Draw();

    if(which_mode==0){
    	h0_pdf.GetXaxis()->SetTitle("#chi^{2}");
    }else if(which_mode==1){
    	h0_pdf.GetXaxis()->SetTitle("#Delta #chi^{2} = #chi^{2}_{H0} - #chi^{2}_{H1} ");
    }
	h0_pdf.GetYaxis()->SetTitle("PDF");

	if(draw_both) h1_pdf.Draw("hist same");	
	cp->Write();	
	cp->SaveAs(("SBNfit_Cls_"+tag+".pdf").c_str(),"pdf");	
	fp->Close();
	
	return 0;

}



/*double SBNcls::NPCalculator(SBNchi * chi_H0, SBNchi *chi_H1, SBNspec * spec_H0, SBNspec * spec_H1){
    return 0;
}*/
