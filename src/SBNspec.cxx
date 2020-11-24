#include "SBNspec.h"
#include <cassert>
using namespace sbn;

SBNspec::SBNspec(std::string whichxml, int which_universe, bool isverbose) : SBNconfig(whichxml,isverbose){

	//Initialise all the things
	//for every multisim, create a vector of histograms, one for every subchannel we want
	std::cout << "check in SBNspec" << __LINE__ << std::endl;
	int ctr=0;
	for(auto fn: fullnames){
		for(int c=0; c<channel_names.size(); c++){
			if(fn.find("_"+channel_names[c]+"_")!=std::string::npos){
				double * tbins =&bin_edges[c][0];
				std::string thisname;
				if(which_universe<0){
					thisname = fn;
				}else{
					thisname = fn+"_MS"+std::to_string(which_universe);
				}
				TH1D thischan(thisname.c_str(),"",num_bins[c], tbins );
                //thischan.AddDirectory(false);
				hist.push_back(thischan);
				//hist.emplace_back(TH1D(thisname.c_str(), thisname.c_str(), num_bins[c], tbins ));
				//auto it = hist.begin()+ctr;
				//map_hist[fn] = &(*it);
				map_hist[fn] = ctr;

				ctr++;
			}
		}
	}
	std::cout << "check in SBNspec" << __LINE__ << std::endl;

	has_been_scaled = false;
	this->CollapseVector();

}

SBNspec::SBNspec(std::string whichxml): SBNspec(whichxml,-1,true){}
SBNspec::SBNspec(std::string whichxml, int which_universe): SBNspec(whichxml,which_universe, true){}

SBNspec::SBNspec(std::string rootfile, std::string whichxml) : SBNspec(rootfile, whichxml, true){}
SBNspec::SBNspec(std::string rootfile, std::string whichxml, bool isverbose) : SBNconfig(whichxml, isverbose) {
	//Contruct from a prexisting histograms that exist in a rootfile
	std::cout << "SBNspec || Opening root file " << rootfile << std::endl;
	TFile *f = new TFile(rootfile.c_str(),"read");

	//Loop over all filenames that should be there, and load up the histograms.
	int n=0;
	for(auto fn: fullnames){
		hist.push_back(*((TH1D*)f->Get(fn.c_str())));
		map_hist[fn] = n;
		n++;
	}

	has_been_scaled=false;


	f->Close();


}

SBNspec::SBNspec(std::vector<double> input_full_vec, std::string whichxml) : SBNspec(input_full_vec, whichxml, false){ };
SBNspec::SBNspec(std::vector<double> input_full_vec, std::string whichxml, bool isverbose) : SBNspec(input_full_vec,whichxml,-1,isverbose){};
SBNspec::SBNspec(std::vector<double> input_full_vec, std::string whichxml, int universe, bool isverbose) : SBNspec(whichxml,universe,isverbose){

	std::cout << "check in SBNspec" << __LINE__ << std::endl;
	for(int i=0; i< input_full_vec.size(); i++){

		int which_hist = GetHistNumber(i);

		int exact_bin = i;
		for(int b=0; b<which_hist; b++){
			exact_bin -= hist.at(b).GetNbinsX();
		}
		hist.at(which_hist).SetBinContent(exact_bin+1, input_full_vec.at(i));

	}


	this->CalcFullVector();
}







int SBNspec::Add(std::string which_hist, TH1 * histin){
	//Addes all hists toGether
	if(map_hist.count(which_hist) <= 0){
		std::cout<<"SBNspec::Add || ERROR the passed in histgram name is not one defined in the xml! passsed in: "<<which_hist<<" "<<map_hist.count(which_hist)<<std::endl;

		for(std::map<std::string, int >::iterator  it = map_hist.begin(); it != map_hist.end(); ++it){
			std::cout<<"SBNspec::Add || "<<it->first<<" @ "<<it->second<<std::endl;
		}


		exit(EXIT_FAILURE);
	}

	int h=map_hist[which_hist];

	if(hist.at(h).GetNbinsX() != histin->GetNbinsX()){
		std::cout<<"SBNspec::Add || ERROR the passed in histgram has different number of bins!"<<std::endl;
		std::cout<<"SBNspec::Add || "<<which_hist<<" : "<<hist.at(h).GetNbinsX()<<std::endl;
		std::cout<<"SBNspec::Add || Inputted hist : "<<histin->GetNbinsX()<<std::endl;
		exit(EXIT_FAILURE);
	}

	hist.at(h).Add(histin);

	this->CollapseVector();
	return 0;
}






int SBNspec::Add(SBNspec *in){
	//Addes all hists toGether
	if(xmlname != in->xmlname){ std::cout<<"ERROR: SBNspec::Add, trying to add differently configured SBNspecs!"<<std::endl; exit(EXIT_FAILURE);}

	for(int i=0; i< hist.size(); i++){
		hist[i].Add( &(in->hist[i]));
	}

	this->CollapseVector();
	return 0;
}

int SBNspec::SetAsGaussian(double mean, double sigma, int ngen){
	TRandom3 *seedGetter = new TRandom3(0);
	int seed = seedGetter->Integer(1000000);

	for(auto &h: hist){
		TRandom3 *rangen = new TRandom3(seed);
		h.Reset();
		for(int i=0; i<ngen; i++){
			double eve = rangen->Gaus(mean,sigma);
			h.Fill( eve );
		}
	}

	return 0;

}

int SBNspec::SetAsFlat(double val){
	for(auto &h: hist){
		for(int i=0; i<h.GetSize(); i++){
			h.SetBinContent(i, val );
		}
	}
}



//All scaling functions are quite self explanatory
int SBNspec::ScalePoisson(){
	TRandom3 *rangen = new TRandom3(0);
	for(auto &h: hist){
		for(int i=0; i<h.GetSize(); i++){
			h.SetBinContent(i, rangen->Poisson( h.GetBinContent(i)    ));
		}
	}
	return 0;
}

int SBNspec::ScalePoisson(TRandom3* rangen){
	for(auto &h: hist){
		for(int i=0; i<h.GetSize(); i++){
			h.SetBinContent(i, rangen->Poisson( h.GetBinContent(i)    ));
		}
	}
	return 0;
}



int SBNspec::ScaleRandom(){
	TRandom3 *rangen    = new TRandom3(0);

	for(auto& h: hist){
		h.Scale(rangen->Uniform(0,2));

	}
	return 0;
}


int SBNspec::Scale(std::string name, TF1 * func){
	for(auto& h: hist){
		std::string test = h.GetName();
		if(test.find(name)!=std::string::npos ){
			for(int b=0; b<=h.GetNbinsX(); b++){
				//std::cout<<h.GetBinContent(b)<<" "<<h.GetBinCenter(b)<<" "<<func->Eval(h.GetBinCenter(b) )<<std::endl;
				h.SetBinContent(b, h.GetBinContent(b)*func->Eval(h.GetBinCenter(b) ) );
			}
		}

	}

	return 0;
}


int SBNspec::ScaleAll(double sc){
	for(auto& h: hist){
		h.Scale(sc);
	}
	this->CollapseVector();

	return 0;
}

int SBNspec::Scale(std::string name, double val){
	for(auto& h: hist){
		std::string test = h.GetName();

		if(test.find(name)!=std::string::npos){
			//	std::cout<<name<<". found in: "<<test<<" at "<<test.find(name)<<std::endl;
			h.Scale(val);

            //std::cout<<"scaled "<<name<<" by "<<val<<std::endl;
		}

	}

	has_been_scaled = true;
	scale_hist_name =name;
	scale_hist_val = val;

	this->CollapseVector();
	return 0;
}

int SBNspec::NormAll(double n){

	for(auto& h: hist) {
		h.Scale(n/h.GetSumOfWeights());
	}
	return 0;
}

int SBNspec::Norm(std::string name, double val){
	for(auto& h: hist){
		std::string test = h.GetName();

		if(test.find(name)!=std::string::npos ){
			//std::cout<<name<<". found in: "<<test<<" at "<<test.find(name)<<std::endl;
			h.Scale(val/h.GetSumOfWeights());
		}

	}
	return 0;

}





int SBNspec::CalcFullVector(){
  full_vector.clear();

  full_vector.resize(num_bins_total);

  int hoffset = 0;
  for(size_t hid=0; hid<hist.size(); ++hid) {
    const auto& h =  hist[hid];
    for(int i = 1; i < (h.GetSize()-1); ++i){
      full_vector[hoffset + i - 1] = h.GetBinContent(i);
    }
    hoffset += (h.GetSize()-2);
  }
    
  assert (hoffset == num_bins_total);

  return 0;
}

int SBNspec::CollapseVector(){

	collapsed_vector.clear();
	f_collapsed_vector.clear();
	CalcFullVector();

	for(int im = 0; im < num_modes; im++){
		for(int id =0; id < num_detectors; id++){
			int edge = id*num_bins_detector_block + num_bins_mode_block*im; // This is the starting index for this detector

			for(int ic = 0; ic < num_channels; ic++){
				int corner=edge;

				for(int j=0; j< num_bins.at(ic); j++){

					double tempval=0;

					for(int sc = 0; sc < num_subchannels.at(ic); sc++){

						//std::cout<<im<<"/"<<num_modes<<" "<<id<<"/"<<num_detectors<<" "<<ic<<"/"<<num_channels<<" "<<j<<"/"<<num_bins[ic]<<" "<<sc<<"/"<<num_subchannels[ic]<<std::endl;
						tempval += full_vector.at(j+sc*num_bins.at(ic)+corner);
						edge +=1;	//when your done with a channel, add on every bin you just summed
					}
					//we can size this vector beforehand and get rid of all push_back()
					collapsed_vector.push_back(tempval);
					f_collapsed_vector.push_back((float)tempval);
				}
			}
		}
	}
	return 0;
}

double SBNspec::GetTotalEvents(){
	double ans =0;
	this->CalcFullVector();

	for(double d: full_vector){
		ans+=d;
	}

	return ans;


}

int SBNspec::PrintFullVector(){
	for(double d: full_vector){
		std::cout<<d<<" ";
	}
	std::cout<<std::endl;
	return 0;
}

int SBNspec::PrintCollapsedVector(){
	for(double d: collapsed_vector){
		std::cout<<d<<" ";
	}
	std::cout<<std::endl;
	return 0;
}


int SBNspec::WriteOut(std::string tag){
	//kWhite  = 0,   kBlack  = 1,   kGray    = 920,  kRed    = 632,  kGreen  = 416,
	//kBlue   = 600, kYellow = 400, kMagenta = 616,  kCyan   = 432,  kOrange = 800,
	//kSpring = 820, kTeal   = 840, kAzure   =  860, kViolet = 880,  kPink   = 900

	std::vector<int> mycol = {kGreen+1, kRed-7, kBlue-4, kOrange-3, kMagenta+1, kCyan-3,kYellow, kGreen-3 };
	int colindex =0;
	TFile *f2 = new TFile((tag+".SBNspec.root").c_str(),"recreate");
    f2->cd();
    std::cout<<"This file has "<<hist.size()<<" histograms and "<<fullnames.size()<<" fullnames "<<std::endl;
	for(auto& h: hist){
        std::cout<<"WriteOut "<<h.Integral()<<std::endl;
		h.Write();
	}
	f2->Close();

	TFile *f = new TFile(("SBNfit_spectrum_plots_"+tag+".root").c_str(),"RECREATE");
	f->cd();

	std::vector<TH1D> temp_hists = hist;


	for(int im = 0; im <mode_names.size(); im++){
		for(int id = 0; id <detector_names.size(); id++){
			for(int ic = 0; ic <channel_names.size(); ic++){


				std::string canvas_name = mode_names.at(im)+"_"+detector_names.at(id)+"_"+channel_names.at(ic);


				bool this_run = false;

				TCanvas* Cstack= new TCanvas(canvas_name.c_str(),canvas_name.c_str());
				Cstack->SetFixedAspectRatio();

				Cstack->cd();
				THStack * hs = new THStack(canvas_name.c_str(),  canvas_name.c_str());
				TLegend legStack(0.65,0.25,0.89,0.89);
				legStack.SetLineWidth(0);
				legStack.SetLineColor(kWhite);
				int n=0;
				int nc=0;
				TH1D *hsum;

				//Ok to sort the histograms based on "size" will need a few tricks
				std::vector<double> integral_sorter;
				std::vector<TH1*> to_sort;
				std::vector<std::string> l_to_sort;

				for(auto &h : temp_hists){
					std::string test = h.GetName();
					if(test.find(canvas_name)!=std::string::npos ){

						double total_events = h.GetSumOfWeights();
						h.Sumw2(false);
						h.Scale(1,"width,nosw2");
						h.GetYaxis()->SetTitle("Events/GeV");
						h.SetMarkerStyle(20);
						h.SetMarkerColor(mycol[n]);
						h.SetFillColor(mycol[n]);
						h.SetLineColor(kBlack);
						h.SetTitle(h.GetName());
						//h.Write();

						if(!this_run){
							hsum = (TH1D*)h.Clone(("sum_"+canvas_name).c_str());
							hsum->Reset();
						}

						std::ostringstream out;
						out <<std::fixed<< std::setprecision(3) << total_events;
						std::string hmm = " \t ";
						std::string tmp =h.GetName() +hmm+ out.str();
						//std::string tmp = map_subchannel_plotnames.at(h.GetName()) +hmm+ out.str();
						//legStack.AddEntry(&h, tmp.c_str() , "f");
						hsum->Add(&h);
						//hs->Add(&h);
						n++;

						this_run=true;

						to_sort.push_back(&h);
						l_to_sort.push_back(tmp);
						integral_sorter.push_back(total_events);

					}
				}
				//Sort!
				for (int i: SortIndexes(integral_sorter)) {
					hs->Add(to_sort.at(i));
					legStack.AddEntry(to_sort.at(i), l_to_sort.at(i).c_str(),"f");
				}



				if(this_run){
					double plot_pot=5e19;

					double title_size_ratio=0.1;
					double label_size_ratio=0.1;
					double title_offSet_ratioY = 0.3 ;
					double title_offSet_ratioX = 1.1;

					double title_size_upper=0.15;
					double label_size_upper=0.05;
					double title_offSet_upper = 1.45;

					Cstack->cd();
					hs->Draw();

					hs->GetYaxis()->SetTitle(("Events/"+channel_units.at(ic)).c_str());
					hs->GetXaxis()->SetTitle(channel_units.at(ic).c_str());


					//hcomp->Draw("hist same");
					hs->SetMaximum(hs->GetMaximum()*1.1);
					hs->SetMinimum(0.001);

					Cstack->Update();
					legStack.Draw();

					f->cd();

					Cstack->Write(canvas_name.c_str() );

				}





			}
		}
	}

	f->Close();

	return 0;
}


//given collapsed systematic covariance matrix, compare SBNspec with another one with error bars shown
int SBNspec::CompareSBNspecs(TMatrixT<double> collapse_covar, SBNspec * compsec, std::string tag){
	std::vector<int> mycol = {kAzure -9, kRed-7, kGreen-3, kBlue-6, kMagenta-3, kYellow-7,  kOrange-3, kBlue, kBlue+2,  kGreen+1,kBlue-7, kPink, kViolet, kCyan,kMagenta,kAzure};

	if(collapse_covar.GetNcols() != this->num_bins_total_compressed){
                std::cout << "ERROR: dimension of covariance matrix " << collapse_covar.GetName() << " is diffrent from num_bins_total_compressed" << std::endl;
                exit(EXIT_FAILURE);
        }



        TFile *f = new TFile(("SBNfit_compare_plots_"+tag+".root").c_str(),"RECREATE");
        f->cd();


        std::vector<TH1D> temp = hist;
        std::vector<TH1D> temp_comp = compsec->hist;

        for(int k=0; k< fullnames.size(); k++){
                TCanvas *ctmp = new TCanvas((tag+"_"+std::to_string(k)+"_"+fullnames.at(k)).c_str(), (std::to_string(k)+"_"+fullnames.at(k)).c_str(),1200,1200);
                ctmp->cd();

        if(map_hist.count(fullnames.at(k)) <=0 ){
                std::cout<<"WARNING:  hist "<<fullnames[k]<<" is not in map_hist"<<std::endl;
        }
                TH1D * h1 = (TH1D*) temp.at(map_hist[fullnames.at(k)]).Clone((std::to_string(k)+fullnames.at(k)+"_1").c_str());
                TH1D * h2 = (TH1D*) temp_comp.at(map_hist[fullnames.at(k)]).Clone((std::to_string(k)+fullnames.at(k)+"_2").c_str());
	
		h1->SetLineColor(kRed);
                h1->SetLineWidth(2);
                h1->Draw("hist");

                h2->SetLineColor(kBlue);
                h2->SetLineStyle(5);
                h2->SetLineWidth(2);
                h2->Draw("hist same");

                h1->SetMaximum(std::max(h1->GetMaximum(), h2->GetMaximum())*1.10);

                ctmp->Write();
        }



        int error_bin=0;

        for(auto m: mode_names){
                for(auto d: detector_names){
                        for(auto c: channel_names){

                                std::string canvas_name = m+"_"+d+"_"+c;

                                bool this_run = false;
                                bool this_run_comp = false;

                                TCanvas* Cstack= new TCanvas((tag+"_"+canvas_name).c_str(),canvas_name.c_str(),1200,1200);
                                //Cstack->SetFixedAspectRatio();
                                Cstack->cd();
                                THStack * hs = new THStack(canvas_name.c_str(),  canvas_name.c_str());
				TLegend legStack(0.11, 0.7, 0.80, 0.89);
                                legStack.SetNColumns(2);
                                legStack.SetLineWidth(0);
                                legStack.SetLineColor(kWhite);
                                int n=0;
                                int nc=0;
                                TH1D * hcomp;
                                TH1D *hsum;
				double hcomp_sum=0;
                                double hsum_sum=0;


                                for(auto &h : temp_comp){
                                        std::string test = h.GetName();
                                        if(test.find(canvas_name)!=std::string::npos){
                                                double total_events = h.GetSumOfWeights();
						hcomp_sum += total_events;
						//h.Scale(1,"width,nosw2");   //not care about the statistics
						h.Scale(1,"width");
						//h.GetYaxis()->SetTitle("Events/GeV");
						h.SetLineColor(kBlack);
						if(!this_run_comp){
                                                        hcomp = (TH1D*)h.Clone(("comp_"+canvas_name).c_str());
                                                        hcomp->Reset();
                                                }

                                                std::ostringstream out;
                                                out << std::setprecision(3) << total_events;
                                                std::string hmm = "\t";
                                                std::string tmp = h.GetName() +hmm+ out.str();


                                                hcomp->Add(&h);
                                                nc++;

                                                this_run_comp=true;

                                        }
                                }


                                //Ok to sort the histograms based on "size" will need a few tricks
                                std::vector<double> integral_sorter;
                                std::vector<TH1*> to_sort;
                                std::vector<std::string> l_to_sort;

                                for(auto &h : temp){
                                        std::string test = h.GetName();
                                        if(test.find(canvas_name)!=std::string::npos ){

                                                double total_events = h.GetSumOfWeights();
						hsum_sum += total_events;
						//h.Scale(1,"width,nosw2");
						h.Scale(1,"width");
						h.GetYaxis()->SetTitle("Events/GeV");
                                                h.SetMarkerStyle(20);
                                                h.SetMarkerColor(mycol[n]);
                                                h.SetFillColor(mycol[n]);
                                                h.SetLineColor(kBlack);
                                                h.SetTitle(h.GetName());
						if(!this_run){
                                                        hsum = (TH1D*)h.Clone(("sum_"+canvas_name).c_str());
                                                        hsum->Reset();
                                                }

                                                std::ostringstream out;
                                                out <<std::fixed<< std::setprecision(0) << total_events;
                                                std::string hmm = " | ";
                                                std::string tmp = h.GetName()+hmm+ out.str();
						hsum->Add(&h);
						n++;

                                                this_run=true;

                                                to_sort.push_back(&h);
                                                l_to_sort.push_back(tmp);
                                                integral_sorter.push_back(total_events);

                                        }
                                }
                                //Sort!
                                for (int i: SortIndexes(integral_sorter)) {
                                        hs->Add(to_sort.at(i), "HIST");
                                        legStack.AddEntry(to_sort.at(i), l_to_sort.at(i).c_str(),"f");
                                }

				legStack.AddEntry(hsum, Form("MC Stack | %.2f", hsum_sum), "fl");
                                legStack.AddEntry(hcomp, Form("Data | %.2f", hcomp_sum), "flp");

                                //set the error of 'this' spec according to the covariance matrix
                                for(int i=0; i<hsum->GetNbinsX(); i++){
					double xbin_width = hsum->GetXaxis()->GetBinWidth(i+1);
					//systematic errors from diag element need to be scaled by 1/width too
					//if covariance matrix only include systematic error; add stats error and syst error in quadrature
                                        double error = sqrt(pow(hsum->GetBinError(i+1), 2.0) + collapse_covar(error_bin+i, error_bin+i)/pow(xbin_width, 2.0));
                                        //double error = hsum->GetBinError(i+1) + sqrt(collapse_covar(error_bin+i, error_bin+i))/xbin_width;
					//if covariance matrix includes both stats and syst
                                        //double error = sqrt(collapse_covar(error_bin+i, error_bin+i))/xbin_width;

                                        //std::cout << "Stats only error: "<< hsum->GetBinError(i+1) << ", Syst+Stat one: " << error << std::endl;
					//std::cout << "Bin content: " << hsum->GetBinContent(i+1) << std::endl;
                                        hsum->SetBinError(i+1, error);
                                }
                                error_bin +=hsum->GetNbinsX();


				//start drawing plots!!!
				if(this_run && this_run_comp){
                                        double plot_pot=5e19;

                                        double title_size_ratio=0.115;
                                        double label_size_ratio=0.11;
                                        double title_offSet_ratioY = 0.35;
                                        double title_offSet_ratioX = 1.1;

                                        double title_size_upper=0.05;
                                        double label_size_upper=0.05;
                                        double title_offSet_upper = 0.85;


                                        Cstack->cd();
                                        //gStyle->SetErrorX(0);
                                        TPad *pad0top = new TPad(("pad0top_"+canvas_name).c_str(), ("pad0top_"+canvas_name).c_str(), 0, 0.35, 1, 1.0);
                                        pad0top->SetBottomMargin(0); // Upper and lower plot are joined
					pad0top->Draw();             // Draw the upper pad: pad2top
					pad0top->cd();               // pad2top becomes the current pad
					hs->Draw();
                                        //draw error bar for 'this' SBNspec
                                        hsum->SetFillColor(kBlack);
                                        hsum->SetFillStyle(3145);
                                        hsum->Draw("E2 same");
                                        hs->GetYaxis()->SetTitle("Events/GeV");
					hs->GetYaxis()->SetTitleSize(title_size_upper);
					hs->GetYaxis()->SetLabelSize(label_size_upper);
					hs->GetYaxis()->SetTitleOffset(title_offSet_upper);
                                        //hs->GetYaxis()->SetTitle("Events");
                                        
					//draw compared histogram
					hcomp->SetMarkerStyle(8);
                                        hcomp->SetMarkerColor(kBlack);
                                        hcomp->SetMarkerSize(1.5);
                                        hcomp->SetLineWidth(3);
                                        hcomp->SetLineColor(kBlack);
                                        hcomp->Draw("EP same");

                                        hs->SetMaximum(std::max(hs->GetMaximum(), hcomp->GetMaximum())*1.4);
                                        //hs->SetMaximum(std::max(hs->GetMaximum(), hcomp->GetMaximum())*1.1);
                                        hs->SetMinimum(0.001);

                                        Cstack->Update();
                                        legStack.Draw();

                                        Cstack->cd();
                                        gStyle->SetOptStat(0);
                                        TPad *pad0bot = new TPad(("padbot_"+canvas_name).c_str(),("padbot_"+canvas_name).c_str(), 0, 0.05, 1, 0.35);
                                        pad0bot->SetTopMargin(0);
                                        pad0bot->SetBottomMargin(0.351);
                                        pad0bot->SetGridx(); // vertical grid
					pad0bot->Draw();
                                        pad0bot->cd();       // pad0bot becomes the current pad
	
						
					//spectra ratio 
					TH1* ratpre = (TH1*)hcomp->Clone(("ratio_"+canvas_name).c_str());
					//ratpre->Sumw2();
                                        ratpre->Divide(hsum);
                                        ratpre->SetStats(false);

					//to draw the 1 sigma error band on the ratio plot
					TH1* h_err = (TH1*)hsum->Clone("error_band");
					h_err->Divide(hsum);
					//h->Divide() give sqrt(2) of the actual error size, so setup error bars manually
					for(int i=0; i<h_err->GetNbinsX(); i++){
						h_err->SetBinError(i+1, hsum->GetBinError(i+1)/hsum->GetBinContent(i+1));
						ratpre->SetBinError(i+1, hcomp->GetBinError(i+1)/hsum->GetBinContent(i+1));
					}

					ratpre->Draw("E");
                                        //ratpre->Draw("hist");
  					h_err->Draw("E2 same");
                                        ratpre->SetFillColor(kWhite);
                                        ratpre->SetFillStyle(0);
                                        ratpre->SetLineWidth(2);

                                        gStyle->SetOptStat(0);
                                        TLine *line = new TLine(ratpre->GetXaxis()->GetXmin(),1.0,ratpre->GetXaxis()->GetXmax(),1.0 );
                                        line->Draw("same");
                                        ratpre->SetLineColor(kBlack);
                                        ratpre->SetTitle("");
                                        ratpre->GetYaxis()->SetTitle("Ratio");
                                        ratpre->GetXaxis()->SetTitleOffset(title_offSet_ratioX);
                                        ratpre->GetYaxis()->SetTitleOffset(title_offSet_ratioY);
                                        ratpre->SetMinimum(std::min(0.5, ratpre->GetMinimum())*0.8);
                                        //ratpre->SetMinimum(0.0);
                                        ratpre->SetMaximum(std::max(1.0, ratpre->GetMaximum())*1.2);
                                        //ratpre->SetMaximum(2.0);
					ratpre->GetYaxis()->SetNdivisions(505, kTRUE);   //change the label division in y axis
                                        ratpre->GetYaxis()->SetTitleSize(title_size_ratio);
                                        ratpre->GetXaxis()->SetTitleSize(title_size_ratio);
                                        ratpre->GetYaxis()->SetLabelSize(label_size_ratio);
                                        ratpre->GetXaxis()->SetLabelSize(label_size_ratio);
                                        ratpre->GetXaxis()->SetTitle("Reconstructed Energy [GeV]");
                                        //ratpre->GetXaxis()->SetTitle("cos(#theta_{#pi})"); // for NC pi0 fit
					Cstack->Write(canvas_name.c_str() );
                    Cstack->SaveAs((tag+"_"+canvas_name+".pdf").c_str(),"pdf");

                                }

                        }
                }
        }

        f->Close();

        return 0;
}


int SBNspec::CompareSBNspecs(SBNspec * compsec, std::string tag){
	//kWhite  = 0,   kBlack  = 1,   kGray    = 920,  kRed    = 632,  kGreen  = 416,
	//kBlue   = 600, kYellow = 400, kMagenta = 616,  kCyan   = 432,  kOrange = 800,
	//kSpring = 820, kTeal   = 840, kAzure   =  860, kViolet = 880,  kPink   = 900
	std::vector<int> mycol = {kAzure -9, kRed-7, kGreen-3, kBlue-6, kMagenta-3, kYellow-7,  kOrange-3, kBlue, kBlue+2,  kGreen+1,kBlue-7, kPink, kViolet, kCyan,kMagenta,kAzure};


	TFile *f = new TFile(("SBNfit_compare_plots_"+tag+".root").c_str(),"RECREATE");
	f->cd();


	std::vector<TH1D> temp = hist;
	std::vector<TH1D> temp_comp = compsec->hist;

	for(int k=0; k< fullnames.size(); k++){
		TCanvas *ctmp = new TCanvas((tag+"_"+std::to_string(k)+"_"+fullnames.at(k)).c_str(), (std::to_string(k)+"_"+fullnames.at(k)).c_str(),1200,1200);
		ctmp->cd();

        if(map_hist.count(fullnames.at(k)) <=0 ){
                std::cout<<"WARNING:  hist "<<fullnames[k]<<" is not in map_hist"<<std::endl;
        }
		TH1D * h1 = (TH1D*) temp.at(map_hist[fullnames.at(k)]).Clone((std::to_string(k)+fullnames.at(k)+"_1").c_str());
		TH1D * h2 = (TH1D*) temp_comp.at(map_hist[fullnames.at(k)]).Clone((std::to_string(k)+fullnames.at(k)+"_2").c_str());

		h1->Scale(1,"width");
		h2->Scale(1,"width");

		h1->SetLineColor(kRed);
		h1->SetLineWidth(2);
		h1->Draw("hist");

		h2->SetLineColor(kBlue);
		h2->SetLineStyle(5);
		h2->SetLineWidth(2);
		h2->Draw("hist same");

		h1->SetMaximum(std::max(h1->GetMaximum(), h2->GetMaximum())*1.10);

		ctmp->Write();
	}




	for(auto m: mode_names){
		for(auto d: detector_names){
			for(auto c: channel_names){

				std::string canvas_name = m+"_"+d+"_"+c;

				bool this_run = false;
				bool this_run_comp = false;

				TCanvas* Cstack= new TCanvas((tag+"_"+canvas_name).c_str(),canvas_name.c_str(),1200,1200);
				//Cstack->SetFixedAspectRatio();

				Cstack->cd();
				THStack * hs = new THStack(canvas_name.c_str(),  canvas_name.c_str());
				TLegend legStack(0.11, 0.7, 0.80, 0.89);
				legStack.SetNColumns(2);
                                legStack.SetLineWidth(0);
				legStack.SetLineColor(kWhite);
				int n=0;
				int nc=0;
				TH1D * hcomp;
				TH1D *hsum;



				for(auto &h : temp_comp){
					std::string test = h.GetName();
					if(test.find(canvas_name)!=std::string::npos){
						double total_events = h.GetSumOfWeights();


						h.Sumw2(false);
						h.Scale(1,"width,nosw2");
						//h.GetYaxis()->SetTitle("Events/GeV");
						//h.SetMarkerStyle(20);
						//h.SetMarkerColor(mycol[n]);
						//h.SetFillColor(mycol[n]);
						h.SetLineColor(kBlack);
						//h.SetTitle(h.GetName());


						if(!this_run_comp){
							hcomp = (TH1D*)h.Clone(("comp_"+canvas_name).c_str());
							hcomp->Reset();
						}

						std::ostringstream out;
						out << std::setprecision(3) << total_events;
						std::string hmm = "\t";
						std::string tmp = h.GetName() +hmm+ out.str();


						hcomp->Add(&h);
						nc++;

						this_run_comp=true;

					}
				}


				//Ok to sort the histograms based on "size" will need a few tricks
				std::vector<double> integral_sorter;
				std::vector<TH1*> to_sort;
				std::vector<std::string> l_to_sort;

				for(auto &h : temp){
					std::string test = h.GetName();
					if(test.find(canvas_name)!=std::string::npos ){

						double total_events = h.GetSumOfWeights();
						h.Sumw2(false);
						h.Scale(1,"width,nosw2");
						h.GetYaxis()->SetTitle("Events/GeV");
						h.SetMarkerStyle(20);
						h.SetMarkerColor(mycol[n]);
						h.SetFillColor(mycol[n]);
						h.SetLineColor(kBlack);
						h.SetTitle(h.GetName());
						//h.Write();

						if(!this_run){
							hsum = (TH1D*)h.Clone(("sum_"+canvas_name).c_str());
							hsum->Reset();
						}

						std::ostringstream out;
						out <<std::fixed<< std::setprecision(0) << total_events;
						std::string hmm = " | ";
						std::string tmp = h.GetName()+hmm+ out.str();
						//std::string tmp = map_subchannel_plotnames.at(h.GetName()) +hmm+ out.str(); //fix this
						//legStack.AddEntry(&h, tmp.c_str() , "f");
						hsum->Add(&h);
						//hs->Add(&h);
						n++;

						this_run=true;

						to_sort.push_back(&h);
						l_to_sort.push_back(tmp);
						integral_sorter.push_back(total_events);

					}
				}
				//Sort!
				for (int i: SortIndexes(integral_sorter)) {
					hs->Add(to_sort.at(i));
					legStack.AddEntry(to_sort.at(i), l_to_sort.at(i).c_str(),"f");
				}

				legStack.AddEntry(hcomp, "Data", "fl");




				/****Not sure why but this next line seg faults...******
				 *	hs->GetYaxis()->SetTitle("Events/GeV");
				 ******************************************************/
				if(this_run && this_run_comp){
					double plot_pot=5e19;

					double title_size_ratio=0.1;
					double label_size_ratio=0.1;
					double title_offSet_ratioY = 0.3 ;
					double title_offSet_ratioX = 1.1;

					double title_size_upper=0.15;
					double label_size_upper=0.05;
					double title_offSet_upper = 1.45;


					Cstack->cd();
					TPad *pad0top = new TPad(("pad0top_"+canvas_name).c_str(), ("pad0top_"+canvas_name).c_str(), 0, 0.35, 1, 1.0);
					pad0top->SetBottomMargin(0); // Upper and lower plot are joined
					pad0top->Draw();             // Draw the upper pad: pad2top
					pad0top->cd();               // pad2top becomes the current pad
					hs->Draw();
					hcomp->SetLineColor(kBlack);
					hcomp->SetLineWidth(2);
					hcomp->SetFillColor(kBlack);
					hcomp->SetFillStyle(3244);

					hs->GetYaxis()->SetTitle("Events/GeV");
					hcomp->Draw("E2 same");
					TH1* hcompclone = (TH1*)hcomp->Clone("dangit");
					hcompclone->SetFillStyle(0);

					hcompclone->SetLineWidth(3);
					hcompclone->SetLineStyle(9);
					hcompclone->Draw("same hist");


					//hcomp->Draw("hist same");
					hs->SetMaximum(std::max(hs->GetMaximum(), hcomp->GetMaximum())*1.4);
					hs->SetMinimum(0.001);

					Cstack->Update();
					legStack.Draw();

					Cstack->cd();
					gStyle->SetOptStat(0);
					TPad *pad0bot = new TPad(("padbot_"+canvas_name).c_str(),("padbot_"+canvas_name).c_str(), 0, 0.05, 1, 0.35);
					pad0bot->SetTopMargin(0);
					pad0bot->SetBottomMargin(0.351);
					pad0bot->SetGridx(); // vertical grid
					pad0bot->Draw();
					pad0bot->cd();       // pad0bot becomes the current pad

					TH1* ratpre = (TH1*)hcomp->Clone(("ratio_"+canvas_name).c_str());
					ratpre->Divide(hsum);
					ratpre->SetStats(false);
					ratpre->Draw("hist");
					ratpre->SetFillColor(kWhite);
					ratpre->SetFillStyle(0);
					ratpre->SetLineWidth(2);

					gStyle->SetOptStat(0);
					TLine *line = new TLine(ratpre->GetXaxis()->GetXmin(),1.0,ratpre->GetXaxis()->GetXmax(),1.0 );
					line->Draw("same");
					ratpre->SetLineColor(kBlack);
					ratpre->SetTitle("");
					ratpre->GetYaxis()->SetTitle("Ratio");
					ratpre->GetXaxis()->SetTitleOffset(title_offSet_ratioX);
					ratpre->GetYaxis()->SetTitleOffset(title_offSet_ratioY);
					ratpre->SetMinimum(ratpre->GetMinimum()*0.97);
					ratpre->SetMaximum(ratpre->GetMaximum()*1.03);
					ratpre->GetYaxis()->SetTitleSize(title_size_ratio);
					ratpre->GetXaxis()->SetTitleSize(title_size_ratio);
					ratpre->GetYaxis()->SetLabelSize(label_size_ratio);
					ratpre->GetXaxis()->SetLabelSize(label_size_ratio);
					ratpre->GetXaxis()->SetTitle("Reconstructed Energy [GeV]");

					Cstack->Write(canvas_name.c_str() );
                    Cstack->SaveAs((tag+canvas_name+".pdf").c_str(),"pdf");

				}

			}
		}
	}

	f->Close();

	return 0;
}



int SBNspec::GetHistNumber(int f){
	//Get which histogram (in hist) that a full_vector entry comes from	
	int counter = 0;	
	for(int i=0; i<hist.size(); i++){
		if(f >= counter && f < counter + hist.at(i).GetNbinsX() ) return i;
		counter+=hist.at(i).GetNbinsX();
	}

	return -99;
}


int SBNspec::GetLocalBinNumber(double invar, int which_hist)
{
	int localbin = hist.at(which_hist).GetXaxis()->FindBin(invar);
	double bin = localbin-1;

	if(localbin==0 || localbin > hist.at(which_hist).GetNbinsX() ){bin = -99;}
	return bin;
}


int SBNspec::GetGlobalBinNumber(double invar, int which_hist)
{
	int localbin = hist.at(which_hist).GetXaxis()->FindBin(invar);
	double bin = localbin-1;

	for(int i=0; i<which_hist; i++){
		bin += hist.at(i).GetNbinsX();
	}

	if(localbin==0 || localbin > hist.at(which_hist).GetNbinsX() ){bin = -99;}
	return bin;
}
