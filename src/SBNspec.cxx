#include "SBNspec.h"
#include <cassert>
using namespace sbn;

SBNspec::SBNspec(std::string whichxml, int which_universe, bool isverbose, bool useuniverse) : SBNconfig(whichxml,isverbose, useuniverse){
	//Initialise all the things
	//for every multisim, create a vector of histograms, one for every subchannel we want
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

    has_been_scaled = false;
    m_bool_use_wire_bayes = false;
    this->CollapseVector();

}

SBNspec::SBNspec(std::string whichxml): SBNspec(whichxml,-1,true, true){}
SBNspec::SBNspec(std::string whichxml, int which_universe): SBNspec(whichxml,which_universe, true, true){}
SBNspec::SBNspec(std::string whichxml, int which_universe, bool isverbose): SBNspec(whichxml, which_universe, isverbose, true){}

SBNspec::SBNspec(std::string rootfile, std::string whichxml) : SBNspec(rootfile, whichxml, true){}
SBNspec::SBNspec(std::string rootfile, std::string whichxml, bool isverbose) : SBNconfig(whichxml, isverbose, true) {
	//Contruct from a prexisting histograms that exist in a rootfile
	TFile *f = new TFile(rootfile.c_str(),"read");

	//Loop over all filenames that should be there, and load up the histograms.
	int n=0;
	for(auto fn: fullnames){
		hist.push_back(*((TH1D*)f->Get(fn.c_str())));
		map_hist[fn] = n;
		n++;
	}

	has_been_scaled=false;
	m_bool_use_wire_bayes = false;
	f->Close();

	this->CalcFullVector();
	this->CalcErrorVector();
}

SBNspec::SBNspec(std::vector<double> input_full_vec, std::string whichxml) : SBNspec(input_full_vec, whichxml, false){ };
SBNspec::SBNspec(std::vector<double> input_full_vec, std::string whichxml, bool isverbose) : SBNspec(input_full_vec,whichxml,-1,isverbose){};
SBNspec::SBNspec(std::vector<double> input_full_vec, std::string whichxml, int universe, bool isverbose) : SBNspec(whichxml,universe,isverbose){

    for(int i=0; i< input_full_vec.size(); i++){

        int which_hist = GetHistNumber(i);

        int exact_bin = i;
        for(int b=0; b<which_hist; b++){
            exact_bin -= hist.at(b).GetNbinsX();
        }
        hist.at(which_hist).SetBinContent(exact_bin+1, input_full_vec.at(i));


    }
	m_bool_use_wire_bayes = false;
	this->CalcFullVector();
	this->CalcErrorVector();
}



SBNspec::SBNspec(std::vector<double> input_full_vec, std::vector<double> input_full_err, std::string whichxml) : SBNspec(input_full_vec, input_full_err, whichxml, false){ };
SBNspec::SBNspec(std::vector<double> input_full_vec, std::vector<double> input_full_err, std::string whichxml, bool isverbose) : SBNspec(input_full_vec, input_full_err, whichxml,-1,isverbose){};
SBNspec::SBNspec(std::vector<double> input_full_vec, std::vector<double> input_full_err, std::string whichxml, int universe, bool isverbose) : SBNspec(whichxml,universe,isverbose){

        if(input_full_vec.size() != input_full_err.size()){
                std::cerr << "SBNspec::SBNspec||\t Input full vector has different size from input full error vector: "<< input_full_vec.size() << ", " << input_full_err.size() << std::endl;
                exit(EXIT_FAILURE);
        }

        for(int i=0; i< input_full_vec.size(); i++){

                int which_hist = GetHistNumber(i);

                int exact_bin = i;
                for(int b=0; b<which_hist; b++){
                        exact_bin -= hist.at(b).GetNbinsX();
                }
                hist.at(which_hist).SetBinContent(exact_bin+1, input_full_vec.at(i));
                hist.at(which_hist).SetBinError(exact_bin+1, input_full_err.at(i));

        }

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
	this->CalcErrorVector();
	return 0;



int SBNspec::Add(SBNspec *in){
    //Addes all hists toGether
    if(xmlname != in->xmlname){ std::cout<<"ERROR: SBNspec::Add, trying to add differently configured SBNspecs!"<<std::endl; exit(EXIT_FAILURE);}

    for(int i=0; i< hist.size(); i++){
        hist[i].Add( &(in->hist[i]));
    }

	this->CollapseVector();
	this->CalcErrorVector();
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
			h.SetBinError(i, sqrt(val) );
		}
	}
}


//All scaling functions are quite self explanatory
int SBNspec::ScalePoisson(){
	TRandom3 *rangen = new TRandom3(0);
	for(auto &h: hist){
		for(int i=1; i<h.GetSize()-1; i++){
			h.SetBinContent(i, rangen->Poisson( h.GetBinContent(i)    ));
			h.SetBinError(i, sqrt( h.GetBinContent(i)));
		}
	}
	return 0;
}

int SBNspec::ScalePoisson(TRandom3* rangen){
	for(auto &h: hist){
		for(int i=1; i<h.GetSize()-1; i++){
			h.SetBinContent(i, rangen->Poisson( h.GetBinContent(i)    ));
			h.SetBinError(i, sqrt( h.GetBinContent(i)));
		}
	}
	return 0;
}



int SBNspec::ScaleRandom(){
    TRandom3 *rangen    = new TRandom3(0);

    for(auto& h: hist){
        h.Scale(rangen->Uniform(0,2));

	}
	this->CalcErrorVector();
	return 0;
}


int SBNspec::Scale(std::string name, TF1 * func){
	for(auto& h: hist){
		std::string test = h.GetName();
		if(test.find(name)!=std::string::npos ){
			for(int b=0; b<=h.GetNbinsX(); b++){
				double old_bin_error = h.GetBinError(b);
				//std::cout<<h.GetBinContent(b)<<" "<<h.GetBinCenter(b)<<" "<<func->Eval(h.GetBinCenter(b) )<<std::endl;
				h.SetBinContent(b, h.GetBinContent(b)*func->Eval(h.GetBinCenter(b) ) );
				h.SetBinError(b, old_bin_error*func->Eval(h.GetBinCenter(b) ) );
			}
		}

	}
	this->CollapseVector();
	this->CalcErrorVector();

	return 0;
}


int SBNspec::ScaleAll(double sc){
	for(auto& h: hist){
		h.Scale(sc);
	}
	this->CollapseVector();
	this->CalcErrorVector();
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
	this->CalcErrorVector();
	return 0;
}

int SBNspec::NormAll(double n){

	for(auto& h: hist) {
		h.Scale(n/h.GetSumOfWeights());
	}

 	this->CollapseVector();
        this->CalcErrorVector();
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

	this->CollapseVector();
        this->CalcErrorVector();
	return 0;

}


int SBNspec::RemoveMCError(){

  full_error.clear();
  full_error.resize(num_bins_total);

  int hoffset = 0;
  for(size_t hid=0; hid<hist.size(); ++hid) {
    auto& h =  hist[hid];
    for(int i = 1; i < (h.GetSize()-1); ++i){
      full_error[hoffset + i - 1] = 0;
      h.SetBinError(i,0.0);
    }
    hoffset += (h.GetSize()-2);
  }
  assert (hoffset == num_bins_total);
  return 0;
}

int SBNspec::Keep(std::string name, double val){
	for(auto& h: hist){
		std::string test = h.GetName();

		if(test.find(name)!=std::string::npos ) h.Scale(val);
		else h.Scale(0.0);
	}
	this->CollapseVector();
	this->CalcErrorVector();
	return 0;
}



int SBNspec::CalcFullVector(){
  full_vector.clear();
  full_error.clear();
  
  full_vector.resize(num_bins_total);
  full_error.resize(num_bins_total);

  int hoffset = 0;
  for(size_t hid=0; hid<hist.size(); ++hid) {
    const auto& h =  hist[hid];
    for(int i = 1; i < (h.GetSize()-1); ++i){
      full_vector[hoffset + i - 1] = h.GetBinContent(i);
      full_error[hoffset + i - 1] =  h.GetBinError(i);
    }
    hoffset += (h.GetSize()-2);
  }
    
  assert (hoffset == num_bins_total);

  if(m_bool_use_wire_bayes){
   
    // Loop over all channels: for each bin in this channel, get 
    // means:     vector<double> of GetBinContent()
    // sigmas2s:  vector<double> of GetBinError() //maybe square of this
    // pots: vector<double> of overall POT scale factor for that sample (for use when 0 bins)?

    // So this needs to be done for a collapsed bin. we only have assigned error for collapsed bins.. Hmm, but we sample from full. 
    // Probably need to not have this here then. 

    std::vector<double> means;
    std::vector<double> sigma2s;
    std::vector<double> pots;
    

    std::vector<double> ans = LEEana::wireWrapperBayes(means, sigma2s, pots);

    //the overall bin error covariance is then ans.back(), so errar is sqrt();
  }

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


int SBNspec::CalcErrorVector(){

	full_err_vector.clear();
	collapsed_err_vector.clear();

	int index=0;
	for(int im=0; im < num_modes ; im++){
	    for(int id = 0; id< num_detectors; id++){
		for(int ic= 0; ic< num_channels; ic++){
		    TH1D* hsum= nullptr;   //holds the summed histogram for one CHANNEL

		    for(int is=0; is< num_subchannels[ic]; is++){
			TH1D& h = hist[index];

			for(int ibin=0; ibin < h.GetNbinsX(); ibin++) full_err_vector.push_back(h.GetBinError(ibin+1));

			if(hsum == nullptr){
				hsum = (TH1D*)h.Clone();
				hsum->Reset();
			}	
			hsum->Add(&h);
			index++;
		    }

		   for(int ibin=0; ibin < hsum->GetNbinsX(); ibin++)
			collapsed_err_vector.push_back(hsum->GetBinError(ibin+1));

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
	std::cout.precision(2);
	std::cout.setf( std::ios::fixed, std:: ios::floatfield );
	for(double d: full_vector){
		std::cout<<d<<" ";
	}
	std::cout<<std::endl;
	return 0;
}

int SBNspec::PrintCollapsedVector(){
	std::cout.precision(2);
	std::cout.setf( std::ios::fixed, std:: ios::floatfield );
	for(double d: collapsed_vector){
		std::cout << d<<" ";
	}
	std::cout<<std::endl;
	return 0;
}


int SBNspec::WriteOut(std::string tag){
    //kWhite  = 0,   kBlack  = 1,   kGray    = 920,  kRed    = 632,  kGreen  = 416,
    //kBlue   = 600, kYellow = 400, kMagenta = 616,  kCyan   = 432,  kOrange = 800,
    //kSpring = 820, kTeal   = 840, kAzure   =  860, kViolet = 880,  kPink   = 900

	std::vector<int> mycol = {kAzure -9, kRed-7, kGreen-3, kBlue-6, kMagenta-3, kYellow-7,  kOrange-3, kBlue, kBlue+2,  kGreen+1,kBlue-7, kPink, kViolet, kCyan,kMagenta,kAzure};
	int colindex =0;


	bool gLEE_plot = true;  //control the style of the stacked histograms
	if(gLEE_plot){
		mycol.clear();

		std::map<std::string, int> color_channel_map;
		std::map<std::string, std::vector<double>> rgb_channel_map={
			{"NCDelta", {255./255.,255./255.,153./255.}},
			{"NCDeltaLEE", {0.97,0.75,0}},
			{"NCPi0Coh", {255./255,189./255.,189./255.}},
			{"NCPi0NotCoh", {1,0.4,0.4}},
			{"NCMultiPi0", {0.9,0.9,1.0}},
			{"CC1Pi0", {0.4,0.4,1.0}},
			{"BNBOther", {0.6,0.8,1.0}},
			{"Nue",{0.9,0.5,0.9}},
			{"Dirt", {0.6,0.4,0.2}},
			{"OTPCExtra", {0.2,0.5,0.2}},
			{"BNBext", {0.2,0.8,0.2}}
		};

		std::map<std::string, std::vector<double>>::iterator iter;
		std::map<std::string, int>::iterator iter_int;
		TColor* t_col = NULL;
		for(iter = rgb_channel_map.begin(); iter!= rgb_channel_map.end(); ++iter){
			int color_index = TColor::GetFreeColorIndex();
			t_col = new TColor(color_index, iter->second.at(0),iter->second.at(1),iter->second.at(2));	
			color_channel_map.insert({iter->first, t_col->GetNumber()});
		}
		
		for(int is = 0; is <subchannel_names[0].size(); is++){
			std::string isubchannel_name = subchannel_names[0][is];
			iter_int = color_channel_map.find(isubchannel_name);
			if(iter_int == color_channel_map.end()){
				std::cout << "Color of " << isubchannel_name << " is not defined, choose a random color" << std::endl;
				mycol.push_back(is);
			}
			else{
				std::cout << "Successfully find the color predefined for "<< isubchannel_name << std::endl;
				mycol.push_back(iter_int->second);
			}
		}
		
	} 


	TFile *f2 = new TFile((tag+".SBNspec.root").c_str(),"recreate");
    f2->cd();
    std::cout<<"This file has "<<hist.size()<<" histograms and "<<fullnames.size()<<" fullnames "<<std::endl;
    for(auto& h: hist){
        h.Write();
    }
    f2->Close();

	std::vector<TH1D> temp_hists = hist;

	TFile *f = new TFile(("SBNfit_spectrum_plots_"+tag+".root").c_str(),"RECREATE");
	f->cd();

    for(int im = 0; im <mode_names.size(); im++){
        for(int id = 0; id <detector_names.size(); id++){
            for(int ic = 0; ic <channel_names.size(); ic++){


                std::string canvas_name = mode_names.at(im)+"_"+detector_names.at(id)+"_"+channel_names.at(ic);


                bool this_run = false;

                TCanvas* Cstack= new TCanvas(canvas_name.c_str(),canvas_name.c_str());
                Cstack->SetFixedAspectRatio();

				Cstack->cd();
				THStack * hs = new THStack(canvas_name.c_str(),  canvas_name.c_str());
				TLegend* legStack = new TLegend(0.11,0.7,0.89,0.89);
				legStack->SetNColumns(2);
				legStack->SetLineWidth(0);
				legStack->SetLineColor(kWhite);
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
						//h.Scale(1,"width,nosw2");
						h.GetYaxis()->SetTitle("Events/GeV");
						h.SetMarkerStyle(20);
						h.SetMarkerColor(mycol[n]);
						h.SetFillColor(mycol[n]);
						if(gLEE_plot & (test.find("BNBext")!=std::string::npos)) h.SetFillStyle(3333);
						h.SetLineColor(kBlack);
						h.SetTitle(h.GetName());
						//h.Write();

                        if(!this_run){
                            hsum = (TH1D*)h.Clone(("sum_"+canvas_name).c_str());
                            hsum->Reset();
                        }

						std::ostringstream out;
						out <<std::fixed<< std::setprecision(2) << total_events;
						std::string hmm = " \t ";
						//std::string tmp =h.GetName() +hmm+ out.str();
						std::string tmp_name = h.GetName();
						std::string tmp = map_subchannel_plotnames[tmp_name] +hmm+ out.str();
						//legStack.AddEntry(&h, tmp.c_str() , "f");
						hsum->Add(&h);
						//hs->Add(&h);
						n++;

                        this_run=true;

                        to_sort.push_back(&h);
                        l_to_sort.push_back(tmp);
                        integral_sorter.push_back(total_events);

						if(gLEE_plot){
							hs->Add(&h);
							legStack->AddEntry(&h, tmp.c_str(),"f");
						}
					}
				}
	
				//Sort!
				if(!gLEE_plot){
					for (int i: SortIndexes(integral_sorter)) {
						hs->Add(to_sort.at(i));
						legStack->AddEntry(to_sort.at(i), l_to_sort.at(i).c_str(),"f");
					}
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

					//hs->GetYaxis()->SetTitle(("Events/"+channel_units.at(ic)).c_str());
					hs->GetYaxis()->SetTitle("Events");
					hs->GetXaxis()->SetTitle(channel_units.at(ic).c_str());


					//hcomp->Draw("hist same");
					hs->SetMaximum(hs->GetMaximum()*1.3);
					hs->SetMinimum(0.001);

					Cstack->Update();
					legStack->Draw();

                    f->cd();

                    Cstack->Write(canvas_name.c_str() );

				}
			}
		}
	}
    f->Close();

    return 0;
}

int SBNspec::CompareSBNspecsMark(TMatrixT<double> collapse_covar, SBNspec * compsec, std::string tag){
    //kWhite  = 0,   kBlack  = 1,   kGray    = 920,  kRed    = 632,  kGreen  = 416,
    //kBlue   = 600, kYellow = 400, kMagenta = 616,  kCyan   = 432,  kOrange = 800,
    //kSpring = 820, kTeal   = 840, kAzure   =  860, kViolet = 880,  kPink   = 900
    std::vector<int> mycol = {kAzure -9, kRed-7, kGreen-3, kBlue-6, kMagenta-3, kYellow-7,  kOrange-3, kBlue, kBlue+2,  kGreen+1,kBlue-7, kPink, kViolet, kCyan,kMagenta,kAzure};


    //first check the dimension of covariance matrix
    if(collapse_covar.GetNcols() != this->num_bins_total_compressed){
        std::cout << "ERROR: dimension of covariance matrix " << collapse_covar.GetName() << " is diffrent from num_bins_total_compressed" << std::endl;
        exit(EXIT_FAILURE);
    }


    bool gLEE_plot = true;  //control the style of the stacked histograms
    if(gLEE_plot){
        mycol.clear();

        std::map<std::string, int> color_channel_map;
        std::map<std::string, std::vector<double>> rgb_channel_map={
            {"NCDeltaRadOverlaySM", {255./255.,255./255.,153./255.}},
            {"NCDeltaRadOverlayLEE", {0.97,0.75,0}},
            {"NCPi0Coh", {255./255,189./255.,189./255.}},
            {"NCPi0NotCoh", {1,0.4,0.4}},
            {"NCMultiPi0", {0.9,0.9,1.0}},
            {"CC1Pi0", {0.4,0.4,1.0}},
            {"BNBOther", {0.6,0.8,1.0}},
            {"NueOverlays",{0.9,0.5,0.9}},
            {"Dirt", {0.6,0.4,0.2}},
            {"BNBext", {0.2,0.8,0.2}}
        };

        std::map<std::string, std::vector<double>>::iterator iter;
        std::map<std::string, int>::iterator iter_int;
        TColor* t_col = NULL;
        for(iter = rgb_channel_map.begin(); iter!= rgb_channel_map.end(); ++iter){
            int color_index = TColor::GetFreeColorIndex();
            t_col = new TColor(color_index, iter->second.at(0),iter->second.at(1),iter->second.at(2));	
            color_channel_map.insert({iter->first, t_col->GetNumber()});
        }

        for(int is = 0; is <subchannel_names[0].size(); is++){
            std::string isubchannel_name = subchannel_names[0][is];
            iter_int = color_channel_map.find(isubchannel_name);
            if(iter_int == color_channel_map.end()){
                std::cout << "Color of " << isubchannel_name << " is not defined, choose a random color" << std::endl;
                mycol.push_back(is);
            }
            else{
                //std::cout << "Successfully find the color predefined for "<< isubchannel_name << std::endl;
                mycol.push_back(iter_int->second);
            }
        }

    } 


    TFile *f = new TFile(("SBNfit_compare_plots_"+tag+".root").c_str(),"RECREATE");
    f->cd();


    std::vector<TH1D> temp = hist;
    std::vector<TH1D> temp_comp = compsec->hist;

    for(int k=0; k< fullnames.size(); k++){
        TCanvas *ctmp = new TCanvas((tag+"_"+std::to_string(k)+"_"+fullnames.at(k)).c_str(), (tag+" | "+std::to_string(k)+"_"+fullnames.at(k)).c_str(),1200,1200);
        ctmp->cd();

        if(map_hist.count(fullnames.at(k)) <=0 ){
            std::cout<<"WARNING:  hist "<<fullnames[k]<<" is not in map_hist"<<std::endl;
        }
        TH1D * h1 = (TH1D*) temp.at(map_hist[fullnames.at(k)]).Clone((std::to_string(k)+fullnames.at(k)+"_1").c_str());
        TH1D * h2 = (TH1D*) temp_comp.at(map_hist[fullnames.at(k)]).Clone((std::to_string(k)+fullnames.at(k)+"_2").c_str());

        //h1->Scale(1,"width");
        //h2->Scale(1,"width");

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

    for(int im = 0; im <mode_names.size(); im++){
        for(int id = 0; id <detector_names.size(); id++){
            for(int ic = 0; ic <channel_names.size(); ic++){

                std::string canvas_name = mode_names.at(im)+"_"+detector_names.at(id)+"_"+channel_names.at(ic);

                bool this_run = false;
                bool this_run_comp = false;

                TCanvas* Cstack= new TCanvas((tag+"_"+canvas_name).c_str(),(tag+" | "+canvas_name).c_str(),1200,1200);
                //Cstack->SetFixedAspectRatio();

                Cstack->cd();
                THStack * hs = new THStack(canvas_name.c_str(),  canvas_name.c_str());
                //TLegend legStack(0.50, 0.5, 0.89, 0.89);
                TLegend legStack(0.11, 0.69, 0.89, 0.89);
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

                        //h.Sumw2(false);
                        h.Scale(1,"width");
                        //h.Scale(1,"width,nosw2");
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
                        out << std::setprecision(2) << total_events;
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

                        //h.Sumw2(false);
                        //h.Scale(1,"width,nosw2");
                        h.Scale(1,"width");
                        h.GetYaxis()->SetTitle("Events/GeV");
                        h.SetMarkerStyle(20);
                        h.SetMarkerColor(mycol[n]);
                        h.SetFillColor(mycol[n]);
                        if(gLEE_plot & (test.find("BNBext")!=std::string::npos)) h.SetFillStyle(3333);
                        h.SetLineColor(kBlack);
                        std::string namo = (canvas_name+" | "+tag);
                        h.SetTitle(namo.c_str());
                        //h.Write();

                        if(!this_run){
                            hsum = (TH1D*)h.Clone(("sum_"+canvas_name).c_str());
                            hsum->Reset();
                        }

                        std::ostringstream out;
                        out <<std::fixed<< std::setprecision(2) << total_events;
                        std::string hmm = " | ";
                        std::string tmp = h.GetName()+hmm+ out.str();
                        std::string tmp_name = h.GetName();
                        //std::string tmp = map_subchannel_plotnames[tmp_name] +hmm+ out.str();
                        //legStack.AddEntry(&h, tmp.c_str() , "f");
                        hsum->Add(&h);
                        //hs->Add(&h);
                        n++;

                        this_run=true;

                        to_sort.push_back(&h);
                        l_to_sort.push_back(tmp);
                        integral_sorter.push_back(total_events);

                        if(gLEE_plot){
                            hs->Add(&h, "HIST");
                            legStack.AddEntry(&h, tmp.c_str(),"f");
                        }

                    }
                }

                if(!gLEE_plot){
                    //Sort!
                    for (int i: SortIndexes(integral_sorter)) {
                        hs->Add(to_sort.at(i), "HIST");
                        legStack.AddEntry(to_sort.at(i), l_to_sort.at(i).c_str(),"f");
                    }
                }




                //set the error of 'this' spec according to the covariance matrix
                for(int i=0; i<hsum->GetNbinsX(); i++){
                    double xbin_width = hsum->GetXaxis()->GetBinWidth(i+1);
                    double error = sqrt(pow(hsum->GetBinError(i+1), 2.0) + collapse_covar(error_bin+i, error_bin+i)/pow(xbin_width, 2.0));
                    //error = sqrt(collapse_covar(error_bin+i, error_bin+i))/xbin_width; //uncomment this is you want JUST systematic erros

                    //double error = hsum->GetBinError(i+1) + sqrt(collapse_covar(error_bin+i, error_bin+i));
                    //std::cout << collapse_covar(error_bin+i, error_bin+i) << std::endl;
                    std::cout << "previous error: "<< hsum->GetBinError(i+1) << ", later one: " << error << std::endl;
                    hsum->SetBinError(i+1, error);
                }
                error_bin +=hsum->GetNbinsX();

                legStack.AddEntry(hsum, Form("MC Stack | %.2f", hsum_sum), "fl");
                legStack.AddEntry(hcomp, Form("Compared Point | %.2f", hcomp_sum), "flp");

                /****Not sure why but this next line seg faults...******
                 *	hs->GetYaxis()->SetTitle("Events/GeV");
                 ******************************************************/
                if(this_run && this_run_comp){

                    double title_size_ratio=0.11;
                    double label_size_ratio=0.11;
                    double title_offSet_ratioY = 0.38;
                    double title_offSet_ratioX = 1.1;

                    double title_size_upper=0.048;
                    double label_size_upper=0.05;
                    double title_offSet_upper = 0.85;

                    Cstack->cd();
                    //gStyle->SetErrorX(0);
                    TPad *pad0top = new TPad(("pad0top_"+canvas_name).c_str(), ("pad0top_"+canvas_name).c_str(), 0, 0.35, 1, 1.0);
                    pad0top->SetBottomMargin(0); // Upper and lower plot are joined
                    pad0top->Draw();             // Draw the upper pad: pad2top
                    pad0top->cd();               // pad2top becomes the current pad
                    //gStyle->SetHatchesLineWidth(1);

                    hs->Draw();
                    //draw error bar for 'this' SBNspec
                    hsum->SetFillColor(kBlack);
                    hsum->SetFillStyle(3354);
                    hsum->SetLineWidth(3);
                    hsum->SetMarkerSize(0);
                    hsum->Draw("E2 same");
                    hsum->DrawCopy("E2 same");
                    hs->GetYaxis()->SetTitle("Events/BinWidth");
                    hs->GetYaxis()->SetTitleSize(title_size_upper);
                    hs->GetYaxis()->SetLabelSize(label_size_upper);
                    hs->GetYaxis()->SetTitleOffset(title_offSet_upper*1.2);
                    //hs->GetYaxis()->SetTitle("Events");

                    //draw rectangular and dashed line for compared spectrum
                    /*hcomp->SetLineColor(kBlack);
                      hcomp->SetLineWidth(2);
                      hcomp->SetFillColor(kBlack);
                      hcomp->SetFillStyle(3244);

                      hcomp->Draw("E2 same");

                      TH1* hcompclone = (TH1*)hcomp->Clone("dangit");
                      hcompclone->SetFillStyle(0);

                      hcompclone->SetLineWidth(3);
                      hcompclone->SetLineStyle(9);
                      hcompclone->Draw("same hist");
                      */

                    gStyle->SetEndErrorSize(5);
                    hcomp->SetLineWidth(2);
                    hcomp->SetMarkerStyle(20);
                    hcomp->SetMarkerSize(1.5);
                    hcomp->SetLineColor(kBlack);
                    hcomp->Draw("E1P same");

                    //hcomp->Draw("hist same");
                    hs->SetMaximum(std::max(hs->GetMaximum(), hcomp->GetMaximum())*1.6);
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


                    TH1* ratpre = (TH1*)hcomp->Clone(("ratio_"+canvas_name).c_str());
                    ratpre->Divide(hsum);
                    ratpre->SetStats(false);
                    //to draw the 1 sigma error band on the ratio plot
                    TH1* h_err = (TH1*)hsum->Clone("error_band");
                    h_err->Divide(hsum);
                    for(int i=0; i<h_err->GetNbinsX(); i++){
                        h_err->SetBinError(i+1, hsum->GetBinError(i+1)/hsum->GetBinContent(i+1));
                        ratpre->SetBinError(i+1, hcomp->GetBinError(i+1)/hsum->GetBinContent(i+1));
                    }


                    ratpre->Draw("E1");
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
                    ratpre->GetYaxis()->SetTitle("Data/Prediction");
                    ratpre->GetXaxis()->SetTitleOffset(title_offSet_ratioX);
                    ratpre->GetYaxis()->SetTitleOffset(title_offSet_ratioY);
                    ratpre->SetMinimum(std::min(0.0, ratpre->GetMinimum())*0.8);
                    //ratpre->SetMinimum(ratpre->GetMinimum()*0.97);
                    ratpre->SetMaximum(std::max(2.0, ratpre->GetMaximum())*1.2);
                    //ratpre->SetMaximum(ratpre->GetMaximum()*1.03);
                    ratpre->GetYaxis()->SetNdivisions(505, kTRUE);   //change the label division in y axis
                    ratpre->GetYaxis()->SetTitleSize(title_size_ratio);
                    ratpre->GetXaxis()->SetTitleSize(title_size_ratio);
                    ratpre->GetYaxis()->SetLabelSize(label_size_ratio);
                    ratpre->GetXaxis()->SetLabelSize(label_size_ratio);
                    //ratpre->GetXaxis()->SetTitle("Reconstructed Energy [GeV]");
                    ratpre->GetXaxis()->SetTitle(channel_units.at(ic).c_str()); 
                    Cstack->SaveAs((canvas_name+"_"+tag+".pdf").c_str(),"pdf");
                    Cstack->Write(canvas_name.c_str() );

                }

            }
        }
    }

    f->Close();

    return 0;
}
//compare 'this' spec with the compsec with tag
//collapsed covariance matrix 'collapse_covar' gives the error of 'this' spec
int SBNspec::CompareSBNspecs(TMatrixT<double> collapse_covar, SBNspec * compsec, bool inbool, std::string tag){
	//kWhite  = 0,   kBlack  = 1,   kGray    = 920,  kRed    = 632,  kGreen  = 416,
	//kBlue   = 600, kYellow = 400, kMagenta = 616,  kCyan   = 432,  kOrange = 800,
	//kSpring = 820, kTeal   = 840, kAzure   =  860, kViolet = 880,  kPink   = 900
	std::vector<int> mycol = {kAzure -9, kRed-7, kGreen-3, kBlue-6, kMagenta-3, kYellow-7,  kOrange-3, kBlue, kBlue+2,  kGreen+1,kBlue-7, kPink, kViolet, kCyan,kMagenta,kAzure};

	bool is_compare_data=inbool;//whether the compared spec should be treated as real data or not

	//first check the dimension of covariance matrix
	if(collapse_covar.GetNcols() != this->num_bins_total_compressed){
		std::cout << "ERROR: dimension of covariance matrix " << collapse_covar.GetName() << " is diffrent from num_bins_total_compressed" << std::endl;
		exit(EXIT_FAILURE);
	}


	bool gLEE_plot = true;  //control the style of the stacked histograms
	if(gLEE_plot){
		mycol.clear();

		std::map<std::string, int> color_channel_map;
		std::map<std::string, std::vector<double>> rgb_channel_map={
			{"NCDelta", {255./255.,255./255.,153./255.}},
			{"NCDeltaLEE", {0.97,0.75,0}},
			{"NCPi0Coh", {255./255,189./255.,189./255.}},
			{"NCPi0NotCoh", {1,0.4,0.4}},
			{"NCMultiPi0", {0.9,0.9,1.0}},
			{"CC1Pi0", {0.4,0.4,1.0}},
			{"BNBOther", {0.6,0.8,1.0}},
			{"Nue",{0.9,0.5,0.9}},
			{"Dirt", {0.6,0.4,0.2}},
			{"OTPCExtra", {0.2,0.5,0.2}},
			{"BNBext", {0.2,0.8,0.2}}
		};

		std::map<std::string, std::vector<double>>::iterator iter;
		std::map<std::string, int>::iterator iter_int;
		TColor* t_col = NULL;
		for(iter = rgb_channel_map.begin(); iter!= rgb_channel_map.end(); ++iter){
			int color_index = TColor::GetFreeColorIndex();
			t_col = new TColor(color_index, iter->second.at(0),iter->second.at(1),iter->second.at(2));	
			color_channel_map.insert({iter->first, t_col->GetNumber()});
		}
		
		for(int is = 0; is <subchannel_names[0].size(); is++){
			std::string isubchannel_name = subchannel_names[0][is];
			iter_int = color_channel_map.find(isubchannel_name);
			if(iter_int == color_channel_map.end()){
				std::cout << "Color of " << isubchannel_name << " is not defined, choose a random color" << std::endl;
				mycol.push_back(is);
			}
			else{
				//std::cout << "Successfully find the color predefined for "<< isubchannel_name << std::endl;
				mycol.push_back(iter_int->second);
			}
		}
		
	} 

    std::vector<TH1D> temp = hist;
    std::vector<TH1D> temp_comp = compsec->hist;


	TFile *f = new TFile(("SBNfit_compare_plots_"+tag+".root").c_str(),"RECREATE");
	f->cd();


	for(int k=0; k< fullnames.size(); k++){
		TCanvas *ctmp = new TCanvas((tag+"_"+std::to_string(k)+"_"+fullnames.at(k)).c_str(), (std::to_string(k)+"_"+fullnames.at(k)).c_str(),1200,1200);
		ctmp->cd();

        if(map_hist.count(fullnames.at(k)) <=0 ){
            std::cout<<"WARNING:  hist "<<fullnames[k]<<" is not in map_hist"<<std::endl;
        }
        TH1D * h1 = (TH1D*) temp.at(map_hist[fullnames.at(k)]).Clone((std::to_string(k)+fullnames.at(k)+"_1").c_str());
        TH1D * h2 = (TH1D*) temp_comp.at(map_hist[fullnames.at(k)]).Clone((std::to_string(k)+fullnames.at(k)+"_2").c_str());

		//h1->Scale(1,"width");
		//h2->Scale(1,"width");

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
                TLegend legStack(0.65,0.25,0.89,0.89);
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

	

	int error_bin=0;

         for(int im = 0; im <mode_names.size(); im++){
                for(int id = 0; id <detector_names.size(); id++){
                        for(int ic = 0; ic <channel_names.size(); ic++){
		
				std::string canvas_name = mode_names.at(im)+"_"+detector_names.at(id)+"_"+channel_names.at(ic);

                        to_sort.push_back(&h);
                        l_to_sort.push_back(tmp);
                        integral_sorter.push_back(total_events);

				TCanvas* Cstack= new TCanvas((tag+"_"+canvas_name).c_str(),(tag+" | "+canvas_name).c_str(),1450,1200);
				//Cstack->SetFixedAspectRatio();

				Cstack->cd();
				THStack * hs = new THStack(canvas_name.c_str(),  canvas_name.c_str());
				//TLegend legStack(0.50, 0.5, 0.89, 0.89);
				TLegend legStack(0.11, 0.58, 0.89, 0.89);
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

						//h.Sumw2(false);
						h.Scale(1,"width");
						//h.Scale(1,"width,nosw2");
						//h.GetYaxis()->SetTitle("Events/GeV");
						//h.SetMarkerStyle(20);
						//h.SetMarkerColor(mycol[n]);
						//h.SetFillColor(mycol[n]);
						h.SetLineColor(kBlack);
						//h.SetTitle(h.GetName());

                    double title_size_upper=0.15;
                    double label_size_upper=0.05;
                    double title_offSet_upper = 1.45;


						std::ostringstream out;
						out << std::setprecision(2) << total_events;
						std::string hmm = "\t";
						std::string tmp = h.GetName() +hmm+ out.str();

                    hs->GetYaxis()->SetTitle("Events/GeV");
                    hcomp->Draw("E2 same");
                    TH1* hcompclone = (TH1*)hcomp->Clone("dangit");
                    hcompclone->SetFillStyle(0);

                    hcompclone->SetLineWidth(3);
                    hcompclone->SetLineStyle(9);
                    hcompclone->Draw("same hist");


                    //hcomp->Draw("hist same");
                    hs->SetMaximum(std::max(hs->GetMaximum(), hcomp->GetMaximum())*1.1);
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

						double total_events = h.GetSumOfWeights();
						hsum_sum += total_events;

						//h.Sumw2(false);
						//h.Scale(1,"width,nosw2");
						h.Scale(1,"width");
						h.GetYaxis()->SetTitle("Events/GeV");
						h.SetMarkerStyle(20);
						h.SetMarkerColor(mycol[n]);
						h.SetFillColor(mycol[n]);
						if(gLEE_plot & (test.find("BNBext")!=std::string::npos)) h.SetFillStyle(3333);
						h.SetLineColor(kBlack);
						h.SetTitle(h.GetName());
						//h.Write();

						if(!this_run){
							hsum = (TH1D*)h.Clone(("sum_"+canvas_name).c_str());
							hsum->Reset();
						}

						std::ostringstream out;
						out <<std::fixed<< std::setprecision(2) << total_events;
						std::string hmm = " | ";
						//std::string tmp = h.GetName()+hmm+ out.str();
						std::string tmp_name = h.GetName();
						std::string tmp = map_subchannel_plotnames[tmp_name] +hmm+ out.str();
						//legStack.AddEntry(&h, tmp.c_str() , "f");
						hsum->Add(&h);
						//hs->Add(&h);
						n++;

						this_run=true;

						to_sort.push_back(&h);
						l_to_sort.push_back(tmp);
						integral_sorter.push_back(total_events);

						if(gLEE_plot){
							hs->Add(&h, "HIST");
							legStack.AddEntry(&h, tmp.c_str(),"f");
						}

					}
				}

				if(!gLEE_plot){
				//Sort!
				for (int i: SortIndexes(integral_sorter)) {
					hs->Add(to_sort.at(i), "HIST");
					legStack.AddEntry(to_sort.at(i), l_to_sort.at(i).c_str(),"f");
				}
				}




				//set the error of 'this' spec according to the covariance matrix
				for(int i=0; i<hsum->GetNbinsX(); i++){
					double xbin_width = hsum->GetXaxis()->GetBinWidth(i+1);
					double error = sqrt(collapse_covar(error_bin+i, error_bin+i))/xbin_width;
					//std::cout << collapse_covar(error_bin+i, error_bin+i) << std::endl;
					std::cout << "MC intrinsic error: "<< hsum->GetBinError(i+1) << ", overall syst error: " << error << std::endl;
					hsum->SetBinError(i+1, error);
	
					//treat compared spec as real data
					if(is_compare_data){
					   double comp_error = sqrt(hcomp->GetBinContent(i+1)/xbin_width);
					   hcomp->SetBinError(i+1, comp_error);
					}
				}
				error_bin +=hsum->GetNbinsX();

				legStack.AddEntry(hsum, Form("MC Stack | %.2f", hsum_sum), "fl");
				legStack.AddEntry(hcomp, Form("Data | %.2f", hcomp_sum), "flp");

				if(this_run && this_run_comp){
					double plot_pot=5e19;

					double title_size_ratio=0.11;
                                        double label_size_ratio=0.11;
                                        double title_offSet_ratioY = 0.45;
                                        double title_offSet_ratioX = 1.1;

                                        double title_size_upper=0.048;
                                        double label_size_upper=0.05;
                                        double title_offSet_upper = 0.85;

					Cstack->cd();
					//gStyle->SetErrorX(0);
					TPad *pad0top = new TPad(("pad0top_"+canvas_name).c_str(), ("pad0top_"+canvas_name).c_str(), 0, 0.35, 1, 1.0);
					pad0top->SetBottomMargin(0); // Upper and lower plot are joined
					pad0top->Draw();             // Draw the upper pad: pad2top
					pad0top->cd();               // pad2top becomes the current pad
					//gStyle->SetHatchesLineWidth(1);
					hs->Draw();
					//draw error bar for 'this' SBNspec
					hsum->SetFillColor(kBlack);
					hsum->SetFillStyle(3354);
					hsum->SetMarkerSize(0);
					hsum->SetLineWidth(3);
					hsum->Draw("E2 same");
					//hsum->DrawClone("E2 same");
					hs->GetYaxis()->SetTitle(("Events/"+channel_units.at(ic)).c_str());
					hs->GetYaxis()->SetTitleSize(title_size_upper);
                                        hs->GetYaxis()->SetLabelSize(label_size_upper);
                                        hs->GetYaxis()->SetTitleOffset(title_offSet_upper*1.2);
					//hs->GetYaxis()->SetTitle("Events");
					
					//draw rectangular and dashed line for compared spectrum
					/*hcomp->SetLineColor(kBlack);
					hcomp->SetLineWidth(2);
					hcomp->SetFillColor(kBlack);
					hcomp->SetFillStyle(3244);

					hcomp->Draw("E2 same");

					TH1* hcompclone = (TH1*)hcomp->Clone("dangit");
					hcompclone->SetFillStyle(0);

					hcompclone->SetLineWidth(3);
					hcompclone->SetLineStyle(9);
					hcompclone->Draw("same hist");
					*/

					gStyle->SetEndErrorSize(5);
					hcomp->SetMarkerStyle(20);
					hcomp->SetMarkerColor(kBlack);
					hcomp->SetMarkerSize(1.5);
					hcomp->SetLineWidth(3);
					hcomp->SetLineColor(kBlack);
					hcomp->Draw("E1P same");

					//hcomp->Draw("hist same");
					hs->SetMaximum(std::max(hs->GetMaximum(), hcomp->GetMaximum())*2.5);
					//hs->SetMaximum(std::max(hs->GetMaximum(), hcomp->GetMaximum())*1.1);
					hs->SetMinimum(0.001);

				        Cstack->Update();
					legStack.Draw();
					TText t_text(0.7*(hs->GetXaxis()->GetXmax()), 1.4*hs->GetMaximum(), (channel_names.at(ic)).c_str());
					t_text.SetTextSize(0.08);
					t_text.Draw();

					Cstack->cd();
					//gStyle->SetOptStat(0);
					TPad *pad0bot = new TPad(("padbot_"+canvas_name).c_str(),("padbot_"+canvas_name).c_str(), 0, 0.05, 1, 0.35);
					pad0bot->SetTopMargin(0);
					pad0bot->SetBottomMargin(0.351);
					pad0bot->SetGridx(); // vertical grid
					pad0bot->Draw();
					pad0bot->cd();       // pad0bot becomes the current pad


					TH1* ratpre = (TH1*)hcomp->Clone(("ratio_"+canvas_name).c_str());
					ratpre->Divide(hsum);
					ratpre->SetStats(false);
					//to draw the 1 sigma error band on the ratio plot
					TH1* h_err = (TH1*)hsum->Clone("error_band");
                                        h_err->Divide(hsum);
					for(int i=0; i<h_err->GetNbinsX(); i++){
                                                h_err->SetBinError(i+1, hsum->GetBinError(i+1)/hsum->GetBinContent(i+1));
						ratpre->SetBinError(i+1, hcomp->GetBinError(i+1)/hsum->GetBinContent(i+1));
                                        }


					ratpre->Draw("E1");
					//ratpre->Draw("hist");
					
					ratpre->SetFillColor(kWhite);
					//ratpre->SetFillStyle(0);
					ratpre->SetLineWidth(2);
					h_err->Draw("E2same");
					Cstack->Update();
					gStyle->SetOptStat(0);
					TLine *line = new TLine(ratpre->GetXaxis()->GetXmin(),1.0,ratpre->GetXaxis()->GetXmax(),1.0 );
					line->Draw("same");
					ratpre->SetLineColor(kBlack);
					ratpre->SetTitle("");
					ratpre->GetYaxis()->SetTitle("Data/Prediction");
					ratpre->GetXaxis()->SetTitleOffset(title_offSet_ratioX);
					ratpre->GetYaxis()->SetTitleOffset(title_offSet_ratioY);
					ratpre->SetMinimum(std::min(0.5, ratpre->GetMinimum())*0.8);
					//ratpre->SetMinimum(ratpre->GetMinimum()*0.97);
					ratpre->SetMaximum(std::max(1.5, ratpre->GetMaximum())*1.2);
					//ratpre->SetMaximum(ratpre->GetMaximum()*1.03);
					ratpre->GetYaxis()->SetNdivisions(505, kTRUE);   //change the label division in y axis
					ratpre->GetYaxis()->SetTitleSize(title_size_ratio);
					ratpre->GetXaxis()->SetTitleSize(title_size_ratio);
					ratpre->GetYaxis()->SetLabelSize(label_size_ratio);
					ratpre->GetXaxis()->SetLabelSize(label_size_ratio);
					//ratpre->GetXaxis()->SetTitle("Reconstructed Energy [GeV]");
					ratpre->GetXaxis()->SetTitle(channel_units.at(ic).c_str()); 
					Cstack->Update();
                    			Cstack->SaveAs((tag+"_"+canvas_name+".pdf").c_str(),"pdf");
					Cstack->Write(canvas_name.c_str() );
					//Cstack->Print((tag+canvas_name+".pdf").c_str());

                }

		/*		delete Cstack;
				delete hcomp;
				delete hsum;
				delete hs;
		*/
			}
		}
	}

    f->Close();

    return 0;
}




int SBNspec::CompareSBNspecs(SBNspec * compsec, std::string tag){
	TMatrixT<double> covar(this->num_bins_total_compressed, this->num_bins_total_compressed);
        covar.Zero();
	this->CalcErrorVector();
	// add MC intrinsic error
	for(int i=0 ; i< this->num_bins_total_compressed; i++) covar(i,i) = pow(collapsed_err_vector[i], 2.0);

	this->CompareSBNspecs(covar, compsec, true, tag);
}

int SBNspec::CompareSBNspecs(SBNspec * compsec, bool inbool, std::string tag){
	TMatrixT<double> covar(this->num_bins_total_compressed, this->num_bins_total_compressed);
	covar.Zero();
	this->CalcErrorVector();
	// add MC intrinsic error
	for(int i=0 ; i< this->num_bins_total_compressed; i++) covar(i,i) = pow(collapsed_err_vector[i], 2.0);

	return this->CompareSBNspecs(covar,compsec, inbool, tag);
}


int SBNspec::CompareSBNspecs(TMatrixT<double> collapse_covar, SBNspec * compsec, std::string tag){
	return this->CompareSBNspecs(collapse_covar,compsec, true, tag);
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


std::vector<int> SBNspec::GetIndiciesFromSubchannel(std::string const & subchannel){
    std::vector<int> ans;
    if(map_hist.count(subchannel)!=1){
        std::cout<<"Error! You asked for a subchannel that doesnt exist in xml: "<<subchannel<<std::endl;
        exit(EXIT_FAILURE);
    }
    int hist_num =this->map_hist[subchannel];

    this->CalcFullVector();
    for(int i=0; i<this->full_vector.size(); i++){
        if( GetHistNumber(i)==hist_num) ans.push_back(i); 
    }

    return ans;
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

int SBNspec::GetGlobalBinNumber(int local_bin, std::string histname){
	
	int which_hist = map_hist.at(histname);
	int bin = local_bin -1;

    for(int i=0; i<which_hist; i++){
        bin += hist.at(i).GetNbinsX();
    }

    return bin;
}


std::map<int,std::vector<int>> SBNspec::GetCollapsedChannelIndicies(){

    std::map<int,std::vector<int>> ans;
    int nc =0;

	for(int im = 0; im < num_modes; im++){
		for(int id =0; id < num_detectors; id++){
            int edge = id*num_bins_detector_block + num_bins_mode_block*im; // This is the starting index for this detector
			for(int ic = 0; ic < num_channels; ic++){
				int corner=edge;
                int end = corner+num_bins.at(ic);
                std::cout<<"This channel "<<ic<<" starts at "<<corner<<" and ends at "<<end<<std::endl;
                edge += num_bins.at(ic);
                std::vector<int> tmp = {corner,end};
                ans[nc] = tmp;
                nc++;
			}
		}
	}
	return ans ;
}


std::vector<TH1D> SBNspec::GetBlankChannelHists(){

    std::vector<TH1D> ans;
    std::vector<int> vchan;
    int nc =0;

	for(int im = 0; im < num_modes; im++){
		for(int id =0; id < num_detectors; id++){
			for(int ic = 0; ic < num_channels; ic++){
                //channel num nc
                vchan.push_back(nc);
                nc+= num_subchannels.at(ic);
			}
		}
	}
    std::cout<<"We have : "<<vchan.size()<< " unique channels"<<std::endl;
    for(int i=0; i< vchan.size();i++){
        std::cout<<vchan[i]<<" "<<std::endl;
        ans.push_back(hist[vchan[i]]);
        ans.back().Reset();
    }
	return ans ;
}

