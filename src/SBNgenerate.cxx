#include "SBNgenerate.h"
//#include "MCEventWeight.h"

using namespace sbn;



SBNgenerate::SBNgenerate(std::string xmlname) {
    NeutrinoModel nullModel(0,0,0,true);
    SBNgenerate(xmlname, nullModel);
}

SBNgenerate::SBNgenerate(std::string xmlname, NeutrinoModel inModel ) : SBNconfig(xmlname), nu_model(inModel) {




    TRandom3 *rangen = new TRandom3(0);
    
    bool m_use_eventweight = false;


    //	gROOT->ProcessLine("#include <map>");
    //	gROOT->ProcessLine("#include <vector>");
    //	gROOT->ProcessLine("#include <string>");

    //gSystem->Load("../src/libranch_weightsMapDict.so");

    //	std::string dict_location = "../../dict/AutoDict_map_string__vector_double____cxx.so";
    //	gSystem->Load(  (dict_location).c_str());

    //	gSystem->Load("/uboone/app/users/markrl/sbnfit/whipping_star/src/mdict_h.so");
    //
    //std::cout<<"Trying to load dictionary: "<<dict_location<<std::endl;
    //

    std::map<std::string, int> parameter_sims;

    //Initialise the central value SBNspec.
    SBNspec tm(xmlname,-1,false);
    spec_central_value = tm;
    spec_osc_sin  = tm;
    spec_osc_sinsq = tm;

    int num_files = montecarlo_file.size();
    montecarlo_additional_weight.resize(num_files,1.0);
    montecarlo_additional_weight_formulas.resize(num_files);   


    for(auto &fn: montecarlo_file){
        files.push_back(new TFile(fn.c_str()));
        if(files.back()->IsZombie() || !files.back()->IsOpen()){
            std::cout<<"SBNgenerate || ERROR! Failed to oben the file "<<fn<<std::endl;
            exit(EXIT_FAILURE);
        }
    }


    for(int i=0; i<montecarlo_name.size(); i++){
        std::cout<<"Getting TTree "<<montecarlo_name[i]<<" from file "<<montecarlo_file[i]<<std::endl;
        trees.push_back((TTree*)files.at(i)->Get(montecarlo_name.at(i).c_str()) );
        std::cout<<"--TTree has "<<trees.back()->GetEntries()<<" entries. "<<std::endl;
    }

    for(int i=0; i<montecarlo_file.size(); i++){
	const auto& fn = montecarlo_file.at(i);
	auto montecarlo_file_friend_treename_iter = montecarlo_file_friend_treename_map.find(fn);
	if (montecarlo_file_friend_treename_iter != montecarlo_file_friend_treename_map.end()) {
            std::cout<<" Detected friend trees "<<std::endl;

            auto montecarlo_file_friend_iter = montecarlo_file_friend_map.find(fn);
            if (montecarlo_file_friend_iter == montecarlo_file_friend_map.end()) {
                std::stringstream ss;
                ss << "Looked for filename=" << fn << " in fnmontecarlo_file_friend_iter, but could not be found... bad config?" << std::endl;
                throw std::runtime_error(ss.str());
            }

            for(int k=0; k < (*montecarlo_file_friend_iter).second.size(); k++){

                std::string treefriendname = (*montecarlo_file_friend_treename_iter).second.at(k);
                std::string treefriendfile = (*montecarlo_file_friend_iter).second.at(k);

                std::cout <<" Adding a friend tree:  " <<treefriendname<<" from file: "<< treefriendfile <<std::endl;

                trees[i]->AddFriend(treefriendname.c_str(),treefriendfile.c_str());
            }
        }

    }

    std::vector<int> nentries;
    for(auto &t: trees){
        nentries.push_back(t->GetEntries());
    }

    f_weights.resize(num_files,nullptr);

    for(int i=0; i< num_files; i++){

        double pot_scale = 1.0;
        if(montecarlo_pot[i]!=-1){
            pot_scale = this->plot_pot/montecarlo_pot[i];
        }

        montecarlo_scale[i] = montecarlo_scale[i]*pot_scale;

        std::cout << " TFile::Open() file=" << files[i]->GetName() << " @" << files[i] << std::endl;
        std::cout << " Has POT " <<montecarlo_pot[i] <<" and "<<nentries[i] <<" entries "<<std::endl;


        //if(m_use_eventweight)  trees[i]->SetBranchAddress("eventweights", &(f_weights[i]) );
	    if(m_use_eventweight)  trees.at(i)->SetBranchAddress(montecarlo_eventweight_branch_names[i].c_str(), &(f_weights[i]));
        //delete f_weights->at(i);	f_weights->at(i) = 0;
        
        for(int k=0; k<branch_variables.at(i).size(); k++){
	    const auto branch_variable = branch_variables.at(i).at(k);
            std::cout<<"Setting Branch: "<<branch_variable->name<<std::endl;
            //trees.at(i)->SetBranchAddress( branch_variables.at(i).at(k)->name.c_str(), branch_variables.at(i).at(k)->GetValue() );
	    branch_variable->branch_formula =  new TTreeFormula(("branch_form"+std::to_string(i)).c_str(), branch_variable->name.c_str(), trees[i]);

            if(branch_variable->GetOscillate()){
                std::cout<<"Setting true branch variables"<<std::endl;
                //trees.at(i)->SetBranchAddress( branch_variable->true_param_name.c_str(), branch_variable->GetTrueValue() );
                //trees.at(i)->SetBranchAddress( branch_variable->true_L_name.c_str(), branch_variable->GetTrueL() );

		//update to use TTreeFormula! 
		branch_variable->branch_true_value_formula = new TTreeFormula(("branch_true_value_form"+std::to_string(i)).c_str(),branch_variable->true_param_name.c_str(), trees[i]);
		branch_variable->branch_true_L_formula = new TTreeFormula(("branch_true_L_form"+std::to_string(i)).c_str(), branch_variable->true_L_name.c_str(), trees[i]);
            }
        }

        if(montecarlo_additional_weight_bool[i]){
            //we have an additional weight we want to apply at run time, otherwise its just set at 1. 
	    std::cout<<"Setting Additional weight of : "<< montecarlo_additional_weight_names[i].c_str()<<std::endl;
            //trees[i]->SetBranchAddress(montecarlo_additional_weight_names[i].c_str(), &montecarlo_additional_weight[i]); 
	    montecarlo_additional_weight_formulas[i] =  new TTreeFormula(("a_w"+std::to_string(i)).c_str(),montecarlo_additional_weight_names[i].c_str(),trees[i]);
        }


    }
    std::string     bnbcorrection_str = "bnbcorrection_FluxHist";


    std::cout<<"SBNgenerate::SBNgenerate\t|| -------------------------------------------------------------\n";
    std::cout<<"SBNgenerate::SBNgenerate\t|| -------------------------------------------------------------\n";
    std::vector<double> base_vec (spec_central_value.num_bins_total,0.0);

    for(int j=0;j<num_files;j++){


        for(int i=0; i< std::min(  montecarlo_maxevents.at(j)  ,nentries.at(j)); i++){
            trees.at(j)->GetEntry(i);
            std::map<std::string,std::vector<eweight_type>>* thisfWeight;
            if(m_use_eventweight) thisfWeight = f_weights[j];

            if(i%100==0) std::cout<<"SBNgenerate::SBNgenerate\t|| On event: "<<i<<" of "<<nentries[j]<<" from File: "<<montecarlo_file[j]<<std::endl;

            double global_weight = 1.0;
	    if( montecarlo_additional_weight_bool[j]){
		    montecarlo_additional_weight_formulas[j]->GetNdata();
		    global_weight = montecarlo_additional_weight_formulas[j]->EvalInstance();
            };//this will be 1.0 unless specified
            global_weight = global_weight*montecarlo_scale[j];

            if(m_use_eventweight){
                if(thisfWeight->count("bnbcorrection_FluxHist")>0){
                     global_weight = global_weight*thisfWeight->at("bnbcorrection_FluxHist").front();
                }
            }


            if(std::isinf(global_weight) || global_weight != global_weight){
                std::cout<<"SBNgenerate::SBNgenerate\t|| ERROR  error @ "<<i<<" in File "<<montecarlo_file.at(j)<<" as its either inf/nan: "<<global_weight<<std::endl;
                exit(EXIT_FAILURE);
            }
	

            if( this->EventSelection(j) ){
                
                for(int t=0; t<branch_variables[j].size();t++){
                    //std::cout<<"Starting branch : "<<branch_variables.at(j).at(t)->name<<" "<<branch_variables.at(j).at(t)->associated_hist<<std::endl;
                    //Need the histogram index, the value, the global bin...

		            const auto branch_variable = branch_variables[j][t];
                    int ih = spec_central_value.map_hist.at(branch_variable->associated_hist);
		            branch_variable->GetFormula()->GetNdata();
		            double reco_var = branch_variable->GetFormula()->EvalInstance();
                    //double reco_var = *(static_cast<double*>(branch_variables[j][t]->GetValue()));
                    int reco_bin = spec_central_value.GetGlobalBinNumber(reco_var,ih);

                    //std::cout<<ih<<" "<<reco_var<<" "<<reco_bin<<" JJ"<<std::endl;
                    //Find if this event should be oscillated
                    if(branch_variables[j][t]->GetOscillate()){
                        //Working
                        double true_var = *(static_cast<double*>(branch_variables[j][t]->GetTrueValue()));
                        double true_L = *(static_cast<double*>(branch_variables[j][t]->GetTrueL()));

                        double osc_Probability_sin = nu_model.oscProbSin(true_var, true_L);
                        double osc_Probability_sinsq = nu_model.oscProbSinSq(true_var, true_L);

                        std::cout<<ih<<" YARP "<<spec_osc_sin.hist[ih].GetName()<<" "<<true_var<<" "<<true_L<<" "<<osc_Probability_sin<<" "<<osc_Probability_sinsq<<" gWei "<<global_weight<<std::endl;

                        spec_osc_sinsq.hist[ih].Fill(reco_var, global_weight*osc_Probability_sinsq);
                        spec_osc_sin.hist[ih].Fill(reco_var, global_weight*osc_Probability_sin);
                        spec_central_value.hist[ih].Fill(reco_var,global_weight);
                        
                        std::cout<<"Reco: "<<reco_var<<" True: "<<true_var<<" L: "<<true_L<<" "<<osc_Probability_sin<<" "<<osc_Probability_sinsq<<" glob: "<<global_weight<<std::endl;
                    }else{
                        spec_central_value.hist[ih].Fill(reco_var,global_weight);
                        spec_osc_sinsq.hist[ih].Fill(reco_var, global_weight);
                        spec_osc_sin.hist[ih].Fill(reco_var, global_weight);
                        
                        //	std::cout<<reco_var<<" "<<std::endl;
                    }
                }
            }
        } //end of entry loop
    }//end of file loop


    /***************************************************************
     *		Now some clean-up and Writing
     * ************************************************************/

}

SBNgenerate::~SBNgenerate(){
    std::cout<<"~Closing SBNgenerate"<<std::endl;
    for(auto &f: files){
        f->Close();
    }
}


/***************************************************************
 *		Some virtual functions for selection and histogram filling
 * ************************************************************/

int SBNgenerate::WritePrecomputedOscSpecs(std::string tag){

    std::cout<<"SBNGenerate::WritePrecomputedOscSpecs()\t\t||\t\tWriting out "<<tag<<"SINXX_dm_"<<nu_model.mass_tag<<std::endl;
    spec_osc_sinsq.WriteOut(tag+"_SINSQ_dm_"+nu_model.mass_tag);
    spec_osc_sin.WriteOut(tag+"_SIN_dm_"+nu_model.mass_tag);

    return 0;
}

int SBNgenerate::WriteCVSpec(std::string tag){

    std::cout<<"SBNGenerate::WriteCVSpec()\t\t||\t\tWriting out "<<tag<<std::endl;
    spec_central_value.WriteOut(tag+"_CV");
    return 0;
}



bool SBNgenerate::EventSelection(int which_file){
    return true;
}

int SBNgenerate::FillHistograms(int file, int uni, double wei){
    return 0;
}
