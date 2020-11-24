#include "SBNcovariance.h"
#include <stdexcept>
#include <sstream>
#include <cassert>

using namespace sbn;

// for single photon: when we use root files with systematics already applied.
// usually 'useuniverse' is set to 'false' in this case.
SBNcovariance::SBNcovariance(std::string xmlname, bool useuniverse) : SBNconfig(xmlname, true, useuniverse) {
    otag = "SBN covariance::SBNcovariance\t||\t";

    std::cout <<otag<<"Start in DetSys mode. " << std::endl;

    universes_used = 0;
    tolerence_positivesemi = 1e-5;
    is_small_negative_eigenvalue = false;
    abnormally_large_weight = 1e3;//1e20;//20.0;

    bool restrict_variations = false;

    //Initialise the central value SBNspec.
    spec_central_value = SBNspec(xmlname,-1,false, useuniverse);

    int num_files = montecarlo_file.size();


    //unique variations used in the files
    variations.clear();
    std::vector<std::string> variations_tmp(systematic_name.begin(), systematic_name.end());
    std::sort(variations_tmp.begin(),variations_tmp.end());
    auto unique_iter = std::unique(variations_tmp.begin(), variations_tmp.end());
    variations.insert(variations.begin(),variations_tmp.begin(),unique_iter);

    std::cout<<otag<<" Construct for num_files=" << num_files << std::endl;

    std::vector<int> nentries(num_files,0);
    std::vector<int> used_montecarlos(num_files,0);

    files.resize(num_files,nullptr);
    trees.resize(num_files,nullptr);
    f_weights.resize(num_files,nullptr);

    montecarlo_additional_weight.resize(num_files,1.0);
    montecarlo_additional_weight_formulas.resize(num_files);

    int good_event = 0;

    for(int fid=0; fid < num_files; ++fid) {
        const auto& fn = montecarlo_file.at(fid);


        files[fid] = TFile::Open(fn.c_str());
        trees[fid] = (TTree*)(files[fid]->Get(montecarlo_name.at(fid).c_str()));
        nentries[fid]= (int)trees.at(fid)->GetEntries();

        //Some POT counting
        double pot_scale = 1.0;
        if(montecarlo_pot[fid]!=-1){
            pot_scale = this->plot_pot/montecarlo_pot[fid];
        }

        montecarlo_scale[fid] = montecarlo_scale[fid]*pot_scale;


        std::cout << otag<<"" << std::endl;
        std::cout << otag<<" TFile::Open() file=" << files[fid]->GetName() << " @" << files[fid] << std::endl;
        std::cout << otag<<" Has POT " <<montecarlo_pot[fid] <<" and "<<nentries[fid] <<" entries "<<std::endl;

        auto montecarlo_file_friend_treename_iter = montecarlo_file_friend_treename_map.find(fn);
        if (montecarlo_file_friend_treename_iter != montecarlo_file_friend_treename_map.end()) {
            std::cout<<otag<<" Detected friend trees "<<std::endl;

            auto montecarlo_file_friend_iter = montecarlo_file_friend_map.find(fn);
            if (montecarlo_file_friend_iter == montecarlo_file_friend_map.end()) {
                std::stringstream ss;
                ss << "Looked for filename=" << fn << " in fnmontecarlo_file_friend_iter, but could not be found... bad config?" << std::endl;
                throw std::runtime_error(ss.str());
            }

            for(int k=0; k < (*montecarlo_file_friend_iter).second.size(); k++){

                std::string treefriendname = (*montecarlo_file_friend_treename_iter).second.at(k);
                std::string treefriendfile = (*montecarlo_file_friend_iter).second.at(k);

                std::cout << otag<<" Adding a friend tree:  " <<treefriendname<<" from file: "<< treefriendfile <<std::endl;

                trees[fid]->AddFriend(treefriendname.c_str(),treefriendfile.c_str());
            }
        }


        for(const auto branch_variable : branch_variables[fid]) {
            //quick check that this branch associated subchannel is in the known chanels;
            int is_valid_subchannel = 0;
            for(const auto &name: fullnames){
                if(branch_variable->associated_hist==name){
                    std::cout<<otag<<" Found a valid subchannel for this branch: " <<name<<std::endl;
                    is_valid_subchannel++;
                }
            }
            if(is_valid_subchannel==0){
                std::cout<<otag<<" ERROR ERROR: This branch did not match one defined in the .xml : " <<branch_variable->associated_hist<<std::endl;
                std::cout<<otag<<" ERROR ERROR: There is probably a typo somehwhere in xml! "<<std::endl;
                exit(EXIT_FAILURE);

            }else if(is_valid_subchannel>1){
                std::cout<<otag<<" ERROR ERROR: This branch matched more than 1 subchannel!: " <<branch_variable->associated_hist<<std::endl;
                exit(EXIT_FAILURE);
            }
            //set the address of branch , old way required a branch. new is better
            //trees.at(fid)->SetBranchAddress(branch_variable->name.c_str(),    branch_variable->GetValue());
            branch_variable->branch_formula =  new TTreeFormula(("branch_form"+std::to_string(fid)).c_str(), branch_variable->name.c_str(), trees[fid]);
        }


        // set the address of 'reco_weight'.
        if(montecarlo_additional_weight_bool[fid]){
            //we have an additional weight we want to apply at run time, otherwise its just set at 1.
            std::cout<<"Setting Additional weight of : "<< montecarlo_additional_weight_names[fid].c_str()<<std::endl; 
            //trees[fid]->SetBranchAddress(montecarlo_additional_weight_names[fid].c_str(), &montecarlo_additional_weight[fid]); 
            montecarlo_additional_weight_formulas[fid] =  new TTreeFormula(("a_w"+std::to_string(fid)).c_str(), montecarlo_additional_weight_names[fid].c_str(), trees[fid]);
        }

        std::cout<<"Total Entries: "<<trees.at(fid)->GetEntries()<<" good event "<<good_event<<std::endl;
        trees.at(fid)->GetEntry(good_event);

    } // end fid


    // make a map and start filling, before filling find if already in map, if it is check size.
    std::cout << otag<<" Found " << variations.size() << " unique variations: " << std::endl;

    map_universe_to_var.clear();
    vec_universe_to_var.clear();
    num_universes_per_variation.clear();

    for(size_t vid=0; vid<variations.size(); ++vid) {
        const auto &v =  variations[vid];
        int in_n_files = 0;


        int n_cv_files = 0;

        for(int fid=0; fid < num_files; ++fid) {

            for(const auto branch_variable : branch_variables[fid]) {
                if(branch_variable->associated_systematic == v && branch_variable->central_value){
                    n_cv_files++;
                }
            }
        }

        int variation_length = std::count(systematic_name.begin(), systematic_name.end(), v)-n_cv_files;

        std::cout<<otag<<" "<<v<<" with length "<<variation_length<<" and an additional "<<n_cv_files <<" which are CV"<<std::endl;

        if(variation_length == -1){
            std::cout << otag << "Can't find this variation in the systematic variation vector!" << std::endl;
            exit(EXIT_FAILURE);
        }else if(variation_length == 0){
            std::cout << otag << "Can only find 1 central value file, can't build covarice matrix!" << std::endl;
        }


        //Variation length is 1. I.e we have a CV and 1 detsys shift
        variation_length = 1;

        for(int p=0; p < variation_length; p++){
            map_universe_to_var[num_universes_per_variation.size()] = v;
            vec_universe_to_var.push_back(vid);
            num_universes_per_variation.push_back(variation_length);// universe stuff
        }

        map_var_to_num_universe[v] = variation_length;
    }  //end looping over variations 


    m_variation_modes.resize(variations.size(),0);

    universes_used = num_universes_per_variation.size();

    std::cout << otag<<" -------------------------------------------------------------" << std::endl;
    std::cout << otag<<" Initilizing " << universes_used << " universes." << std::endl;
    std::cout << otag<<" -------------------------------------------------------------" << std::endl;

    std::vector<double> base_vec (spec_central_value.num_bins_total,0.0);

    std::cout << otag<<" Full concatanated vector has : " << spec_central_value.num_bins_total << std::endl;

    multi_vecspec.clear();
    multi_vecspec.resize(universes_used,base_vec);

    std::cout << otag<<" multi_vecspec now initilized of size :" << multi_vecspec.size() << std::endl;
    std::cout << otag<<" Reading the data files" << std::endl;
    watch.Reset();
    watch.Start();


    double reco_weight;
    for( size_t vid =0; vid < variations.size() ; vid++){
        const auto &v =  variations[vid];
        std::cout <<otag<<" On Variation " << v << ", total # of files: " << num_files << std::endl;
        for(int fid=0; fid< num_files; fid++){

            int nevents = std::min(montecarlo_maxevents[fid], nentries[fid]);
            for(int t=0; t< branch_variables[fid].size(); t++){
                const auto branch_var_jt = branch_variables[fid][t];

                //std::cout<<fid<<" "<<t<<" "<<branch_var_jt->associated_systematic<<" "<<branch_var_jt->associated_hist<<" on "<<v<<std::endl;
                if(branch_var_jt->associated_systematic == v){  // if we find the branch with right systematic variation
                    int ih = spec_central_value.map_hist.at(branch_var_jt->associated_hist);
                    double reco_var;
                    int reco_bin;

                    if(branch_var_jt->central_value == true ){
                        for(int i=0; i< nevents; i++){
                            trees[fid]->GetEntry(i);
                //            reco_var = *(static_cast<double*>(branch_var_jt->GetValue()));
                            branch_var_jt->GetFormula()->GetNdata();
                            double reco_var = branch_var_jt->GetFormula()->EvalInstance();


                            // calculate the reconstructed weight
                            reco_weight = 1.0;
                            if(montecarlo_additional_weight_bool[fid]){
                                montecarlo_additional_weight_formulas[fid]->GetNdata();
                                reco_weight = montecarlo_additional_weight_formulas[fid]->EvalInstance();
                                if(reco_weight!= reco_weight || reco_weight <0){
                                    std::cout<<"ERROR! the additional wight is nan or negative "<<reco_weight<<" Breakign!"<<std::endl;
                                    exit(EXIT_FAILURE);
                                }
                            }
                            reco_weight *= montecarlo_scale[fid]; 				
                           // std::cout<<"CV "<<reco_var<<" "<<reco_weight<<" "<<ih<<" "<<montecarlo_additional_weight_formulas[fid]->EvalInstance()<<" "<<montecarlo_scale[fid]<<std::endl;
                            spec_central_value.hist[ih].Fill(reco_var, reco_weight);
                        }
                    }else{
                        for(int i=0; i< nevents ; i++){
                            trees[fid]->GetEntry(i);
                            //reco_var = *(static_cast<double*>(branch_var_jt->GetValue()));
                            branch_var_jt->GetFormula()->GetNdata();
                            double reco_var = branch_var_jt->GetFormula()->EvalInstance();


                            // use CV spec to get the global bin number, even for systematically varied histograms
                            reco_bin = spec_central_value.GetGlobalBinNumber(reco_var,ih);
                            // calculate the reconstructed weight
                            reco_weight = 1.0;          
 
                            if(montecarlo_additional_weight_bool[fid]){
                                montecarlo_additional_weight_formulas[fid]->GetNdata();
                                reco_weight = montecarlo_additional_weight_formulas[fid]->EvalInstance();
                                if(reco_weight!= reco_weight || reco_weight <0){
                                    std::cout<<"ERROR! the additional wight is nan or negative "<<reco_weight<<" Breakign!"<<std::endl;
                                    exit(EXIT_FAILURE);
                                }
                            }
                            reco_weight *= montecarlo_scale[fid]; 				
                            //std::cout<<"SYS "<<reco_var<<" "<<reco_weight<<" "<<ih<<" "<<montecarlo_additional_weight_formulas[fid]->EvalInstance()<<" "<<montecarlo_scale[fid]<<std::endl;
                            if(reco_bin<0)continue; // if its not to be plotted
                            multi_vecspec[vid][reco_bin] += reco_weight;
                        }
                    }		
                } else continue;
            }// end of branch loop
        } // end of file loop
    } //end of variation loop



    watch.Stop();
    std::cout << otag<<" done CpuTime=" << watch.CpuTime() << " RealTime=" << watch.RealTime() << std::endl;

    /***************************************************************
     *		Now some clean-up and Writing
     * ************************************************************/

    for(auto f: files){
        std::cout << otag<<" TFile::Close() file=" << f->GetName() << " @" << f << std::endl;
        //        f->Close();
    }
    std::cout << otag<<" End" << std::endl;

}  //end constructor.



SBNcovariance::SBNcovariance(std::string xmlname) : SBNconfig(xmlname) {
    otag = "SBNcovariance::SBNcovariance\t||\t";

    std::cout <<otag<<"Start in EventWeight Mode." << std::endl;

    universes_used = 0;
    tolerence_positivesemi = 1e-5;
    is_small_negative_eigenvalue = false;
    abnormally_large_weight = 1e3;//1e20;//20.0;
    //Is there a Global weight to be applied to ALL weights, CV and otherwise, inside the eventweight class? 
    bnbcorrection_str = "NAN";//"TunedCentralValue_Genie";//"bnbcorrection_FluxHist";


    bool restrict_variations = false;
    
    variations_to_use = {"expskin_FluxUnisim","horncurrent_FluxUnisim","kminus_PrimaryHadronNormalization","kplus_PrimaryHadronFeynmanScaling","kzero_PrimaryHadronSanfordWang","nucleoninexsec_FluxUnisim","nucleonqexsec_FluxUnisim","nucleontotxsec_FluxUnisim","piminus_PrimaryHadronSWCentralSplineVariation","pioninexsec_FluxUnisim","pionqexsec_FluxUnisim","piontotxsec_FluxUnisim","piplus_PrimaryHadronSWCentralSplineVariation","genie_ccresAxial_Genie","genie_ncresAxial_Genie","genie_qema_Genie","genie_NC_Genie","genie_NonResRvbarp1pi_Genie","genie_NonResRvbarp2pi_Genie","genie_NonResRvp1pi_Genie","genie_NonResRvp2pi_Genie","genie_NonResRvbarp1piAlt_Genie","genie_NonResRvbarp2piAlt_Genie","genie_NonResRvp1piAlt_Genie","genie_NonResRvp2piAlt_Genie"};
    
    
    variations_to_use = {"AGKYpT1pi_Genie"," AGKYxF1pi_Genie"," AhtBY_Genie"," AxFFCCQEshape_Genie"," BhtBY_Genie"," CV1uBY_Genie"," CV2uBY_Genie"," DecayAngMEC_Genie"," EtaNCEL_Genie"," FrAbs_N_Genie"," FrAbs_pi_Genie"," FrCEx_N_Genie"," FrCEx_pi_Genie"," FrInel_N_Genie"," FrInel_pi_Genie"," FrPiProd_N_Genie"," FrPiProd_pi_Genie"," FracDelta_CCMEC_Genie"," FracPN_CCMEC_Genie"," MFP_N_Genie"," MFP_pi_Genie"," MaCCQE_Genie"    ," MaCCRES_Genie"," MaCOHpi_Genie"," MaNCEL_Genie"," MvCCRES_Genie"," NonRESBGvbarnCC1pi_Genie"," NonRESBGvbarnCC2pi_Genie"," NonRESBGvbarnNC1pi_Genie"," NonRESBGvbarnNC2pi_Genie"," NonRESBGvbarpCC1pi_Genie"," NonRESBGvbarpCC2pi_Genie"," NonRESBGvbarpNC1pi_Genie"," NonRESBGvbarpNC2pi_Genie"," NonRESBGvnCC1pi_Genie"," NonRESBGvnCC2pi_Genie"," NonRESBGvnNC1pi_Genie"," NonRESBGvnNC2pi_Genie"," NonRESBGvpCC1pi_Genie"," NonRESBGvpCC2pi_Genie"," NonRESBGvpNC1pi_Genie"," NonRESBGvpNC2pi_Genie"," NormCCMEC_Genie"," NormNCMEC_Genie"," R0COHpi_Genie"," RDecBR1eta_Genie"," RDecBR1gamma_Genie"," RPA_CCQE_Genie"," Theta_Delta2Npi_Genie"," VecFFCCQEshape_Genie"," XSecShape_CCMEC_Genie"};
    for(auto &s: variations_to_use){
        m_variations_to_use[s]=true;
    }



    std::map<std::string, int> parameter_sims;

    //Initialise the central value SBNspec.
    spec_central_value = SBNspec(xmlname,-1,false);

    int num_files = montecarlo_file.size();

    variations.clear();
    std::vector<std::string> variations_tmp;

    std::cout<<otag<<" Construct for num_files=" << num_files << std::endl;

    std::vector<int> nentries(num_files,0);
    std::vector<int> used_montecarlos(num_files,0);

    files.resize(num_files,nullptr);
    trees.resize(num_files,nullptr);
    f_weights.resize(num_files,nullptr);

    montecarlo_additional_weight.resize(num_files,1.0);
    montecarlo_additional_weight_formulas.resize(num_files);

    int good_event = 0;

    for(int fid=0; fid < num_files; ++fid) {
        const auto& fn = montecarlo_file.at(fid);

        files[fid] = TFile::Open(fn.c_str());
        trees[fid] = (TTree*)(files[fid]->Get(montecarlo_name.at(fid).c_str()));
        nentries[fid]= (int)trees.at(fid)->GetEntries();

        //Some POT counting
        double pot_scale = 1.0;
        if(montecarlo_pot[fid]!=-1){
            pot_scale = this->plot_pot/montecarlo_pot[fid];
        }

        montecarlo_scale[fid] = montecarlo_scale[fid]*pot_scale;


        std::cout << otag<<"" << std::endl;
        std::cout << otag<<" TFile::Open() file=" << files[fid]->GetName() << " @" << files[fid] << std::endl;
        std::cout << otag<<" Has POT " <<montecarlo_pot[fid] <<" and "<<nentries[fid] <<" entries "<<std::endl;

        auto montecarlo_file_friend_treename_iter = montecarlo_file_friend_treename_map.find(fn);
        if (montecarlo_file_friend_treename_iter != montecarlo_file_friend_treename_map.end()) {
            std::cout<<otag<<" Detected friend trees "<<std::endl;

            auto montecarlo_file_friend_iter = montecarlo_file_friend_map.find(fn);
            if (montecarlo_file_friend_iter == montecarlo_file_friend_map.end()) {
                std::stringstream ss;
                ss << "Looked for filename=" << fn << " in fnmontecarlo_file_friend_iter, but could not be found... bad config?" << std::endl;
                throw std::runtime_error(ss.str());
            }

            for(int k=0; k < (*montecarlo_file_friend_iter).second.size(); k++){

                std::string treefriendname = (*montecarlo_file_friend_treename_iter).second.at(k);
                std::string treefriendfile = (*montecarlo_file_friend_iter).second.at(k);

                std::cout << otag<<" Adding a friend tree:  " <<treefriendname<<" from file: "<< treefriendfile <<std::endl;

                trees[fid]->AddFriend(treefriendname.c_str(),treefriendfile.c_str());
            }
        }

        std::cout<<otag<<" Read variations & universe size" << std::endl;

        //trees.at(fid)->SetBranchAddress("mcweight", &(f_weights[fid]));
        trees.at(fid)->SetBranchAddress(montecarlo_eventweight_branch_names[fid].c_str(), &(f_weights[fid]));

        for(const auto branch_variable : branch_variables[fid]) {
            //quick check that this branch associated subchannel is in the known chanels;
            int is_valid_subchannel = 0;
            for(const auto &name: fullnames){
                if(branch_variable->associated_hist==name){
                    std::cout<<otag<<" Found a valid subchannel for this branch: " <<name<<std::endl;
                    is_valid_subchannel++;
                }
            }
            if(is_valid_subchannel==0){
                std::cout<<otag<<" ERROR ERROR: This branch did not match one defined in the .xml : " <<branch_variable->associated_hist<<std::endl;
                std::cout<<otag<<" ERROR ERROR: There is probably a typo somehwhere in xml! "<<std::endl;
                exit(EXIT_FAILURE);

            }else if(is_valid_subchannel>1){
                std::cout<<otag<<" ERROR ERROR: This branch matched more than 1 subchannel!: " <<branch_variable->associated_hist<<std::endl;
                exit(EXIT_FAILURE);
            }
            
            
            //trees.at(fid)->SetBranchAddress(branch_variable->name.c_str(),    branch_variable->GetValue());
            branch_variable->branch_formula =  new TTreeFormula(("branch_form"+std::to_string(fid)).c_str(), branch_variable->name.c_str(), trees[fid]);

        }

        if(montecarlo_additional_weight_bool[fid]){
            //we have an additional weight we want to apply at run time, otherwise its just set at 1.
            std::cout<<"Setting Additional weight of : "<< montecarlo_additional_weight_names[fid].c_str()<<std::endl; 
            //trees[fid]->SetBranchAddress(montecarlo_additional_weight_names[fid].c_str(), &montecarlo_additional_weight[fid]); 
            montecarlo_additional_weight_formulas[fid] =  new TTreeFormula(("a_w"+std::to_string(fid)).c_str(),montecarlo_additional_weight_names[fid].c_str(),trees[fid]);
        }

        std::cout<<"Total Entries: "<<trees.at(fid)->GetEntries()<<" good event "<<good_event<<std::endl;
        trees.at(fid)->GetEntry(good_event);

        const auto f_weight = f_weights[fid];
        if (f_weight == nullptr) {
            std::stringstream ss;
            ss << "Could not read weight branch for file=" << fid << std::endl;
            throw std::runtime_error(ss.str());
        }

        //This bit will calculate how many "universes" the file has. if ALL default is the inputted xml value

        std::cout<<"starting"<<std::endl;
        for(const auto& it : *f_weight) {
            std::cout<<"On : "<<it.first<<std::endl;
            if(it.first == bnbcorrection_str) {
                std::cout<<"Found a variation consistent with "<<bnbcorrection_str<<" . This will be instead applied as a general weight"<<std::endl;
                continue;    
            }



            if(variation_whitelist.size()> 0 ){
                    if(variation_whitelist.count(it.first)==0){
                        std::cout<<"Skipping "<<it.first<<" as its not in the WhiteList!!"<<std::endl;
                        continue;
                    }
            }

            if(variation_blacklist.size()> 0 ){
                    if(variation_blacklist.count(it.first)>0){
                        std::cout<<"Skipping "<<it.first<<" as it is the BlackList!!"<<std::endl;
                        continue;
                    }
            }


            /*
               if(it.first == "genie_all_Genie") {
               std::cout<<otag<<"Skipping genie_all_Genie!"<<std::endl;
               continue;
               }
               if(it.first == "genie_NC_Genie" ||  it.first == "genie_FormZone_Genie" || it.first == "genie_IntraNukeNmfp_Genie"){
            //|| it.first =="genie_IntraNukePImfp_Genie" || it.first == "genie_FormZone_Genie" || it.first == "piplus_PrimaryHadronSWCentralSplineVariation" || it.first == "") {
            //std::cout<<otag<<"Skipping genie_NC_Genie!"<<std::endl;
            continue;
            }
            */

            std::cout <<otag
                << it.first << " has " << it.second.size() << " montecarlos in file " << fid << std::endl;

            used_montecarlos[fid] += it.second.size();

            variations_tmp.push_back(it.first);
    }
        } // end fid

        std::sort(variations_tmp.begin(),variations_tmp.end());
        auto unique_iter = std::unique(variations_tmp.begin(), variations_tmp.end());
        variations.insert(variations.begin(),variations_tmp.begin(),unique_iter);



        //Variation Weight Maps Area
        m_variation_modes.resize(variations.size(),0);
        std::vector<std::string> s_formulas = this->buildWeightMaps();
        m_variation_weight_formulas.resize(num_files, std::vector<TTreeFormula*>(s_formulas.size()));


        for(int fid=0; fid < num_files; ++fid) {
            files[fid]->cd();
            for(int vid = 0; vid < variations.size(); vid++){ 
                m_variation_weight_formulas[fid][vid] =  new TTreeFormula(("weightMapsFormulas_"+std::to_string(fid)+"_"+std::to_string(vid)).c_str(), s_formulas[vid].c_str(),trees[fid]);
            }
        }



        // make a map and start filling, before filling find if already in map, if it is check size.
        std::cout << otag<<" Found " << variations.size() << " unique variations: " << std::endl;

        map_universe_to_var.clear();
        vec_universe_to_var.clear();
        num_universes_per_variation.clear();

        for(size_t vid=0; vid<variations.size(); ++vid) {
            const auto &v =  variations[vid];
            int max_variation_length = 0;
            int in_n_files = 0;

            std::cout<<otag<<" "<<v<<std::endl;
            //Lets loop over all trees

            for(size_t fid=0; fid<num_files; fid++){
                trees[fid]->GetEntry(good_event);

                //is this variation in this tree?
                int is_found = (*(f_weights[fid])).count(v);

                if(is_found==0){
                    std::cout<<otag<<"  WARNING @  variation " <<v<<"  in File " << montecarlo_file.at(fid)<<". Variation does not exist in file! "<<std::endl;
                }else{


                    int thissize = (*(f_weights[fid])).at(v).size(); // map lookup
                    in_n_files++;       
                    max_variation_length = std::max(thissize,max_variation_length);

                }

            }

            std::cout<<otag<<" "<<v<<" is of max length: "<<max_variation_length<<" and in "<<in_n_files<<" of "<<num_files<<" total files"<<std::endl;

            for(int p=0; p < max_variation_length; p++){
                map_universe_to_var[num_universes_per_variation.size()] = v;
                vec_universe_to_var.push_back(vid);
                num_universes_per_variation.push_back(max_variation_length);
            }

            map_var_to_num_universe[v] = max_variation_length;
        }

        std::cout << otag<<" File: 0 | " << montecarlo_file.at(0) << " has " << used_montecarlos.at(0) << " montecarlos" << std::endl;
        for(int i=1; i<num_files; i++){
            std::cout << otag<<" File: " << i <<" |  "<<montecarlo_file[i]<< " has " << used_montecarlos.at(i) << " montecarlos" << std::endl;
            if(used_montecarlos.at(i)!= used_montecarlos.at(i-1)){
                std::cerr << otag<<" Warning, number of universes for are different between files" << std::endl;
                std::cerr << otag<<" The missing universes are Set to weights of 1. Make sure this is what you want!" << std::endl;
                for(int j=0; j<num_files; j++){
                    if(universes_used < used_montecarlos.at(j)) 
                        universes_used = used_montecarlos.at(j);
                    std::cerr <<otag<<"File " << j << " montecarlos: " << used_montecarlos.at(j) << std::endl;
                }
            }
        }

        ///Quick check on minmax
        for(int v=0; v< variations.size(); v++){
            if(m_variation_modes[v]==1 && map_var_to_num_universe[variations[v]]!=2){
                std::cerr <<otag<<"ERROR! variation "<<variations[v]<<" is tagged as minmax mode, but has "<<map_var_to_num_universe[variations[v]]<<" universes (can only be 2)."<<std::endl;
                std::cout <<otag<<"ERROR! variation "<<variations[v]<<" is tagged as minmax mode, but has "<<map_var_to_num_universe[variations[v]]<<" universes (can only be 2)."<<std::endl;
                exit(EXIT_FAILURE);
            }
        }

        //But in reality we want the max universes to be the sum of all max variaitons across all files, NOT the sum over all files max variations.
        universes_used = num_universes_per_variation.size();

        std::cout << otag<<" -------------------------------------------------------------" << std::endl;
        std::cout << otag<<" Initilizing " << universes_used << " universes." << std::endl;
        std::cout << otag<<" -------------------------------------------------------------" << std::endl;

        std::vector<double> base_vec (spec_central_value.num_bins_total,0.0);

        std::cout << otag<<" Full concatanated vector has : " << spec_central_value.num_bins_total << std::endl;

        multi_vecspec.clear();
        multi_vecspec.resize(universes_used,base_vec);

        std::cout << otag<<" multi_vecspec now initilized of size :" << multi_vecspec.size() << std::endl;
        std::cout << otag<<" Reading the data files" << std::endl;
        watch.Reset();
        watch.Start();

        for(int j=0; j < num_files; j++){
            int nevents = std::min(montecarlo_maxevents[j], nentries[j]);
            std::cout << otag<<" Starting @ data file=" << files[j]->GetName() <<" which has "<<nevents<<" Events. "<<std::endl;
            size_t nbytes = 0;
            for(int i=0; i < nevents; i++) {
                if(i%100==0)std::cout<<otag<<" -- uni :"<<i<<" / "<<nevents<<std::endl;
                nbytes+= trees[j]->GetEntry(i);
                ProcessEvent(*(f_weights[j]),j,i);
            } //end of entry loop
            std::cout << otag<<" nbytes read=" << nbytes << std::endl;

        } //end of file loop

        watch.Stop();
        std::cout << otag<<" done CpuTime=" << watch.CpuTime() << " RealTime=" << watch.RealTime() << std::endl;

        /***************************************************************
         *		Now some clean-up and Writing
         * ************************************************************/

        for(auto f: files){
            std::cout << otag<<" TFile::Close() file=" << f->GetName() << " @" << f << std::endl;
            f->Close();
        }
        std::cout << otag<<" End" << std::endl;
    }


    void SBNcovariance::ProcessEvent(
            const std::map<std::string, 
            std::vector<eweight_type> >& thisfWeight,
            size_t fileid,
            int entryid) {

        double global_weight = 1.0;
        if( montecarlo_additional_weight_bool[fileid]){
            montecarlo_additional_weight_formulas[fileid]->GetNdata();
            global_weight = montecarlo_additional_weight_formulas[fileid]->EvalInstance();
        };//this will be 1.0 unless specifi



        global_weight *= montecarlo_scale[fileid];

        double additional_CV_weight = 1.0;

        const auto bnbcorr_iter = thisfWeight.find(bnbcorrection_str);
        if (bnbcorr_iter != thisfWeight.end())    additional_CV_weight *= (*bnbcorr_iter).second.front();

        if(std::isinf(global_weight) or (global_weight != global_weight)){
            std::stringstream ss;
            ss << "SBNcovariance::ProcessEvent\t||\tERROR  error @ " << entryid
                << " in File " << montecarlo_file.at(fileid) 
                << " as its either inf/nan: " << global_weight << std::endl;
            throw std::runtime_error(ss.str());
        }

        if(!EventSelection(fileid)) return;

        // precompute the weight size
        std::vector<double> weights(universes_used,global_weight);

        //Loop over all variations
        std::map<std::string, std::vector<eweight_type> >::const_iterator var_iter;
        int woffset = 0;
        int vid = 0;
        for(const auto& var : variations){


            //check if variation is in this file, if it isn't: then just push back 1's of appropiate number to keep universes consistent
            //this is of length of whatever the maximum length that was found in ANY file
            auto expected_num_universe_sz = map_var_to_num_universe.at(var); 

            //Grab newwer variation specfic weights;
            m_variation_weight_formulas[fileid][vid]->GetNdata();
            double indiv_variation_weight = m_variation_weight_formulas[fileid][vid]->EvalInstance();
            if((indiv_variation_weight!= indiv_variation_weight || indiv_variation_weight <0) && !montecarlo_fake[fileid]){
                std::cout<<"ERROR! the additional wight is nan or negative "<<indiv_variation_weight<<" Breakign!"<<std::endl;
                exit(EXIT_FAILURE);
            }
            //std::cout<<var<<" "<<indiv_variation_weight<<" "<<fileid<<" "<<vid<<std::endl;

            //is  
            var_iter = thisfWeight.find(var);
            int quick_fix = 0;

            if (var_iter == thisfWeight.end()) {
                //This we need to drop this for new version (where we add 1's)
                //std::cout<<var<<" is not in this universe, adding "<<expected_num_universe_sz<<" to woffset "<<woffset<<std::endl;
                //woffset += expected_num_universe_sz;
                //continue;
            }else {
                //first one is what we expect, second is whats directly inside the map.

                if (expected_num_universe_sz != (*var_iter).second.size()) {
                    std::stringstream ss;
                    //std::cout<< "Number of universes is not the max in this file" << std::endl;
                    //throw std::runtime_error(ss.str());

                    if( (*var_iter).second.size() > expected_num_universe_sz && var_iter != thisfWeight.end()){
                        ss << ". REALLY BAD!!  iter.size() " <<(*var_iter).second.size()<<" expected "<<expected_num_universe_sz<<" on var "<<var<<std::endl;
                        throw std::runtime_error(ss.str());
                    }
                }
                //so if this file contains smaller number of variations
                quick_fix = (*var_iter).second.size();
                //std::cout<< "So setting quick fix to the true number of universes in this varaiion : "<<quick_fix<< std::endl;


            }

            if (woffset >= weights.size()) {
                std::stringstream ss;
                ss << "woffset=" << woffset << " weights sz=" << weights.size() << " !" << std::endl;
                throw std::runtime_error(ss.str());
            }

            for(size_t wid0 = 0, wid1 = woffset; wid1 < (woffset + expected_num_universe_sz); ++wid0, ++wid1) {
                double wei = 1.0;
                //            std::cout<<"wid0 "<<wid0<<"/ "<<expected_num_universe_sz<<"  wid1 "<<wid1<<" / "<<weights.size()<<" woffset "<<woffset<<" quick_fix "<<quick_fix<<std::endl;

                if(wid0<quick_fix){
                    wei = (*var_iter).second[wid0];
                }

                bool is_inf = std::isinf(wei);
                bool is_nan = (wei != wei);

                if(is_inf or is_nan){
                    std::stringstream ss;
                    ss << "SBNcovariance::ProcessEvent\t||\t ERROR! Killing :: event # " << entryid
                        << " in File " << montecarlo_file.at(fileid) << " weight: " << wei << " global bnb: " << global_weight << " in " << var << std::endl;
                    //                throw std::runtime_error(ss.str());
                    wei = 1.0;                
                }

                if(wei > abnormally_large_weight){
                    std::cout<<"SBNcovariance::ProcessEvent\t||\tATTENTION!! HUGE weight: "<<wei<<" at "<<var<<" event "<<entryid<<" file "<<fileid<<std::endl;
                    wei=1.0;
                }

		if(montecarlo_fake[fileid]){
                	weights[wid1] *= wei;
		}
                else{
			 weights[wid1] *= wei*indiv_variation_weight;
		}
            }

            woffset += expected_num_universe_sz;
            vid++;
        }//end of all variations

        if (woffset != weights.size()) {
            std::stringstream ss;
            ss << "woffset=" << woffset << " weights sz=" << weights.size() << " !" << std::endl;
            throw std::runtime_error(ss.str());
        }

        //So the size of weights must equal global universes ya?
        if(universes_used != num_universes_per_variation.size()){
            std::stringstream ss;
            ss <<otag<<" ERROR "<<std::endl;
            ss <<"weights.size() "<<weights.size()<<std::endl;
            ss <<"universes_used "<<universes_used<<std::endl;
            ss <<"multi_vecspec.size() "<<multi_vecspec.size()<<std::endl;
            ss <<"num_universes_per_variation.size() "<<num_universes_per_variation.size()<<std::endl;
            throw std::runtime_error(ss.str());
        }

        for(int t=0; t < branch_variables[fileid].size(); t++) {


            const auto branch_var_jt = branch_variables[fileid][t];
            int ih = spec_central_value.map_hist.at(branch_var_jt->associated_hist);

            branch_var_jt->GetFormula()->GetNdata();
            double reco_var = branch_var_jt->GetFormula()->EvalInstance();
           // double reco_varold = *(static_cast<double*>(branch_var_jt->GetValue()));

            int reco_bin = spec_central_value.GetGlobalBinNumber(reco_var,ih);
            spec_central_value.hist[ih].Fill(reco_var, global_weight*additional_CV_weight);
            //std::cout<<reco_var<<" "<<reco_bin<<" "<<ih<<std::endl;

            for(int m=0; m<weights.size(); m++){
                if(reco_bin<0) continue;
                multi_vecspec[m][reco_bin] += weights[m];
            }
        }

        return;
    }


    /***************************************************************
     *		Some virtual functions for selection and histogram filling
     * ************************************************************/

    bool SBNcovariance::EventSelection(int which_file){
        //from here have access to vars_i  and vars_d  to make a selection
        return true;
    }

    int SBNcovariance::FillHistograms(int file, int uni, double wei){
        //double en = vars_d.at(file)[0];
        //multi_sbnspec.at(uni).hist.at(file).Fill(en, wei);
        return 0;
    }

    /***************************************************************
     *		And then form a covariance matrix (well 3)
     * ************************************************************/


    int SBNcovariance::FormCovarianceMatrix(std::string tag){

	// ********************* beginning of histograms normalization **********************

	//std::map<bool, std::string> map_shape_only{{true, "NCDeltaRadOverlaySM"}};
	std::map<bool, std::string> map_shape_only{{true, "NCDeltaRadOverlayLEE"}};
	//std::map<bool, std::string> map_shape_only{{false, "NCPi0NotCoh"}};
	//std::map<bool, std::string> map_shape_only{{true, "NCPi0NotCoh"}};

	for(auto const& imap : map_shape_only){
		bool lshape_only = imap.first;
		if(lshape_only == false) continue;
		std::string lname_subchannel = imap.second;
		if(is_verbose) std::cout << "SBNcovariance::FormCovariancematrix\t||\tSubchannel " << lname_subchannel << " will be constructed as shape-only matrix ? " << lshape_only << std::endl;	



		//save the toal number of events of CV for specific subchannels, and their global bin indices.
	        spec_central_value.CalcFullVector();
		std::vector<double> CV_tot_count;
		std::map<int, std::vector<double>> map_index_global_bin;
		for(auto const& lh:spec_central_value.hist){
			std::string lname = lh.GetName();
			if(lname.find(lname_subchannel) != std::string::npos){
				//store the max/min global bin index for histogram
				std::vector<double> lglobal_bin{spec_central_value.GetGlobalBinNumber(1, lname), spec_central_value.GetGlobalBinNumber(lh.GetNbinsX(), lname)};
				map_index_global_bin.insert( std::pair<int, std::vector<double>>( (int)CV_tot_count.size(), lglobal_bin)  );

				CV_tot_count.push_back(lh.Integral());
				
			}
		}
			
		//start modify 'multi_vecspec'
		for(int l=0; l< universes_used; l++){
		    //now, multi_vecspec[l] is a spectra vector of 1 universe
			std::string var_l = map_universe_to_var.at(l);
			if(var_l.find("UBGenie") == std::string::npos) continue;

			// loop over each histogram that has certain names		    
			for(auto const& lmap:map_index_global_bin){
				std::vector<double> lglobal_bin = lmap.second;

				// total # of events of a subchannel of this universe
				double uni_temp_count = std::accumulate(multi_vecspec[l].begin()+ lglobal_bin[0], multi_vecspec[l].begin()+lglobal_bin[1]+1, 0.0);
				for(int k=lglobal_bin[0]; k <= lglobal_bin[1]; k++)
					multi_vecspec[l][k] *= CV_tot_count[lmap.first]/uni_temp_count;
			}

		}
	}
	// ********************** end of normalization histograms **************************** 

        std::cout<<"SBNcovariance::FormCovariancematrix\t||\tStart" << std::endl;
        full_covariance.ResizeTo(num_bins_total, num_bins_total);
        frac_covariance.ResizeTo(num_bins_total, num_bins_total);
        full_correlation.ResizeTo(num_bins_total, num_bins_total);

        full_covariance.Zero();
        frac_covariance.Zero();
        full_correlation.Zero();

        vec_full_covariance.clear();
        vec_frac_covariance.clear();
        vec_full_correlation.clear();

        vec_full_covariance.resize(variations.size(),full_covariance);
        //auto vec_full_covariance2 = vec_full_covariance;
        vec_frac_covariance.resize(variations.size(),frac_covariance);
        vec_full_correlation.resize(variations.size(),full_correlation);

        int num_variations = variations.size();
        int num_bins_total2 = num_bins_total*num_bins_total;
        for(size_t vid=0; vid < variations.size(); ++vid) {
            const auto& v = variations[vid];
            map_var_to_matrix[v] = vid;
        }

        spec_central_value.CalcFullVector();

        std::vector<double> CV = spec_central_value.full_vector;

        // prepare pointer memory, incase we go to GPU
        double* a_CV = CV.data();
        int* a_num_universes_per_variation = num_universes_per_variation.data();
        int* a_vec_universe_to_var = vec_universe_to_var.data();

        double** a_multi_vecspec = new double*[multi_vecspec.size()];
        for(size_t m=0; m<multi_vecspec.size(); ++m)
            a_multi_vecspec[m] = multi_vecspec[m].data();

        double** a_vec_full_covariance  = new double*[vec_full_covariance.size()];
        //double** a_vec_full_covariance2 = new double*[vec_full_covariance.size()];
        double** a_vec_frac_covariance  = new double*[vec_frac_covariance.size()];
        double** a_vec_full_correlation = new double*[vec_full_correlation.size()];
        for(size_t k=0; k<vec_full_covariance.size(); ++k) {
            a_vec_full_covariance[k]  = (double*) vec_full_covariance[k].GetMatrixArray();
            // a_vec_full_covariance2[k] = (double*) vec_full_covariance2[k].GetMatrixArray();
            a_vec_frac_covariance[k]  = (double*) vec_frac_covariance[k].GetMatrixArray();
            a_vec_full_correlation[k] = (double*) vec_full_correlation[k].GetMatrixArray();
        }

        double* a_full_covariance  = full_covariance.GetMatrixArray();
        double* a_frac_covariance  = frac_covariance.GetMatrixArray();
        double* a_full_correlation = full_correlation.GetMatrixArray();

        std::cout << "SBNcovariance::FormCovariancematrix\t||\tForm variation sz= (" << num_bins_total << "X" << num_bins_total << ") covariance matrix(s)" << std::endl;
        watch.Reset();
        watch.Start();
#pragma acc parallel loop						\
        copy(a_vec_full_covariance[:num_variations][:num_bins_total2])	\
        copyin(a_vec_universe_to_var[:universes_used],			\
                a_CV[:num_bins_total],						\
                a_multi_vecspec[:universes_used][:num_bins_total],		\
                a_num_universes_per_variation[:universes_used])
        for(int k=0; k<universes_used; k++) {
            int varid = a_vec_universe_to_var[k];
            double vec_bot = ((double)a_num_universes_per_variation[k]);
            //next bit probably breaks the acc
            int varmode = m_variation_modes[varid];

            if(varmode==0){ //run as normal. 
#pragma acc loop seq
                for(int i=0; i<num_bins_total; i++) {
#pragma acc loop vector
                    for(int j=0; j<num_bins_total; j++) {
                        double vec_value = (a_CV[i]-a_multi_vecspec[k][i])*(a_CV[j]-a_multi_vecspec[k][j]) / vec_bot;
#pragma acc atomic update
                        a_vec_full_covariance[varid][i*num_bins_total+j] += vec_value;
                    }
                }
            }else if(varmode==1){
                //Instead, assign the covariance to be identicall the difference between this and the next universe (they come in 2's)
                for(int i=0; i<num_bins_total; i++) {
                       for(int j=0; j<num_bins_total; j++){
        	            a_vec_full_covariance[varid][i*num_bins_total+j] = (a_multi_vecspec[k][i]-a_multi_vecspec[k+1][i])* (a_multi_vecspec[k][j]-a_multi_vecspec[k+1][j]);
 	               }
                //a_vec_full_covariance[varid][i*num_bins_total+i] = fabs(a_multi_vecspec[k][i]-a_multi_vecspec[k+1][i]);
                }
                //we will also need to jump th universe count ahead by 1, just to skip the variation on the other side too.
                k++;

            }
        }
        watch.Stop();
        std::cout << "SBNcovariance::FormCovariancematrix\t||\tdone CpuTime=" << watch.CpuTime() << " RealTime=" << watch.RealTime() << std::endl;


        // std::cout << "SBNcovariance::FormCovariancematrix\t|| Calculating on the CPU" << std::endl;
        // watch.Reset();
        // watch.Start();
        // for(int i=0; i<num_bins_total; i++) {
        //   for(int j=0; j<num_bins_total; j++) {
        //     for(int k=0; k<universes_used; k++) {
        // 	int varid = a_vec_universe_to_var[k];
        // 	double vec_value = (a_CV[i]-a_multi_vecspec[k][i])*(a_CV[j]-a_multi_vecspec[k][j]);
        // 	vec_value /= ((double)a_num_universes_per_variation[k]);
        // 	a_vec_full_covariance2[varid][i*num_bins_total+j] += vec_value;
        //     }
        //   }
        // }
        // watch.Stop();
        // std::cout << "SBNcovariance::FormCovariancematrix\t|| done CpuTime=" << watch.CpuTime() << " RealTime=" << watch.RealTime() << std::endl;

        // std::cout << "SBNcovariance::FormCovariancematrix\t|| Comparing CPU & GPU" << std::endl;
        // double EPS = 1e-9;
        // for(int i=0; i<num_bins_total; i++) {
        //   for(int j=0; j<num_bins_total; j++) {
        //     for(int k=0; k<variations.size(); k++) {
        // 	if (std::abs(vec_full_covariance[k](i,j) - vec_full_covariance2[k](i,j)) > EPS)
        // 	  std::cout << "@(" << i << "," << j << "," << k << ") gpu=" << vec_full_covariance[k](i,j) << " cpu=" << vec_full_covariance2[k](i,j) << std::endl;
        //     }
        //   }
        // }
        // std::cout << "SBNcovariance::FormCovariancematrix\t|| done" << std::endl;

        watch.Reset();
        watch.Start();
        std::cout << "SBNcovariance::FormCovariancematrix\t||\tSumming over variations for covariance matrix" << std::endl;
        for(int vid=0; vid<variations.size(); ++vid) {
            full_covariance += vec_full_covariance[vid];
        }
        watch.Stop();
        std::cout << "SBNcovariance::FormCovariancematrix\t||\tdone CpuTime=" << watch.CpuTime() << " RealTime=" << watch.RealTime() << std::endl;

        std::cout<<"SBNcovariance::FormCovariancematrix\t||\tNow calculating fractional covariance and correlation matrix from full covariance."<<std::endl;
        watch.Reset();
        watch.Start();

#pragma acc parallel loop gang collapse(2)				\
        copyin(a_CV[:num_bins_total])						\
        copy(a_full_covariance[:num_bins_total2],				\
                a_frac_covariance[:num_bins_total2],				\
                a_full_correlation[:num_bins_total2],				\
                a_vec_full_covariance[:num_variations][:num_bins_total2],	\
                a_vec_frac_covariance[:num_variations][:num_bins_total2],	\
                a_vec_full_correlation[:num_variations][:num_bins_total2])
        for(int i=0; i < num_bins_total; i++) {
            for(int j=0; j < num_bins_total; j++) {
                a_frac_covariance[i*num_bins_total+j]  = a_full_covariance[i*num_bins_total+j]/(a_CV[i]*a_CV[j]);
//                std::cout<<"HALP "<<i<<" "<<j<<" "<<a_full_covariance[i*num_bins_total+j]<<" "<<a_CV[i]<<" "<<a_CV[j]<<" || "<<a_frac_covariance[i*num_bins_total+j]<<std::endl;
                a_full_correlation[i*num_bins_total+j] = a_full_covariance[i*num_bins_total+j]/(sqrt(a_full_covariance[i*num_bins_total+i])*sqrt(a_full_covariance[j*num_bins_total+j]));
#pragma acc loop
                for(int m=0; m<num_variations; m++){
                    a_vec_frac_covariance[m][i*num_bins_total+j]  = a_vec_full_covariance[m][i*num_bins_total+j]/(a_CV[i]*a_CV[j]) ;
                    a_vec_full_correlation[m][i*num_bins_total+j] = a_vec_full_covariance[m][i*num_bins_total+j]/(sqrt(a_vec_full_covariance[m][i*num_bins_total+i])*sqrt(a_vec_full_covariance[m][j*num_bins_total+j]));
                }
            }
        }
        watch.Stop();
        std::cout << "SBNcovariance::FormCovariancematrix\t||\tdone CpuTime=" << watch.CpuTime() << " RealTime=" << watch.RealTime() << std::endl;

        /************************************************************
         *			Saving to file				    *
         * *********************************************************/
        TFile *fout=new TFile((tag+".SBNcovar.root").c_str(),"RECREATE");
        fout->cd();
        full_covariance.Write("full_covariance",TObject::kWriteDelete);
        frac_covariance.Write("frac_covariance",TObject::kWriteDelete);
        full_correlation.Write("full_correlation",TObject::kWriteDelete);


        SBNchi collapse_chi(xmlname);

        TMatrixT<double > coll_correlation(num_bins_total_compressed,num_bins_total_compressed);
        TMatrixT<double > coll_frac_covariance(num_bins_total_compressed,num_bins_total_compressed);
        TMatrixT<double > coll_covariance(num_bins_total_compressed,num_bins_total_compressed);

        collapse_chi.CollapseModes(full_covariance, coll_covariance);
        spec_central_value.CollapseVector();

        for(int i=0; i<num_bins_total_compressed; i++){
            for(int j=0; j<num_bins_total_compressed; j++){
                coll_frac_covariance(i,j) = coll_covariance(i,j)/(spec_central_value.collapsed_vector.at(i)*spec_central_value.collapsed_vector.at(j)) ;
                coll_correlation(i,j)= coll_covariance(i,j)/(sqrt(coll_covariance(i,i))*sqrt(coll_covariance(j,j)));
            }
        }
        fout->cd();
        coll_covariance.Write("collapsed_covariance",TObject::kWriteDelete);
        coll_frac_covariance.Write("collapsed_frac_covariance",TObject::kWriteDelete);
        coll_correlation.Write("collapsed_correlation",TObject::kWriteDelete);


        TDirectory *individualDir = fout->GetDirectory("individualDir"); 
        if (!individualDir) { 
            individualDir = fout->mkdir("individualDir");       
        }
        fout->cd(); 
        individualDir->cd();

        for(int m=0; m< variations.size();m++){
            vec_full_correlation.at(m).Write( (variations.at(m)+"_full_correlation").c_str(), TObject::kWriteDelete);
            vec_frac_covariance.at(m).Write( (variations.at(m)+"_frac_covariance").c_str(), TObject::kWriteDelete);
            vec_full_covariance.at(m).Write( (variations.at(m)+"_full_covariance").c_str(), TObject::kWriteDelete);
        }




        std::vector<TH2D> h2_corr;
        std::vector<TH2D> h2_cov;
        std::vector<TH2D> h2_fcov;

        /*
           for(int m=0; m< variations.size();m++){
        //		vec_frac_covariance.at(m).Write((variations.at(m)+"_frac_covariance_"+tag).c_str() ,TObject::kWriteDelete);
        //		vec_full_covariance.at(m).Write((variations.at(m)+"_full_covariance_"+tag).c_str() ,TObject::kWriteDelete);
        //		vec_full_correlation.at(m).Write((variations.at(m)+"_full_correlation_"+tag).c_str() ,TObject::kWriteDelete);

        h2_corr.push_back(TH2D(vec_full_correlation.at(m)));
        h2_cov.push_back(TH2D(vec_full_covariance.at(m)));
        h2_fcov.push_back(TH2D(vec_frac_covariance.at(m)));

        h2_fcov.back().SetName((variations.at(m)+"_frac_covariance_"+tag).c_str());
        h2_corr.back().SetName((variations.at(m)+"_full_correlation_"+tag).c_str());
        h2_cov.back().SetName((variations.at(m)+"_full_covariance_"+tag).c_str());

        h2_fcov.back().Write();
        h2_cov.back().Write();
        h2_corr.back().Write();

        }
        */

        fout->Close();

        spec_central_value.WriteOut(tag);

        qualityTesting();

        std::cout<<"SBNcovariance::FormCovariancematrix\t||\tEnd" << std::endl;
        return 0;
    }


    int SBNcovariance::qualityTesting() {

        /************************************************************
         *		Quality Testing Suite			    *
         * *********************************************************/
        std::cout<<"SBNcovariance::qualityTesting\t||-----------------------------------------------------" << std::endl;
        std::cout<<"SBNcovariance::qualityTesting\t||----------------Quality Testing Suite"<<std::endl;
        std::cout<<"SBNcovariance::qualityTesting\t||-----------------------------------------------------" << std::endl;


        std::cout<<"SBNcovariance::qualityTesting\t||\tChecking if generated matrix is indeed a valid covariance matrix." << std::endl;
        std::cout<<"SBNcovariance::qualityTesting\t||\tFirst checking if matrix is symmetric." << std::endl;

        double max_sym_violation = 0;
        for(int i=0; i<num_bins_total; i++){
            for(int j=0; j<num_bins_total; j++){
                double tnp = fabs((full_covariance(j,i)-full_covariance(i,j))/(full_covariance(j,i)+full_covariance(i,j)));
                if(tnp > max_sym_violation) max_sym_violation = tnp;
            }
        }


        if(max_sym_violation < 1e-13){
            std::cout<<"SBNcovariance::qualityTesting\t||\tPASS: Generated covariance matrix is symmetric"<<std::endl;
        }else{
            std::cout<<"SBNcovariance::qualityTesting\t||\tERROR result is not symmetric! "<<max_sym_violation<<std::endl;
            for(int i=0; i<num_bins_total; i++){
                for(int j=0; j<num_bins_total; j++){
                    double tnp = fabs((full_covariance(j,i)-full_covariance(i,j))/(full_covariance(j,i)+full_covariance(i,j)));
                    if(full_covariance(i,j) != full_covariance(j,i)) std::cout<<i<<" "<<j<<" "<<full_covariance(i,j)<<" "<<full_covariance(j,i)<<" "<<tnp<<std::endl;
                }
            }
            std::cout<<"SBNcovariance::qualityTesting\t||\tERROR result is not symmetric!"<<std::endl;
            exit(EXIT_FAILURE);
        }

        std::cout<<"SBNcovariance::qualityTesting\t||\tChecking if generated matrix is positive semi-definite by looking at eigenvalues." << std::endl;
        //if a matrix is (a) real and (b) symmetric (checked above) then to prove positive semi-definite, we just need to check eigenvalues and >=0;
        TMatrixDEigen eigen (full_covariance);
        TVectorD eigen_values = eigen.GetEigenValuesRe();


        for(int i=0; i< eigen_values.GetNoElements(); i++){
            if(eigen_values(i)<0){
                is_small_negative_eigenvalue = true;
                if(fabs(eigen_values(i))> tolerence_positivesemi ){
                    std::cout << "SBNcovariance::qualityTesting\t||\tERROR contains (at least one)  negative eigenvalue: " << eigen_values(i) << std::endl;
                    exit(EXIT_FAILURE);
                }
            }
        }


        if(is_small_negative_eigenvalue){
            std::cout<<"SBNcovariance::qualityTesting\t||\tPASS: Generated covariance matrix is (allmost) positive semi-definite."<<std::endl;
            std::cout<<"SBNcovariance::qualityTesting\t||\tIt did contain small negative values of absolute value <= :"<<tolerence_positivesemi<<std::endl;
        }else{
            std::cout<<"SBNcovariance::qualityTesting\t||\tPASS: Generated covariance matrix is positive semi-definite."<<std::endl;
        }
        std::cout<<"SBNcovariance::qualityTesting\t||\tCongratulations, matrix is indeed a valid covariance matrix." << std::endl;

        return 0;
    }

    int SBNcovariance::PrintVariations(std::string tag){
        TFile *fout = new TFile(("SBNfit_variation_plots_"+tag+".root").c_str(),"recreate");
        fout->cd();

        if (access("variations",F_OK) == -1){
            mkdir("variations",0777);//Create a folder for pdf.
        }

        std::cout << "SBNcovariance::PrintVariations\t||\tStarting to Print all variations, this can take a little. " << std::endl;

        std::vector<TDirectory*> vec_dir;

        std::vector<std::vector<TCanvas*>> vec_canvas;

        for(const auto &v: variations){
            //std::cout<<"SBNcovariance::PrintVariations\t|| Preparing directory and canvases for variation: "<<v<<std::endl;
            fout->cd();
            vec_dir.push_back( fout->GetDirectory(v.c_str()));
            if (!vec_dir.back()) { 
                vec_dir.back() = fout->mkdir(v.c_str());       
            }
            vec_dir.back()->cd();

            std::vector<TCanvas *> tmpc;

            for(int i=0; i< spec_central_value.hist.size(); i++){
                tmpc.push_back(new TCanvas((fullnames.at(i)+"||"+v).c_str()));
                tmpc.back()->cd();
                TH1D * temp_cv_spec = (TH1D*)spec_central_value.hist.at(i).Clone((std::to_string(i)+v).c_str());
                temp_cv_spec->Scale(1,"width");

                tmpc.back()->cd();
                double maxval = temp_cv_spec->GetMaximum();
                if(maxval > 0) 	temp_cv_spec->SetMaximum(maxval*1.45);
                temp_cv_spec->SetStats(false);
                temp_cv_spec->SetLineColor(kBlack);
                temp_cv_spec->SetLineWidth(2);
                temp_cv_spec->GetXaxis()->SetTitle(fullnames.at(i).c_str());
                temp_cv_spec->GetYaxis()->SetTitle("Events/unit");
                temp_cv_spec->SetTitle((v + " || " +fullnames.at(i)).c_str());
                temp_cv_spec->DrawCopy("hist");

                delete temp_cv_spec;
            }
            vec_canvas.push_back(tmpc);	
        }

        time_t start_time = time(0);
        std::cout<<"SBNcovariance::PrintVariations\t||Starting universe loop [This can take a while!] as were plotting "<<universes_used<<" Histograms."<<std::endl;
        TRandom3 *rangen = new TRandom3(20);
        for(int m=0; m < universes_used; m++){

            //are we in a new variation? maybe. 

            if( m%((int)floor((double)universes_used/100.0))==0){
                std::cout<<"SBNcovariance::PrintVariations\t||\t On Universe: "<<m<<"/"<<universes_used<<" ---time elapsed since last notification: " << difftime(time(0), start_time) << " Seconds.\n";
                start_time = time(0);
            }

            std::string var =  map_universe_to_var.at(m);
            int which_matrix = map_var_to_matrix.at(var);
            vec_dir.at(which_matrix)->cd();
            SBNspec temp_spec(multi_vecspec.at(m), xmlname,false);

            for(int i=0; i< temp_spec.hist.size(); i++){
                vec_canvas.at(which_matrix).at(i)->cd();
                temp_spec.hist.at(i).Scale(1,"width");
                temp_spec.hist.at(i).SetLineColor((int)rangen->Uniform(300,1000));	
                temp_spec.hist.at(i).DrawCopy("same hist");

            }	
            //check to see if variation is over. if so
            if(m+1 != universes_used ) {
                if(var != map_universe_to_var[m+1]){
                    fout->cd();

                    for(int i=0; i< spec_central_value.hist.size(); i++){
                        fout->cd();
                        vec_dir.at(which_matrix)->cd();
                        vec_canvas.at(which_matrix).at(i)->cd();
                        TH1D * temp_cv_spec = (TH1D*)spec_central_value.hist.at(i).Clone((std::to_string(i)+var+"tmp2").c_str());
                        temp_cv_spec->Scale(1,"width");
                        temp_cv_spec->SetLineColor(kBlack);
                        temp_cv_spec->SetMarkerStyle(34);
                        temp_cv_spec->SetLineWidth(2);
                        temp_cv_spec->DrawCopy("same hist p");

                        fout->cd();
                        vec_dir.at(which_matrix)->cd();
                        vec_canvas.at(which_matrix).at(i)->Write();
                        vec_canvas.at(which_matrix).at(i)->SaveAs(("variations/Variation_"+tag+"_"+var+"_"+fullnames[i]+"_1D.pdf").c_str(),"pdf");
                        ;
                        delete temp_cv_spec;
                        delete vec_canvas.at(which_matrix).at(i);
                    }

                }
            }



        }//end universe loop

        std::cout << "SBNcovariance::PrintVariations\t||\tFinished. Just tidying up and writing TCanvas. " << std::endl;
        /*
           for(int v =0; v< variations.size(); v++){
           fout->cd();
           vec_dir.at(v)->cd();

           for(int i=0; i< spec_central_value.hist.size(); i++){
           vec_canvas.at(v).at(i)->cd();
           TH1D * temp_cv_spec = (TH1D*)spec_central_value.hist.at(i).Clone((std::to_string(i)+variations.at(v)+"tmp2").c_str());
           temp_cv_spec->Scale(1,"width");
           temp_cv_spec->SetLineColor(kBlack);
           temp_cv_spec->SetMarkerStyle(34);
           temp_cv_spec->SetLineWidth(2);
           temp_cv_spec->DrawCopy("same hist p");

           vec_canvas.at(v).at(i)->Write();
           delete temp_cv_spec;
           }
           }
           */

        fout->Close();
        return 0;
    }
    int SBNcovariance::PrintVariations_2D(std::string tag){
        TFile *fout = new TFile(("SBNfit_variation_plots_2D_"+tag+".root").c_str(),"recreate");
        fout->cd();

        std::cout << "SBNcovariance::PrintVariations_2D\t||\tStarting to Print all variations, this can take a little. " << std::endl;

        std::vector<TDirectory*> vec_dir; // vector of directories in root file to keep variatinos separate

        std::vector<std::vector<TCanvas*>> vec_canvas; //vector of vectors of canvases.  outer variations innner subchannels
        std::vector<std::vector<TH2D*>> vec_hist; //vector of vector of histograms
        //std::vector<std::vector<TCanvas*>> vec_canvas_2D; //vector of vectors of canvases.  outer variations innner subchannels

        variation_scale.open("variation_scale_"+tag+".txt",std::ofstream::trunc);
        sorted_variation_scale.open("sorted_variation_scale_"+tag+".txt",std::ofstream::trunc);
        correlation_scale.open("correlation_scale_"+tag+".txt",std::ofstream::trunc);
        correlation_constraint.open("correlation_constraint_"+tag+".txt",std::ofstream::trunc);

        std::vector<std::string> var_names;
        std::vector<double> var_percent;

        for(const auto &v: variations){  //for each variationss

            //std::cout<<"SBNcovariance::PrintVariations_2D\t|| Preparing directory and canvases for variation: "<<v<<std::endl;
            fout->cd(); //enter directory
            int which_matrix = map_var_to_matrix.at(v); 
            vec_dir.push_back( fout->GetDirectory(v.c_str())); 
            if (!vec_dir.back()) { 
                vec_dir.back() = fout->mkdir(v.c_str());    //make variation folder    
            }
            vec_dir.back()->cd();


            std::vector<TCanvas *> tmpc;
            std::vector<TH2D *> tmpTH2;
            std::vector<TH1D *> tmpTH1;

            for(int i=0; i< spec_central_value.hist.size(); i++){ 
                tmpc.push_back(new TCanvas((fullnames.at(i)+"||"+v).c_str()));
                tmpc.back()->cd();
                TH1D * temp_cv_spec = (TH1D*)spec_central_value.hist.at(i).Clone((std::to_string(i)+v).c_str()); //make clone so we can scale to arbitrary units.  Make clone so not messing up underlying values in case reused
                tmpTH1.push_back(temp_cv_spec);
                temp_cv_spec->Scale(1,"width");  //bunch of scaling and formatting

                double maxval = temp_cv_spec->GetMaximum()+temp_cv_spec->GetMaximum()/5.0;
                double minval= temp_cv_spec->GetMinimum()-temp_cv_spec->GetMinimum()/5.0;
                //int num_hist_bins=((maxval+0.2)-(minval-0.2))/0.2;
                int num_hist_bins=50;
                std::vector<double> bins;
                for(int j=0; j< temp_cv_spec->GetNbinsX(); j++){
                    bins.push_back(temp_cv_spec->GetBinLowEdge(j+1));
                }

                //then add the last edge,
                bins.push_back(  temp_cv_spec->GetBinLowEdge( temp_cv_spec->GetNbinsX())+ temp_cv_spec->GetBinWidth( temp_cv_spec->GetNbinsX()));
                double * tbins= &bins[0];


                //TH2D * temp_var_spec_2D = (TH2D*)spec_central_value.hist_2D.at(i).Clone((std::to_string(i)+v).c_str());

                TH2D* temp_var_spec_2D=new TH2D((fullnames.at(i)+v+"_hist").c_str(),(fullnames.at(i)+v+"_hist").c_str(), temp_cv_spec->GetNbinsX(), tbins, num_hist_bins, minval, maxval);
                tmpTH2.push_back(temp_var_spec_2D);


                if(maxval > 0) 	temp_cv_spec->SetMaximum(maxval*1.45); //set higher maximum for plot
                temp_cv_spec->SetStats(false);
                temp_cv_spec->SetLineColor(kBlack);
                temp_cv_spec->SetLineWidth(2.5);

                temp_cv_spec->GetXaxis()->SetLabelSize(0.03);
                temp_cv_spec->GetYaxis()->SetTitleOffset(1.1);
                temp_cv_spec->GetYaxis()->SetLabelSize(0.03);
                temp_cv_spec->GetYaxis()->SetTitleOffset(0.9);
                temp_cv_spec->GetXaxis()->SetTitle("Unit");
                temp_cv_spec->GetYaxis()->SetTitle("Arbitary Units");
                temp_cv_spec->GetYaxis()->SetTitleSize(0.045);
                temp_cv_spec->GetXaxis()->SetTitleSize(0.045);

                temp_cv_spec->SetTitle((v + " || " +fullnames.at(i)).c_str());
                std::vector<int> temp_indices = spec_central_value.map_tag_to_covariance_index.at(fullnames.at(i));
                int k=1;	
                TH1D * temp_cv_spec_2 = (TH1D*)temp_cv_spec->Clone((std::to_string(i)+v+"_2").c_str());
                //ofstream variation_scale;

                /*SBNchi collapse_chi(xmlname);
                  TMatrixT<double> coll_covariance(num_bins_total_compressed,num_bins_total_compressed);
                  TMatrixT<double> coll_correlation(num_bins_total_compressed,num_bins_total_compressed);
                  collapse_chi.CollapseModes(vec_full_covariance.at(which_matrix), coll_covariance);
                  for(int i=0; i<num_bins_total_compressed; i++){
                  for(int j=0; j<num_bins_total_compressed; j++){
                  coll_correlation(i,j)= coll_covariance(i,j)/(sqrt(coll_covariance(i,i))*sqrt(coll_covariance(j,j)));
                  }
                  }
                  correlation_constraint<<v<<" "<<coll_correlation(0,0)<< " " <<coll_correlation(0,1)<<" " <<coll_correlation(1,0)<<" " <<coll_correlation(1,1)<<std::endl;
                  */

                for(int j=temp_indices.at(0); j<=temp_indices.at(1); j++){
                    //std::cout<<"j="<<j<<std::endl;
                    temp_cv_spec_2->SetBinError(k,sqrt(full_covariance(j,j))/temp_cv_spec->GetBinWidth(k));
                    temp_cv_spec->SetBinError(k,sqrt(vec_full_covariance.at(which_matrix)(j,j))/temp_cv_spec->GetBinWidth(k));
                    //std::cout<<"total error band="<<sqrt(full_covariance(j,j))/temp_cv_spec->GetBinWidth(k)<<std::endl;
                    //std::cout<<"varation="<<fullnames.at(i)<<std::endl;
                    //std::cout<<"variation error band="<<sqrt(vec_full_covariance.at(which_matrix)(j,j))/temp_cv_spec->GetBinWidth(k)<<std::endl;

                    variation_scale<<which_matrix<<" "<<v<<" "<<sqrt(vec_frac_covariance.at(which_matrix)(j,j))*100<<std::endl;
                    var_names.push_back(v);
                    var_percent.push_back(sqrt(vec_frac_covariance.at(which_matrix)(j,j))*100);

                    //collapse full covariance then make correlation
                    //std::cout<<"matrix indices="<<j<<std::endl;
                    //std::cout<<"hist indices="<<k<<std::endl;

                    //int has_run=0;
                    //std::cout<<"temp_indices.at(0)="<<temp_indices.at(0)<<std::endl;
                    //std::cout<<"temp_indices.at(1)="<<temp_indices.at(1)<<std::endl;
                    //std::cout<<"temp_indices.size()="<<temp_indices.size()<<std::endl;
                    //attempt to get off diagonal entries very adhoc ask Mark for better method

                    correlation_scale<<v<<" ";
                    correlation_scale<<sqrt(vec_frac_covariance.at(which_matrix)(0,1))*100<<std::endl;

                    k++;
                }

                temp_cv_spec->SetLineWidth(1.0);
                temp_cv_spec_2->SetLineWidth(1.0);

                //temp_cv_spec->SetLineColorAlpha(kBlack, 0.35);
                //temp_cv_spec_2->SetLineColorAlpha(kRed, 0.35);
                temp_cv_spec->DrawCopy("E1 P");
                temp_var_spec_2D->Draw("colz same");
                temp_cv_spec_2->SetLineColor(kRed);
                temp_cv_spec_2->DrawCopy("same E1 P");
                temp_cv_spec->DrawCopy("same E1 P");
                TLegend *leg = new TLegend(0.61,0.75,0.89,0.89);
                leg->AddEntry(temp_cv_spec_2, "Full uncertainty","l");
                leg->AddEntry(temp_cv_spec, "Variation uncertainty","l");
                leg->Draw("same");

            }
            vec_canvas.push_back(tmpc);	 //vector of vector of canvases variations:signal with empty histogram plotted in each
            vec_hist.push_back(tmpTH2); //vector of vector of histograms (currently empty)
        }

        //end variation loop
        std::cout<<"SBNcovariance::PrintVariations_2D\t||Starting universe loop [This can take a while!] "<<std::endl;
        TRandom3 *rangen = new TRandom3(20);

        std::string last_var = "";
        int number_of_used_universes = 0;


        time_t start_time = time(0);

        for(int m=universes_used-1; m >=0; m--){  //looping over every universe
            std::string var = map_universe_to_var.at(m);  //in built map to universe, given universe which variation are we in
            int which_matrix = map_var_to_matrix.at(var); //integer position in vector canvas
            if(var != last_var){
                last_var = var;
                number_of_used_universes = 0;
            }

            if( m%((int)floor((double)universes_used/100.0))==0){
                std::cout<<"SBNcovariance::PrintVariations_2D\t||\t On Universe: "<<m<<"/"<<universes_used<<" ---time elapsed since last notification: " << difftime(time(0), start_time) << " Seconds.\n";
                start_time = time(0);
            }

            vec_dir.at(which_matrix)->cd();

            SBNspec temp_spec(multi_vecspec.at(m), xmlname,false);

            for(int i=0; i< temp_spec.hist.size(); i++){
                for(int j=0;j<temp_spec.hist.at(i).GetNbinsX(); j++){
                    TH1D * temp_var_spec = (TH1D*)temp_spec.hist.at(i).Clone("temp_var_spec_name");
                    temp_var_spec->Scale(1,"width");
                    double center     = temp_var_spec->GetBinCenter(j+1);
                    double value      = temp_var_spec->GetBinContent(j+1);
                    vec_hist.at(which_matrix).at(i)->Fill(center, value);
                    delete temp_var_spec;
                }


                vec_canvas.at(which_matrix).at(i)->Update();
            }
            number_of_used_universes++;
        }

        std::vector<size_t> sorted_indicies = SortIndexes(var_percent);
        for(int i= var_names.size()-1; i>-1; i--){
            size_t index = sorted_indicies[i];
            sorted_variation_scale<<var_names[index]<<" "<<var_percent[index]<<std::endl;
        }


        std::cout << "SBNcovariance::PrintVariations_2D\t||\tFinished. Just tidying up and writing TCanvas. " << std::endl;
        if (access("variations_2D",F_OK) == -1){
            mkdir("variations_2D",0777);//Create a folder for pdf.
        }

        for(int v =0; v< variations.size(); v++){
            fout->cd();
            vec_dir.at(v)->cd();

            for(int i=0; i< spec_central_value.hist.size(); i++){
                vec_canvas.at(v).at(i)->cd();
                // TH1D * temp_cv_spec = (TH1D*)spec_central_value.hist.at(i).Clone((std::to_string(i)+variations.at(v)+"tmp2").c_str());
                /*
                   TLegend *l = new TLegend(0.61,0.61,0.89,0.89);
                   l->SetLineColor(kWhite);
                   l->SetLineWidth(0);
                   temp_cv_spec->Scale(1,"width");
                   temp_cv_spec->SetLineColor(kBlack);
                   temp_cv_spec->SetMarkerStyle(48);
                   temp_cv_spec->SetLineWidth(3);
                   temp_cv_spec->DrawCopy("same hist");
                   l->AddEntry(temp_cv_spec,"Central Value","l");
                   l->Draw("same");
                   */
                vec_canvas.at(v).at(i)->Write();
                vec_canvas.at(v).at(i)->SaveAs(("variations_2D/Variation_"+tag+"_"+variations[v]+"_"+fullnames[i]+"_2D.pdf").c_str(),"pdf");

                //delete temp_cv_spec;
            }


        }

        fout->Close();
        return 0;
    }


    int SBNcovariance::PrintMatricies(std::string tag) {
        std::cout << "SBNcovariance::PrintMatricies\t||\tStart" << std::endl;

        TFile* fout = new TFile(("SBNfit_covariance_plots_"+tag+".root").c_str(),"recreate");
        fout->cd();

        gStyle->SetOptStat(0);

        this->plot_one(full_correlation, "SBNfit_correlation_matrix_"+tag, fout, true,false);
        this->plot_one(full_covariance,  "SBNfit_covariance_matrix_"+tag, fout, true,false);
        this->plot_one(frac_covariance,  "SBNfit_fractional_covariance_matrix_"+tag, fout, true,false);

        SBNchi collapse_chi(xmlname);

        TMatrixT<double > coll_correlation(num_bins_total_compressed,num_bins_total_compressed);
        TMatrixT<double > coll_frac_covariance(num_bins_total_compressed,num_bins_total_compressed);
        TMatrixT<double > coll_covariance(num_bins_total_compressed,num_bins_total_compressed);

        collapse_chi.CollapseModes(full_covariance, coll_covariance);

        for(int i=0; i<num_bins_total_compressed; i++){
            for(int j=0; j<num_bins_total_compressed; j++){
                coll_frac_covariance(i,j) = coll_covariance(i,j)/(spec_central_value.collapsed_vector.at(i)*spec_central_value.collapsed_vector.at(j)) ;
                coll_correlation(i,j)= coll_covariance(i,j)/(sqrt(coll_covariance(i,i))*sqrt(coll_covariance(j,j)));
            }
        }

        TH2D h2_coll_corr(coll_correlation);
        h2_coll_corr.SetName("coll_corr");
        TCanvas *c_coll_corr = new TCanvas("collapsed correlation matrix");
        c_coll_corr->cd();
        c_coll_corr->SetFixedAspectRatio();
        h2_coll_corr.Draw("colz");
        h2_coll_corr.SetTitle("Collapsed correlation matrix");
        h2_coll_corr.GetXaxis()->SetTitle("Reco Bin i");
        h2_coll_corr.GetYaxis()->SetTitle("Reco Bin j");

        c_coll_corr->SetRightMargin(0.150);

        int use_coll_corr =0;
        for(int im =0; im<num_modes; im++){
            for(int id =0; id<num_detectors; id++){
                for(int ic = 0; ic < num_channels; ic++){
                    TLine *lv = new TLine(0, num_bins.at(ic)+use_coll_corr, num_bins_total_compressed, num_bins.at(ic)+use_coll_corr);
                    TLine *lh = new TLine(num_bins.at(ic)+use_coll_corr,0, num_bins.at(ic)+use_coll_corr, num_bins_total_compressed);
                    lv->SetLineWidth(1.5);
                    lh->SetLineWidth(1.5);
                    use_coll_corr+=num_bins.at(ic);
                    lv->Draw();
                    lh->Draw();
                }
            }
        }
        c_coll_corr->Write();


        TH2D h2_coll_frac(coll_frac_covariance);
        h2_coll_frac.SetName("coll_frac");
        TCanvas *c_coll_frac = new TCanvas("collapsed fractional covariance matrix");
        c_coll_frac->cd();
        c_coll_frac->SetFixedAspectRatio();
        h2_coll_frac.Draw("colz");
        h2_coll_frac.SetTitle("Collapsed fractional covariance matrix");
        h2_coll_frac.GetXaxis()->SetTitle("Reco Bin i");
        h2_coll_frac.GetYaxis()->SetTitle("Reco Bin j");

        c_coll_frac->SetRightMargin(0.150);

        int use_coll_frac =0;
        for(int im =0; im<num_modes; im++){
            for(int id =0; id<num_detectors; id++){
                for(int ic = 0; ic < num_channels; ic++){
                    TLine *lv = new TLine(0, num_bins.at(ic)+use_coll_frac, num_bins_total_compressed, num_bins.at(ic)+use_coll_frac);
                    TLine *lh = new TLine(num_bins.at(ic)+use_coll_frac,0, num_bins.at(ic)+use_coll_frac, num_bins_total_compressed);
                    lv->SetLineWidth(1.5);
                    lh->SetLineWidth(1.5);
                    use_coll_frac+=num_bins.at(ic);
                    lv->Draw();
                    lh->Draw();
                }
            }
        }
        c_coll_frac->Write();


        TH2D h2_coll_full(coll_covariance);
        h2_coll_full.SetName("coll_full");
        TCanvas *c_coll_full = new TCanvas("collapsed covariance matrix");
        c_coll_full->cd();
        c_coll_full->SetFixedAspectRatio();
        h2_coll_full.Draw("colz");
        h2_coll_full.SetTitle("Collapsed covariance matrix");
        h2_coll_full.GetXaxis()->SetTitle("Reco Bin i");
        h2_coll_full.GetYaxis()->SetTitle("Reco Bin j");

        c_coll_full->SetRightMargin(0.150);

        int use_coll_full =0;
        for(int im =0; im<num_modes; im++){
            for(int id =0; id<num_detectors; id++){
                for(int ic = 0; ic < num_channels; ic++){
                    TLine *lv = new TLine(0, num_bins.at(ic)+use_coll_full, num_bins_total_compressed, num_bins.at(ic)+use_coll_full);
                    TLine *lh = new TLine(num_bins.at(ic)+use_coll_full,0, num_bins.at(ic)+use_coll_full, num_bins_total_compressed);
                    lv->SetLineWidth(1.5);
                    lh->SetLineWidth(1.5);
                    use_coll_full+=num_bins.at(ic);
                    lv->Draw();
                    lh->Draw();
                }
            }
        }
        c_coll_full->Write();

        if (access("matrix_plots",F_OK) == -1){
            mkdir("matrix_plots",0777);//Create a folder for pdf.
        }


        for(int m=0; m< variations.size();m++){
            this->plot_one(vec_full_correlation.at(m), "matrix_plots/varplot_Correlation_"+tag+"_"+variations.at(m), fout,true,true);
            this->plot_one(vec_frac_covariance.at(m), "matrix_plots/varplot_Fractional_Covariance_"+tag+"_"+variations.at(m), fout,true,true);
            this->plot_one(vec_full_covariance.at(m), "matrix_plots/varplot_Full_Covariance_"+tag+"_"+variations.at(m), fout,true,true);
        }

        fout->cd();
        fout->Close();

        std::cout << "SBNcovariance::PrintMatricies\t||\tEnd" << std::endl;
        return 0;
    }


    int SBNcovariance::plot_one(TMatrixD matrix, std::string tag, TFile *fin, bool plot_pdf, bool indiv){
        fin->cd();
        if(indiv){
            TDirectory *individualDir = fin->GetDirectory("individualDir"); 
            if (!individualDir) { 
                individualDir = fin->mkdir("individualDir");       
            }
            fin->cd(); 
            individualDir->cd();
        }
        TH2D h2_full(matrix);
        h2_full.SetName((tag+"_th2d").c_str());
        TCanvas *c_full = new TCanvas((tag+"_canvas").c_str());
        TPad *p_full = (TPad*)c_full->cd();
        c_full->SetFixedAspectRatio();
        h2_full.Draw("colz");
        h2_full.SetTitle(tag.c_str());
        h2_full.GetXaxis()->SetTitle("Global Bin Number");
        h2_full.GetYaxis()->SetTitle(" ");
        h2_full.GetYaxis()->SetLabelSize(0);

        c_full->SetFrameFillColor(kWhite);
        c_full->SetFillColor(kWhite);
        p_full->SetFillColor(kWhite);


        c_full->SetRightMargin(0.150);
        c_full->SetLeftMargin(0.250);
        c_full->SetTopMargin(0.10);
        int use_full =0;

        double percent_left = 0.15;
        double nice_shift = num_bins_total*0.02;

        for(int im =0; im<num_modes; im++){
            for(int id =0; id<num_detectors; id++){
                for(int ic = 0; ic < num_channels; ic++){
                    for(int isc = 0; isc < num_subchannels.at(ic); isc++){


                        std::string mode_det = mode_names[im] +" " +detector_names[id];
                        std::string chan_sub = channel_names[ic]+" "+subchannel_names[ic][isc];


                        TText * tmd = new TText(-num_bins_total*percent_left*0.15, use_full+nice_shift*0.5, (mode_det+" "+chan_sub).c_str() );

                        //TText * tmd = new TText(use_full*1.05, num_bins_total*1.015, chan_sub.c_str());
                        //TText * tcs = new TText(use_full*1.05, num_bins_total*1.055, mode_det.c_str());
                        tmd->SetTextColor(kBlack);
                        //tcs->SetTextColor(kBlack);
                        tmd->SetTextSize(0.03);
                        tmd->SetTextAlign(31);
                        //tcs->SetTextSize(0.03);
                        tmd->Draw();
                        //tcs->Draw();


                        /*
                           TText * tlow_bin = new TText(-num_bins_total*percent_left, use_full+nice_shift*0.5, to_string_prec(bin_edges[ic].front(),0).c_str());
                           TText * thigh_bin = new TText(-num_bins_total*percent_left, (use_full+num_bins[ic])-nice_shift*1.4, to_string_prec(bin_edges[ic].back(),0).c_str());
                           tlow_bin->SetTextSize(0.02);
                           thigh_bin->SetTextSize(0.02);
                           tlow_bin->Draw();
                           thigh_bin->Draw();

                           TText * tunit = new TText(-num_bins_total*percent_left, use_full+0.5*num_bins[ic], channel_units[ic].c_str());
                           tunit->SetTextSize(0.03);
                           tunit->Draw();
                           */

                        if(isc<num_subchannels[ic]-1){
                            TLine *lscv = new TLine(-num_bins_total*percent_left, num_bins.at(ic)+use_full, num_bins_total, num_bins.at(ic)+use_full);
                            TLine *lsch = new TLine(num_bins.at(ic)+use_full,0, num_bins.at(ic)+use_full, num_bins_total*1.045);
                            lscv->SetLineWidth(3);
                            lsch->SetLineWidth(3);
                            lscv->SetLineColor(kRed);
                            lsch->SetLineColor(kRed);
                            lscv->SetLineStyle(9);
                            lsch->SetLineStyle(9);

                            lscv->Draw();
                            lsch->Draw();

                            use_full+=num_bins.at(ic);

                        }
                    }
                    TLine *lv = new TLine(-num_bins_total*percent_left, num_bins.at(ic)+use_full, num_bins_total, num_bins.at(ic)+use_full);
                    TLine *lh = new TLine(num_bins.at(ic)+use_full,0, num_bins.at(ic)+use_full, num_bins_total*1.045);
                    lv->SetLineWidth(3);
                    lh->SetLineWidth(3);
                    lv->SetLineColor(kRed);
                    lh->SetLineColor(kRed);
                    use_full+=num_bins.at(ic);
                    lv->Draw();
                    lh->Draw();

                }
            }
        }


        c_full->Write();
        if(plot_pdf) c_full->SaveAs((tag+".pdf").c_str(),"pdf");


        return 0;
    }

    int SBNcovariance::DoConstraint(int which_signal, int which_constraint){
        std::cout<<"----------------Starting covariance Constraint --------------------"<<std::endl;

        SBNchi collapse_chi(xmlname);
        TMatrixT<double> coll_covariance(num_bins_total_compressed,num_bins_total_compressed);
        collapse_chi.CollapseModes(full_covariance, coll_covariance);

        int num_bins_sig=num_bins[which_signal];
        std::cout<<"signal bins="<<num_bins_sig<<std::endl;
        int num_bins_con=num_bins[which_constraint];
        std::cout<<"constraint bins="<<num_bins_con<<std::endl;
        int num_bins_tot=num_bins_sig+num_bins_con;

        spec_central_value.CollapseVector();
        std::vector<double> n_events = spec_central_value.collapsed_vector; 
        std::vector<double> n_sig;
        std::vector<double> n_con;

        for(int i=0; i< num_bins_sig; i++){
            n_sig.push_back(n_events[i]);
        }

        for(int i=num_bins_sig; i< n_events.size(); i++){
            n_con.push_back(n_events[i]);
        }

        std::cout<<"1g1p: "<<std::accumulate(n_sig.begin(),n_sig.end(),0.0)<<std::endl;
        std::cout<<"2g1p: "<<std::accumulate(n_con.begin(),n_con.end(),0.0)<<std::endl;

        for (int i=0; i<num_bins_tot;i++){
            std::cout<<"vector entry "<<i<<"-"<<n_events[i]<<std::endl;
        }
        TMatrixD m_covar_full(num_bins_tot,num_bins_tot);
        TMatrixD m_covar_full_sys = coll_covariance;
        TMatrixD m_inv_full(num_bins_tot,num_bins_tot);
        TMatrixD m_covar_con(num_bins_tot,num_bins_tot);
        TMatrixD m_covar_con_sys(num_bins_tot,num_bins_tot);

        for(int i=0; i<num_bins_tot;i++)
        {
            for(int j=0; j<num_bins_tot;j++)
            {
                if(i==j){
                    m_covar_full(i,j)=coll_covariance(i,j)+n_events[i];
                    m_inv_full(i,j)=coll_covariance(i,j)+n_events[i];

                }
                else{
                    m_covar_full(i,j)=coll_covariance(i,j);
                    m_inv_full(i,j)=coll_covariance(i,j);
                }
            }
        }
        std::cout<<"covar full pre-inversion"<<std::endl;
        m_covar_full.Print();

        m_inv_full.Invert();
        std::cout<<"covar invert"<<std::endl;
        m_inv_full.Print();

        TMatrixD m_inv_con(num_bins_tot,num_bins_tot);
        for(int i=0; i<num_bins_tot;i++){
            for(int j=0; j<num_bins_tot; j++){
                if(i==j  && j>=num_bins_sig){
                    m_inv_con(i,j)=m_inv_full(i,j)+1/n_events[i];
                    m_covar_con(i,j)=m_inv_full(i,j)+1/n_events[i];
                }
                else{
                    m_inv_con(i,j)=m_inv_full(i,j);
                    m_covar_con(i,j)=m_inv_full(i,j);
                }
            } 
        }
        m_covar_con.Invert();

        for(int i=0; i<num_bins_tot;i++)
        {
            for(int j=0; j<num_bins_tot;j++)
            {

                if(i==j){
                    m_covar_con_sys(i,j)=m_covar_con(i,j)-n_events[i];
                }
                else{
                    m_covar_con_sys(i,j)=m_covar_con(i,j);
                }
            }
        }
        std::cout<<"covar con sys"<<std::endl;
        m_covar_con_sys.Print();

        for(int i=0;i<num_bins_tot;i++){
            double uncon=sqrt(m_covar_full_sys(i,i))/n_events[i];
            double con=sqrt(m_covar_con_sys(i,i))/n_events[i];

            std::cout<<"bin="<<i<<"/events="<<n_events[i]<<"/uncon="<<uncon<<"/con="<<con<<" ratio: "<<uncon/con<<std::endl;
        }


        TFile *fnew = new TFile("constraint_test.root","recreate");
        TMatrixD cov_before = coll_covariance.GetSub(0,num_bins_sig-1,0,num_bins_sig-1);
        TMatrixD cov_after = m_covar_con_sys.GetSub(0,num_bins_sig-1,0,num_bins_sig-1);

        cov_before.Write("before");
        cov_after.Write("after");



        return 0;
    }

    std::vector<double> SBNcovariance::DoConstraint(int which_signal, int which_constraint, std::string tag){

        return this->DoConstraint(which_signal, which_constraint, tag, -1);
    }

    std::vector<double> SBNcovariance::DoConstraint(int which_signal, int which_constraint, std::string tag, int which_var){
        std::cout<<"----------------Starting covariance Constraint --------------------"<<std::endl;

        SBNchi collapse_chi(xmlname);
        TMatrixT<double> coll_covariance(num_bins_total_compressed,num_bins_total_compressed);

        if(which_var<0){
            collapse_chi.CollapseModes(full_covariance, coll_covariance);
        }else{
            collapse_chi.CollapseModes(vec_full_covariance.at(which_var), coll_covariance);
        }

        //Does this assume only 2 channels?
        int num_bins_sig = num_bins[which_signal];
        std::cout<<"signal bins="<<num_bins_sig<<std::endl;
        int num_bins_con=num_bins[which_constraint];
        std::cout<<"constraint bins="<<num_bins_con<<std::endl;
        int num_bins_tot=num_bins_sig+num_bins_con;


        spec_central_value.CollapseVector();
        std::vector<double> n_events = spec_central_value.collapsed_vector; 
        std::vector<double> n_sig;
        std::vector<double> n_con;


        //This assumes that the signal is first in order and the background is next AND THAT IT. 
        for(int i=0; i< num_bins_sig; i++){
            n_sig.push_back(n_events[i]);
        }
        for(int i=num_bins_sig; i< n_events.size(); i++){
            n_con.push_back(n_events[i]);
        }

        std::cout<<"Signal: "<<std::accumulate(n_sig.begin(),n_sig.end(),0.0)<<std::endl;
        std::cout<<"Constraint: "<<std::accumulate(n_con.begin(),n_con.end(),0.0)<<std::endl;

        for (int i=0; i<num_bins_tot;i++){
            std::cout<<"Vector entry "<<i<<"-"<<n_events[i]<<std::endl;
        }

        TMatrixD m_covar_full(num_bins_tot,num_bins_tot);
        TMatrixD m_covar_full_sys = coll_covariance;
        TMatrixD m_inv_full(num_bins_tot,num_bins_tot);
        TMatrixD m_covar_con(num_bins_tot,num_bins_tot);
        TMatrixD m_covar_con_sys(num_bins_tot,num_bins_tot);

        for(int i=0; i<num_bins_tot;i++)
        {
            for(int j=0; j<num_bins_tot;j++)
            {
                if(i==j){
                    m_covar_full(i,j)=coll_covariance(i,j)+n_events[i];
                    m_inv_full(i,j)=coll_covariance(i,j)+n_events[i];

                }
                else{
                    m_covar_full(i,j)=coll_covariance(i,j);
                    m_inv_full(i,j)=coll_covariance(i,j);
                }
            }
        }
        std::cout<<"covar full pre-inversion"<<std::endl;
        m_covar_full.Print();

        m_inv_full.Invert();
        std::cout<<"covar invert"<<std::endl;
        m_inv_full.Print();

        std::cout<<"check inversion"<<std::endl;
        if(check_inversion(m_covar_full,m_inv_full,num_bins_tot)){
            std::cout<<"Initial inversion good"<<std::endl;
        }
        else{
            std::cout<<"bad inversion"<<std::endl;
            //throw;
        }


        TMatrixD m_inv_con(num_bins_tot,num_bins_tot);
        for(int i=0; i<num_bins_tot;i++){
            for(int j=0; j<num_bins_tot; j++){
                if(i==j  && j>=num_bins_sig){
                    m_inv_con(i,j)=m_inv_full(i,j)+1/n_events[i];
                    m_covar_con(i,j)=m_inv_full(i,j)+1/n_events[i];
                }
                else{
                    m_inv_con(i,j)=m_inv_full(i,j);
                    m_covar_con(i,j)=m_inv_full(i,j);
                }
            } 
        }
        m_covar_con.Invert();

        if(check_inversion(m_covar_con,m_inv_con,num_bins_tot)){
            std::cout<<"Initial inversion good"<<std::endl;
        }
        else{
            std::cout<<"bad inversion"<<std::endl;
            //throw;
        }

        for(int i=0; i<num_bins_tot;i++)
        {
            for(int j=0; j<num_bins_tot;j++)
            {

                if(i==j){
                    m_covar_con_sys(i,j)=m_covar_con(i,j)-n_events[i];
                }
                else{
                    m_covar_con_sys(i,j)=m_covar_con(i,j);
                }
            }
        }
        std::cout<<"covar con sys"<<std::endl;
        m_covar_con_sys.Print();
        //constraint_table.open(("constraint_table"+tag+"_knob_"+variations[which_var]+".txt").c_str(),std::ofstream::trunc);
        constraint_table.open(("constraint_table_1bin_"+tag+".txt").c_str(),std::ios_base::app);
        double avg_ratio;
        double good_bins=0.0;
        double uncon;
        double con;
        for(int i=0;i<num_bins_sig;i++){
            std::cout<<"start_loop"<<std::endl;
            uncon=sqrt(m_covar_full_sys(i,i))/n_events[i];
            con=sqrt(m_covar_con_sys(i,i))/n_events[i];
            if (con>0 and uncon>0){
                avg_ratio+=uncon/con;
                good_bins+=1.0;
            }
            std::cout<<"bin="<<i<<"/events="<<n_events[i]<<"/uncon="<<uncon<<"/con="<<con<<" ratio: "<<uncon/con<<std::endl;
            constraint_table<<"bin="<<i<<"/events="<<n_events[i]<<"/uncon="<<uncon<<"/con="<<con<<" ratio: "<<uncon/con<<std::endl;
        }

        avg_ratio=avg_ratio/good_bins;
        //avg_ratio=0;

        std::vector<double> out_1bin={n_events[0], uncon, con, avg_ratio};
        TFile *fnew = new TFile(("constraint_test_"+tag+"_knob_"+variations[which_var]+".root").c_str(),"RECREATE");
        TMatrixD cov_before = coll_covariance.GetSub(0,num_bins_sig-1,0,num_bins_sig-1);
        TMatrixD cov_after = m_covar_con_sys.GetSub(0,num_bins_sig-1,0,num_bins_sig-1);
        //constraint_table.close();
        cov_before.Write("before");
        cov_after.Write("after");

        if(good_bins>0.0){
            return out_1bin;
        }
        else{
            std::vector<double> bad_out={-999,-999,-999,-999};
            return bad_out;
        }
    }



    std::vector<std::string> SBNcovariance::buildWeightMaps(){

        std::cout<<"SBNcovariance::buildWeightMaps()\t\t||\t\t Starting to Build Weight Maps of TTreeFormulas "<<std::endl;
        int n_wei = weightmaps_patterns.size();

        //this is what to return, this to be made into a TTreeFormula (for every single file!)
        std::vector<std::string> variation_weight_formulas(variations.size(),"1");

        for(int i=0; i< n_wei; i++){

            for(int v=0; v< variations.size(); v++){

                // Check to see if pattern is in this variation
                if (variations[v].find(weightmaps_patterns[i]) != std::string::npos) {
                    std::cout << "Variation "<<variations[v]<<" is a match for pattern "<<weightmaps_patterns[i]<<std::endl;
                    variation_weight_formulas[v] = variation_weight_formulas[v] + "*(" + weightmaps_formulas[i]+")";
                    std::cout<<" -- weight is thus "<<variation_weight_formulas[v]<<std::endl;
                    std::cout<<" -- mode is "<<weightmaps_mode[i]<<std::endl;
                    if(weightmaps_mode[i]=="multisim") m_variation_modes[v] = 0;
                    if(weightmaps_mode[i]=="minmax") m_variation_modes[v] = 1;

                }
            }
        }

        return variation_weight_formulas;
    }


