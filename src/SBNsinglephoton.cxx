#include "SBNsinglephoton.h"

using namespace sbn;

SBNsinglephoton::SBNsinglephoton(std::string xmlname, std::string intag):SBNconfig(xmlname), tag(intag){
    m_bool_cv_spectrum_generated = false;
    m_bool_cv_spectrum_loaded = false;
    m_bool_data_spectrum_loaded = false;

    m_full_fractional_covariance_matrix = nullptr;
    m_full_but_genie_fractional_covariance_matrix = nullptr;
    m_genie_fractional_covariance_matrix = nullptr;

    //initialize m_chi
    m_chi = new SBNchi(this->xmlname);
    m_chi->is_stat_only = false;  //do include MC intrinsic error   
}

SBNsinglephoton::SBNsinglephoton(std::string xmlname, std::string intag, NGrid ingrid):SBNsinglephoton(xmlname, intag, ingrid,ingrid, false){}

SBNsinglephoton::SBNsinglephoton(std::string xmlname, std::string intag, NGrid ingrid, NGrid in_polygrid, bool has_polygrid): SBNconfig(xmlname), tag(intag), m_grid(ingrid), m_bool_poly_grid(has_polygrid){

    if(is_verbose) std::cout << "SBNsinglephoton::SBNsinglephoton\t|| Setup grid" << std::endl;
    // setup flat grid
    m_vec_grid = m_grid.GetGrid();
    m_fit_dimension = m_grid.f_num_dimensions;
    m_flat_total_gridpoints = m_grid.f_num_total_points;

    if(is_verbose){
	std::cout << "Flat Grid:" << std::endl;
    	for(auto const & igrid: m_grid.f_dimensions)
	    std::cout << "\t" << igrid.GetName() << ": "<< igrid.GetNPoints()<< " points, ("<<igrid.GetMin() << ", " << igrid.GetMax()<<"), stepsize " << igrid.GetStep() << std::endl; 
    }

    // set up poly grid
    m_poly_total_gridpoints = 1;
    if(m_bool_poly_grid){
	SetPolyGrid(in_polygrid);
    }

    m_total_gridpoints = m_flat_total_gridpoints*m_poly_total_gridpoints;
    m_max_number_iterations = 20;
    m_chi_min_convergance_tolerance = 0.001;
    m_bool_modify_cv = false;
    m_bool_cv_spectrum_generated = false;
    m_bool_cv_spectrum_loaded = false;
    m_bool_data_spectrum_loaded = false;
    m_interpolation_number=100;

    m_full_fractional_covariance_matrix = nullptr;
    m_full_but_genie_fractional_covariance_matrix = nullptr;
    m_genie_fractional_covariance_matrix = nullptr;

    m_file_open=false;

    //initialize m_chi
    m_chi = new SBNchi(this->xmlname);
    //m_chi->is_stat_only = true; //not include MC intrinsic error, to reproduce NCpi0 fit result in technote v6.
    m_chi->is_stat_only = false;  //do include MC intrinsic error
    LocateCVGlobalIndex();
}

int SBNsinglephoton::SetPolyGrid(NGrid in_polygrid){
    if(is_verbose) std::cout << "SBNsinglephoton::SetPolyGrid\t|| Setup polynomial grid " << std::endl;
    m_bool_poly_grid = true;
    m_poly_grid=in_polygrid;
    m_vec_poly_grid = m_poly_grid.GetGrid();
    m_fit_dimension += m_poly_grid.f_num_dimensions;
    m_poly_total_gridpoints = m_poly_grid.f_num_total_points;

    m_total_gridpoints = m_flat_total_gridpoints*m_poly_total_gridpoints;

    if(is_verbose){
	std::cout << "Poly Grid:" << std::endl;
    	for(const auto &pgrid : m_poly_grid.f_dimensions)
            std::cout << "\t" << pgrid.GetName() << ": "<< pgrid.GetNPoints()<< " points, ("<< pgrid.GetMin() << ", " << pgrid.GetMax()<<"), stepsize " << pgrid.GetStep() << std::endl;
    }
    LocateCVGlobalIndex();
    return 0;
}

int SBNsinglephoton::OpenFiles(){

    otag = "SBNsinglephoton::OpenFiles\t|| ";
    num_files = montecarlo_file.size();
    montecarlo_additional_weight.resize(num_files,1.0);
    montecarlo_additional_weight_formulas.resize(num_files);


    if(is_verbose) std::cout<< otag << "Open the files"<< std::endl;
    for(auto &fn: montecarlo_file){
        files.push_back(new TFile(fn.c_str()));
        if(files.back()->IsZombie() || !files.back()->IsOpen()){
            std::cout<<otag << "ERROR! Failed to open the file "<<fn<<std::endl;
            exit(EXIT_FAILURE);
        }
    }


    for(int i=0; i<montecarlo_name.size(); i++){
        std::cout<<"Getting TTree "<<montecarlo_name[i]<<" from file "<<montecarlo_file[i]<<std::endl;
        trees.push_back((TTree*)files.at(i)->Get(montecarlo_name.at(i).c_str()) );
        std::cout<<"--TTree has "<<trees.back()->GetEntries()<<" entries. "<<std::endl;
    }

    if(is_verbose) std::cout << otag <<  "Setup any friendtrees" << std::endl;
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

    for(auto &t: trees){
        nentries.push_back(t->GetEntries());
    }


    for(int i=0; i< num_files; i++){

        double pot_scale = 1.0;
        if(montecarlo_pot[i]!=-1){
            pot_scale = this->plot_pot/montecarlo_pot[i];
        }

        montecarlo_scale[i] = montecarlo_scale[i]*pot_scale;

        std::cout << " TFile::Open() file=" << files[i]->GetName() << " @" << files[i] << std::endl;
        std::cout << " Has POT " <<montecarlo_pot[i] <<" and "<<nentries[i] <<" entries "<<std::endl;

        for(int k=0; k<branch_variables.at(i).size(); k++){
            const auto branch_variable = branch_variables.at(i).at(k);
            std::cout<<"Setting Branch: "<<branch_variable->name<<std::endl;
            branch_variable->branch_formula =  new TTreeFormula(("branch_form"+std::to_string(i)).c_str(), branch_variable->name.c_str(), trees[i]);

            if(branch_variable->GetOscillate()){
                std::cout<<"Setting true branch variables"<<std::endl;
                //if true_param_name is given as a branch name
                //trees.at(i)->SetBranchAddress( branch_variable->true_param_name.c_str(), branch_variable->GetTrueValue() );

                //if true_param_name is given in the form of a formula
                branch_variable->branch_true_value_formula = new TTreeFormula(("branch_true_value_form"+std::to_string(i)).c_str(),branch_variable->true_param_name.c_str(), trees[i]);
            }
        }

        if(montecarlo_additional_weight_bool[i]){
            std::cout<<"Setting Additional weight of : "<< montecarlo_additional_weight_names[i].c_str()<<std::endl;
            montecarlo_additional_weight_formulas[i] =  new TTreeFormula(("a_w"+std::to_string(i)).c_str(),montecarlo_additional_weight_names[i].c_str(),trees[i]);
        }


    }

    m_file_open= true;
    if(is_verbose) std::cout << otag << "Finish opening files and setting up TTrees " << std::endl; 

    return 0;
}


int SBNsinglephoton::CloseFiles(){
    if(m_file_open){
	    if(is_verbose) std::cout<< "SBNsinglephoton::CloseFiles\t|| Closing TFiles..."<< std::endl;
	    for(auto f: files){
		std::cout <<" TFile::Close() file=" << f->GetName() << " @" << f << std::endl;
		f->Close();
	    }
	    m_file_open=false;
    }
    return 0;
}


int SBNsinglephoton::ScaleSpectrum(SBNspec* inspec, double flat_factor, std::vector<double>& param){

    otag="SBNsinglephoton::ScaleSpectrum\t|| ";
    //if(is_verbose) std::cout << otag <<"-----------------------------------------------\n";

    for(int j=0;j<num_files;j++){


        for(int i=0; i< std::min(  montecarlo_maxevents.at(j)  ,nentries.at(j)); i++){
            trees.at(j)->GetEntry(i);

            //if(i%1000==0) std::cout<<otag << "On event: "<<i<<" of "<<nentries[j]<<" from File: "<<montecarlo_file[j]<<std::endl;

            double global_weight = 1.0;
            if( montecarlo_additional_weight_bool[j]){
                    montecarlo_additional_weight_formulas[j]->GetNdata();
                    global_weight = montecarlo_additional_weight_formulas[j]->EvalInstance();
            };//this will be 1.0 unless specified
            global_weight = global_weight*montecarlo_scale[j];



            if(std::isinf(global_weight) || global_weight != global_weight){
                std::cout<<otag << "ERROR  error @ "<<i<<" in File "<<montecarlo_file.at(j)<<" as its either inf/nan: "<<global_weight<<std::endl;
                exit(EXIT_FAILURE);
            }


            for(int t=0; t<branch_variables[j].size();t++){
                    const auto branch_variable = branch_variables[j][t];
                    int ih = inspec->map_hist.at(branch_variable->associated_hist);

		    //grab reco info and determine which reco bin it belongs to
                    //branch_variable->GetFormula()->GetNdata();
                    //double reco_var = branch_variable->GetFormula()->EvalInstance();
		    double reco_var = *(static_cast<double*>(branch_variable->GetValue()));
                    int reco_bin = inspec->GetGlobalBinNumber(reco_var,ih);

                    if(branch_variables[j][t]->GetOscillate()){
			//truth info
  			double true_var = *(static_cast<double*>(branch_variables[j][t]->GetTrueValue()));

                        double prescale_factor = this->ScaleFactor(true_var, flat_factor, param);
			MaskScaleFactor(prescale_factor); // do not allow negative weight
                        inspec->hist[ih].Fill(reco_var, global_weight*prescale_factor);
 		    }
            }
        } //end of entry loop
    } // end of nfile loop

    inspec->CalcFullVector();
    return 0;
}

int SBNsinglephoton::PreScaleSpectrum(std::string xmlname, double flat_factor, std::vector<double>& param){
    otag="SBNsinglephoton::PreScaleSpectrum\t|| ";
    SBNspec tm(xmlname,-1,false);
    SBNspec temp_cv_spectrum = tm;
    SBNspec spec_prescale  = tm;   //prescaled spectrum

    std::cout<<otag << "-----------------------------------------------\n";
    std::cout<<otag << "-----------------------------------------------\n";

    for(int j=0;j<num_files;j++){


        for(int i=0; i< std::min(  montecarlo_maxevents.at(j)  ,nentries.at(j)); i++){
            trees.at(j)->GetEntry(i);

            if(i%100==0) std::cout<<otag << "On event: "<<i<<" of "<<nentries[j]<<" from File: "<<montecarlo_file[j]<<std::endl;

            double global_weight = 1.0;
            if( montecarlo_additional_weight_bool[j]){
                montecarlo_additional_weight_formulas[j]->GetNdata();
                global_weight = montecarlo_additional_weight_formulas[j]->EvalInstance();
            };//this will be 1.0 unless specified
            global_weight = global_weight*montecarlo_scale[j];



            if(std::isinf(global_weight) || global_weight != global_weight){
                std::cout<<otag<<"ERROR  error @ "<<i<<" in File "<<montecarlo_file.at(j)<<" as its either inf/nan: "<<global_weight<<std::endl;
                exit(EXIT_FAILURE);
            }



            for(int t=0; t<branch_variables[j].size();t++){
                const auto branch_variable = branch_variables[j][t];
                int ih = temp_cv_spectrum.map_hist.at(branch_variable->associated_hist);

                //grab reco info and determine which reco bin it belongs to
                //branch_variable->GetFormula()->GetNdata();
                //double reco_var = branch_variable->GetFormula()->EvalInstance();
		double reco_var = *(static_cast<double*>(branch_variable->GetValue()));
                int reco_bin = temp_cv_spectrum.GetGlobalBinNumber(reco_var,ih);

                if(branch_variables[j][t]->GetOscillate()){
                    //truth info
                    double true_var = *(static_cast<double*>(branch_variables[j][t]->GetTrueValue()));

                    double prescale_factor = this->ScaleFactor(true_var, flat_factor, param);
		    MaskScaleFactor(prescale_factor); // do not allow negative weight

                    spec_prescale.hist[ih].Fill(reco_var, global_weight*prescale_factor);
                    if(!m_bool_cv_spectrum_generated)  temp_cv_spectrum.hist[ih].Fill(reco_var,global_weight);
                }else{
                    if(!m_bool_cv_spectrum_generated) temp_cv_spectrum.hist[ih].Fill(reco_var,global_weight);
                }
            }
        } //end of entry loop
    } // end of nfile loop


    if(is_verbose) std::cout<< otag << "Write out spectra" << std::endl;
    // write out prescaled spectra
    std::ostringstream prescale_tag;
    if(flat_factor != 0) prescale_tag << "Flat_" << std::fixed<< std::setprecision(3) << flat_factor << "_PreScaled";
    else prescale_tag << "PreScaled";
    //generate tag for prescaled spectra
    for(int i=0; i< param.size(); i++){
        prescale_tag << "_" << std::fixed<< std::setprecision(3) << param[i]; 
    }   
    if(param.size() != 0) spec_prescale.WriteOut(tag+"_"+prescale_tag.str());

    // write out CV
    if(!m_bool_cv_spectrum_generated){
        m_bool_cv_spectrum_generated = true;
        temp_cv_spectrum.WriteOut(tag+"_CV");
    }
    return 0;
}


int SBNsinglephoton::PreScaleSpectrum(std::string xmlname, std::vector<double>& param){
    return this->PreScaleSpectrum(xmlname, 0.0, param);
}


int SBNsinglephoton::GeneratePreScaledSpectra(){
    otag="SBNsinglephoton::GeneratePreScaledSpectra\t|| ";
    if(!m_file_open) this->OpenFiles();
    if(!m_bool_poly_grid){
        std::cout << otag<< "No grid for polynomial scaling present------------------\n"<< std::endl;
        std::cout << otag<< "Gonna generate only CV, then exiting GeneratePreScaledSpectra ..." << std::endl;
        std::vector<double> empty{};
        this->PreScaleSpectrum(xmlname, empty);
    }else{

	std::cerr << otag << "This is obsolete, please use LoadSpectraApplyFullScaling() funtion to generate scaled SBNspec! " << std::endl;
	std::cerr << otag << "Exiting... " << std::endl;
	exit(EXIT_FAILURE);
 
        for(int i=0; i< m_poly_total_gridpoints; i++){
            std::vector<double> point = m_vec_poly_grid[i];
            this->PreScaleSpectrum(xmlname, point);	
        }
    }

    this->CloseFiles();

    return 0;
}

// this function is safe to use when we have poly grid
// it handles MC intrinsic error correctly
int SBNsinglephoton::LoadSpectraOnTheFly(){
	otag="SBNsinglephoton::LoadSpectraOnTheFly\t|| ";
	if(is_verbose) std::cout << otag << "Calculate and Load scaled spectra on the fly" << std::endl;

	if(!m_file_open) this->OpenFiles();
	int flat_grid_index=-1;
	int flat_grid_npoint=0;
	int period=0;
	for(int i=0; i<m_grid.f_num_dimensions; i++){
	    period *= m_grid.f_dimensions[i].GetNPoints();
	    if(m_grid.f_dimensions[i].GetName() == m_poly_grid.f_dimensions[0].GetName()){
		 flat_grid_index=i;
		 flat_grid_npoint=m_grid.f_dimensions[i].GetNPoints();
		 period =1;
	    }
	}

	//calculate the flat+poly part
	int vec_index=0;
	if(is_verbose) std::cout << otag << "Generating flat+poly part of the spectra"<< std::endl;
	SBNspec tm(xmlname, -1, false);
	std::vector<SBNspec> m_prescale_spec(m_poly_total_gridpoints*flat_grid_npoint, tm);	
	for(int i=0; i<m_poly_total_gridpoints; i++){
	    std::vector<double> ipoint = m_vec_poly_grid[i];
	    std::vector<double> fgrid_point = m_grid.f_dimensions[flat_grid_index].GetPoints();
	    for(int j=0; j< flat_grid_npoint; j++){
		if(is_verbose && (vec_index%50==0)) std::cout << "On Point " << vec_index << std::endl;
		this->ScaleSpectrum(&m_prescale_spec[vec_index], fgrid_point[j], ipoint);
		++vec_index;
	    }
	}

	vec_index=0;
	//get the total scaled spectra
	if(is_verbose) std::cout << otag << "Grab the whole scaled spectra " <<std::endl;
	m_scaled_spec_grid.resize(m_total_gridpoints);
	for(int i=0; i<m_poly_total_gridpoints; i++){
	    for(int j=0; j<m_flat_total_gridpoints; j++){

		if(vec_index% 1000 == 0) std::cout << "On Point "<< vec_index << std::endl;
		
		std::vector<double> jpoint = m_vec_grid[j];

		m_scaled_spec_grid[vec_index]=*m_cv_spectrum;
		for(int k=0; k<jpoint.size(); k++){
		      if(k == flat_grid_index) m_scaled_spec_grid[vec_index].Scale(m_grid.f_dimensions[k].GetName(),0.0);
		      else m_scaled_spec_grid[vec_index].Scale(m_grid.f_dimensions[k].GetName(), jpoint[k]);
		}

		int flat_index = (j%(period*flat_grid_npoint))/period; 
		m_scaled_spec_grid[vec_index].Add(&m_prescale_spec[i*flat_grid_npoint+flat_index]);


 	        if( m_vec_poly_grid[i][0] == -0.2 && jpoint[0] == 0.8 && jpoint[1] == 2 ){
 	        //if( m_vec_poly_grid[i][0] == -0.2 && jpoint == std::vector<double>({0.8, 2})){
		    std::cout << "WriteOut: -0.2, 0.8, 2.0 " << std::endl;
		    m_scaled_spec_grid[vec_index].WriteOut(tag+"_sample_0.8_2.0_neg_0.2");
		}

		++vec_index;
	   }
	}
	return 0;
}



int SBNsinglephoton::LoadSpectraApplyFullScaling(){
    otag="SBNsinglephoton::LoadSpectraApplyFullScaling\t|| ";
    if(!m_bool_poly_grid){

        //SBNspec temp_cv = *m_cv_spectrum;
        SBNspec temp_cv(tag+"_CV.SBNspec.root", this->xmlname, false);

        m_scaled_spec_grid.resize(m_total_gridpoints); 

    
        std::cout<<otag<<"Scale the spectra for the whole grid " << std::endl;
        for(int j=0; j< m_flat_total_gridpoints;++j){
            if(is_verbose && (j%1000 ==0)) std::cout << "On Point " << j <<"/"<<m_total_gridpoints << std::endl;

            std::vector<double> jpoint = m_vec_grid[j];

            m_scaled_spec_grid[j] = temp_cv; //start with genie CV
            for(int jp=0; jp< jpoint.size(); jp++){
                m_scaled_spec_grid[j].Scale(m_grid.f_dimensions[jp].GetName(), jpoint[jp]); // scale corresponding subchannel
                //m_scaled_spec_grid[j].ScaleAll(jpoint[jp]);  
            }

	    //m_scaled_spec_grid[j].Scale("NCDeltaRadOverlayLEE", 0.0); //to reproduce technote v6.0 NCpi0 fit result
        }
    }else{
	std::cout << otag << "Involves poly grid, should load spectra on the fly" << std::endl;
        LoadSpectraOnTheFly(); 
    }

    return 0;
}



int SBNsinglephoton::CalcChiGridScanShapeOnlyFit(){
    otag="SBNsinglephoton::CalcChiGridScanShapeOnlyFit\t|| ";
    CheckCVLoad();
    CheckDataLoad();
    CheckNumScaledSpectra();

    SBNspec background_spectrum;
    SBNspec last_best_spectrum;
    double best_chi, last_best_chi;
    int best_point, last_best_point;
    std::vector<double> vec_chi, last_vec_chi;
    TFile* fout = new TFile(Form("%s_fit_output.root", tag.c_str()), "recreate");  //save matrix plot etc.

    m_full_but_genie_fractional_covariance_matrix->Write("frac_but_genie");
    m_genie_fractional_covariance_matrix->Write("frac_genie"); 

    TMatrixT<double> full_systematic_covariance(num_bins_total, num_bins_total);
    TMatrixT<double> genie_systematic_matrix(num_bins_total, num_bins_total);
    TMatrixT<double> collapsed_full_systematic_matrix(num_bins_total_compressed, num_bins_total_compressed);
    TMatrixT<double> total_covariance_matrix(num_bins_total_compressed, num_bins_total_compressed);
    TMatrixT<double> inversed_total_covariance_matrix(num_bins_total_compressed, num_bins_total_compressed);

    //m_max_number_iterations = 5;
    for(int n_iter =0; n_iter < m_max_number_iterations; n_iter ++){
        std::cout << otag << "On fit iteration "<< n_iter << std::endl;
        //reset everything at the biginning of each iteration
        best_chi = DBL_MAX;
        vec_chi.clear();


        if(n_iter == 0){
            last_best_spectrum = *m_cv_spectrum;
            last_best_spectrum.Scale("NCDeltaRadOverlayLEE", 0.0);  //to reproduce NCpi0 fit result in technote v6.0 
        }else{
            last_best_spectrum = m_scaled_spec_grid[last_best_point];		
        }

	last_best_spectrum.CalcErrorVector();
        background_spectrum = last_best_spectrum;

        //============================Calculate covariance matrix and its invert====================================/

       //The Full except for genie. We will be adding select bins on later
        full_systematic_covariance = m_chi->FillSystMatrix(*m_full_but_genie_fractional_covariance_matrix, last_best_spectrum.full_vector, last_best_spectrum.full_err_vector);

	TMatrixT<double> full_shapeonly_matrix = m_chi->CalcShapeOnlyCovarianceMatrix(*m_genie_fractional_covariance_matrix,&last_best_spectrum,  &last_best_spectrum);
	TMatrixT<double> full_shapemix_matrix = m_chi->CalcShapeMixedCovarianceMatrix(*m_genie_fractional_covariance_matrix, &last_best_spectrum, &last_best_spectrum);

	for(int i=0 ; i<num_bins_total; i++){
		for(int j =0; j<num_bins_total; j++){
			full_shapeonly_matrix(i,j) /= last_best_spectrum.full_vector.at(i)*last_best_spectrum.full_vector.at(j);
			full_shapemix_matrix(i,j) /= last_best_spectrum.full_vector.at(i)*last_best_spectrum.full_vector.at(j);
		}
	}
	fout->cd();
	full_shapeonly_matrix.Write("full_genie_shape_only_systematic_matrix");
	full_shapemix_matrix.Write("full_genie_shape_mixed_systematic_matrix");

        //calculate the shape only covariance matrix for genie uncertainty, to get rid of normalization uncertainty
        for(int i=0; i<m_grid.f_num_dimensions; i++){
            //if(m_grid.f_dimensions[i].GetName() == "NCDeltaLEE") continue;
            SBNspec temp_comp = last_best_spectrum;
            temp_comp.Keep(m_grid.f_dimensions[i].GetName(), 1.0);

            //only genie for the 3 scaled channels
            genie_systematic_matrix = m_chi->FillSystMatrix(*m_genie_fractional_covariance_matrix, temp_comp.full_vector, temp_comp.full_err_vector);
            fout->cd();
            genie_systematic_matrix.Write(Form("full_genie_%s_%d", m_grid.f_dimensions[i].GetName().c_str(), n_iter));	
            
            //This is then the either Shape or Shape plus Mixed genie for 3 scaled channels
            //genie_systematic_matrix = m_chi->CalcShapeOnlyCovarianceMatrix(*m_genie_fractional_covariance_matrix, &temp_comp, &temp_comp);
            genie_systematic_matrix = m_chi->CalcShapeMixedCovarianceMatrix(*m_genie_fractional_covariance_matrix, &temp_comp, &temp_comp);
            genie_systematic_matrix.Write(Form("ShapeMixed_genie_%s_%d", m_grid.f_dimensions[i].GetName().c_str(), n_iter));	
            full_systematic_covariance += genie_systematic_matrix;//this point its flux all+det all + genie only 3 param shape+mixed
            background_spectrum.Scale(m_grid.f_dimensions[i].GetName(), 0.0);   
        }
        //add genie uncertainties for other subchannels to total covariance matrix. Its being Re-Used here. 
        genie_systematic_matrix = m_chi->FillSystMatrix(*m_genie_fractional_covariance_matrix, background_spectrum.full_vector, background_spectrum.full_err_vector);
        full_systematic_covariance += genie_systematic_matrix;//This is all everything EXCEPT 3 param norm

        m_chi->CollapseModes(full_systematic_covariance, collapsed_full_systematic_matrix);

        //use CNP statistical covariance matrix
        total_covariance_matrix = m_chi->AddStatMatrixCNP(&collapsed_full_systematic_matrix, last_best_spectrum.collapsed_vector, m_data_spectrum->collapsed_vector);
        fout->cd();
        genie_systematic_matrix.Write(Form("full_genie_otherbkg_%d", n_iter));
        full_systematic_covariance.Write(Form("syst_uncollapsed_matrix_%d", n_iter));
        collapsed_full_systematic_matrix.Write(Form("syst_collapsed_matrix_%d", n_iter));
        total_covariance_matrix.Write(Form("total_collapsed_matrix_%d", n_iter));

        inversed_total_covariance_matrix= m_chi->InvertMatrix(total_covariance_matrix);


        if(is_verbose && m_bool_data_spectrum_loaded && n_iter ==0) m_chi->DrawComparisonIndividual(last_best_spectrum, *m_data_spectrum, collapsed_full_systematic_matrix, tag+"_CVvsData", true);	
        //============================Done calculating covariance matrix ============================================/


        for(int i=0 ;i <m_total_gridpoints; i++){
            if(i%1000 == 0) std::cout<< "On Point " << i << "/" << m_total_gridpoints << std::endl;
            double temp_chi = m_chi->CalcChi(inversed_total_covariance_matrix, m_scaled_spec_grid[i].collapsed_vector, m_data_spectrum->collapsed_vector);
            vec_chi.push_back(temp_chi);

            if(temp_chi < best_chi){
                best_chi = temp_chi;
                best_point = i;
            }
        }

        if(is_verbose) std::cout << otag << "chi2 minimum :" << best_chi << " at point " << best_point << "/" << m_total_gridpoints << std::endl;

        if(n_iter != 0){
            if(fabs(best_chi - last_best_chi) < m_chi_min_convergance_tolerance){
                std::cout << otag << "chi2 has converged with best chi2 value " << best_chi << " at point " << best_point << "/" << m_total_gridpoints << std::endl;
                break;
            }
        }

        last_best_chi = best_chi;
        last_best_point = best_point;
        last_vec_chi = vec_chi;

	if(n_iter ==0){
		std::map<int, std::vector<double>> temp_map={{best_point, vec_chi}};
		this->PrintOutFitInfo(temp_map, "SBNsinglephoton::CalcChiGridScanShapeOnlyFit\t CV||"+tag, false);
	}
    }//end loop for iteration

    fout->Close();
    if(is_verbose ){
        //SBNspec temp_best_spec = this->GeneratePointSpectra(best_point);
	//m_chi->DrawComparisonIndividual(temp_best_spec, *m_data_spectrum, collapsed_full_systematic_matrix, tag+"_BFvsData");
	m_chi->DrawComparisonIndividual(m_scaled_spec_grid[best_point], *m_data_spectrum, collapsed_full_systematic_matrix, tag+"_BFvsData");
    }
	
    m_map={{best_point, vec_chi}};
    this->PrintOutFitInfo(m_map, otag+"BF - "+tag, true);
    this->WriteOutInfo(m_map);
    this->CloseFiles();
    return 0;
}



int SBNsinglephoton::CalcChiGridScan(){
    otag="SBNsinglephoton::CalcChiGridScan\t|| ";
    CheckCVLoad();
    CheckDataLoad();
    CheckNumScaledSpectra();

    SBNspec last_best_spectrum;
    double best_chi, last_best_chi;
    int best_point, last_best_point;
    std::vector<double> vec_chi, last_vec_chi;
    TFile* fout = new TFile(Form("%s_fit_output.root", tag.c_str()), "recreate");  //save matrix plot etc.

    m_full_fractional_covariance_matrix->Write("full_fractional_matrix");

    TMatrixT<double> collapsed_full_systematic_matrix(num_bins_total_compressed, num_bins_total_compressed);
    TMatrixT<double> total_covariance_matrix(num_bins_total_compressed, num_bins_total_compressed);
    TMatrixT<double> inversed_total_covariance_matrix(num_bins_total_compressed, num_bins_total_compressed);

    TMatrixT<double>* collapse_frac;

    for(int n_iter =0; n_iter < m_max_number_iterations; n_iter ++){
        std::cout << otag << "On fit iteration "<< n_iter << std::endl;
        //reset everything at the biginning of each iteration
        best_chi = DBL_MAX;
        vec_chi.clear();

        if(n_iter == 0){
            last_best_spectrum = *m_cv_spectrum;
            last_best_spectrum.Scale("NCDeltaLEE", 0.0);
        }else{
            last_best_spectrum = m_scaled_spec_grid[last_best_point];		
        }

	last_best_spectrum.CalcErrorVector();

        //============================Calculate covariance matrix and its invert====================================/
        //full systematic covariance matrix, except genie uncertainty

        //for(int i=0; i<last_best_spectrum.full_vector.size(); i++)
        //	std::cout << last_best_spectrum.full_vector[i] << ",";
        //std::cout << std::endl;

        //SBNspec temp_signal = last_best_spectrum; temp_signal.Keep("NCDelta", 1.0);
        //SBNspec temp_other = last_best_spectrum; temp_other.Scale("NCDelta", 0.0);
        //collapsed_full_systematic_matrix = m_chi->FillSystMatrix(*m_full_fractional_covariance_matrix, temp_signal.full_vector, true) + m_chi->FillSystMatrix(*m_full_fractional_covariance_matrix, temp_other.full_vector, true);
	//TMatrixT<double> uncollapsed_full_systematic_matrix = m_chi->FillSystMatrix(*m_full_fractional_covariance_matrix, last_best_spectrum.full_vector,false);

	// since we are using intrinsic error here, need to make sure m_scaled_spec_grid is setup correctly
        collapsed_full_systematic_matrix = m_chi->FillSystMatrix(*m_full_fractional_covariance_matrix, last_best_spectrum.full_vector, last_best_spectrum.full_err_vector, true);
        collapse_frac = (TMatrixT<double>*)collapsed_full_systematic_matrix.Clone();
        for(int i=0; i<num_bins_total_compressed ;i ++){
            for(int j=0 ; j< num_bins_total_compressed; j++)
                (*collapse_frac)(i,j) /= last_best_spectrum.collapsed_vector[i]* last_best_spectrum.collapsed_vector[j];
        }

        //use CNP statistical covariance matrix
        //total_covariance_matrix = m_chi->AddStatMatrix(&collapsed_full_systematic_matrix, m_data_spectrum->collapsed_vector);
        total_covariance_matrix = m_chi->AddStatMatrixCNP(&collapsed_full_systematic_matrix, last_best_spectrum.collapsed_vector, m_data_spectrum->collapsed_vector);

        inversed_total_covariance_matrix= m_chi->InvertMatrix(total_covariance_matrix);

        fout->cd();
	//uncollapsed_full_systematic_matrix.Write(Form("syst_uncollapsed_matrix_%d", n_iter));
        collapsed_full_systematic_matrix.Write(Form("syst_collapsed_matrix_%d", n_iter));
        total_covariance_matrix.Write(Form("total_collapsed_matrix_%d", n_iter));
        collapse_frac->Write(Form("collapsed_fractional_matrix_%d", n_iter));
        inversed_total_covariance_matrix.Write(Form("inversed_matrix_%d", n_iter));

	//print out comparison
        if(is_verbose && m_bool_data_spectrum_loaded && (n_iter==0)){
	    m_chi->DrawComparisonIndividual(last_best_spectrum, *m_data_spectrum, collapsed_full_systematic_matrix, tag+"_CVvsData", true);
	}
        //============================Done calculating covariance matrix ============================================/


        for(int i=0 ;i <m_total_gridpoints; i++){
            if(i%1000 ==0) std::cout<< "On Point " << i << "/" << m_total_gridpoints << std::endl;
            double temp_chi;
            temp_chi = m_chi->CalcChi(inversed_total_covariance_matrix, m_scaled_spec_grid[i].collapsed_vector, m_data_spectrum->collapsed_vector, false);
            vec_chi.push_back(temp_chi);

            if(temp_chi < best_chi){
                best_chi = temp_chi;
                best_point = i;
            }
        }

        if(is_verbose) std::cout << otag << "chi2 minimum :" << best_chi << " at point " << best_point << "/" << m_total_gridpoints << std::endl;

        if(n_iter != 0){
            if(fabs(best_chi - last_best_chi) < m_chi_min_convergance_tolerance){
                std::cout << otag << "chi2 has converged with best chi2 value " << best_chi << " at point " << best_point << "/" << m_total_gridpoints << std::endl;
                break;
            }
        }

        last_best_chi = best_chi;
        last_best_point = best_point;
        last_vec_chi = vec_chi;

	if(n_iter ==0){
		std::map<int, std::vector<double>> temp_map={{best_point, vec_chi}};
		this->PrintOutFitInfo(temp_map, otag+ "CV - "+tag, false);
	}
    }//end loop for iteration

    fout->Close();

    //best-fit vs data comparison	
    if(is_verbose && m_bool_data_spectrum_loaded){
    //if(is_verbose ){
        //SBNspec temp_best_spec = this->GeneratePointSpectra(best_point);
	//m_chi->DrawComparisonIndividualFracMatrix(temp_best_spec, *m_data_spectrum, *m_full_fractional_covariance_matrix, tag+"_BFvsData", true);
	m_chi->DrawComparisonIndividualFracMatrix(m_scaled_spec_grid[best_point], *m_data_spectrum, *m_full_fractional_covariance_matrix, tag+"_BFvsData", true);
    }

    std::cout << " BEST-FIT \t\t\t Data \t\t\t CV " << std::endl;
    for(int i=0 ; i< m_data_spectrum->collapsed_vector.size(); i++){
        std::cout <<i << ": " << m_scaled_spec_grid[best_point].collapsed_vector[i] << "\t\t\t " << m_data_spectrum->collapsed_vector[i] << " \t\t\t " << m_cv_spectrum->collapsed_vector[i] << std::endl;
    }
    m_map={{best_point, vec_chi}};
    this->PrintOutFitInfo(m_map, otag+"BF - "+tag, true);
    this->WriteOutInfo(m_map);
    this->CloseFiles();
    return 0;
}



int SBNsinglephoton::CalcChiGridScanVaryMatrices(){
    CheckCVLoad();
    CheckDataLoad();
    CheckNumScaledSpectra();

    otag="SBNsinglephoton::CalcChiGridScanVaryMatrices\t|| ";
    std::cout << otag<<  "Start grid scan, with covariance matrix varied at each grid point" << std::endl; 

    double best_chi = DBL_MAX;
    int best_point;
    std::vector<double> vec_chi;
    TFile* fout = new TFile(Form("%s_fit_output.root", tag.c_str()), "recreate");  //save matrix plot etc.

    m_full_fractional_covariance_matrix->Write("full_fractional_matrix");

    TMatrixT<double> collapsed_full_systematic_matrix(num_bins_total_compressed, num_bins_total_compressed);
    TMatrixT<double> total_covariance_matrix(num_bins_total_compressed, num_bins_total_compressed);
    TMatrixT<double> inversed_total_covariance_matrix(num_bins_total_compressed, num_bins_total_compressed);

    for(int i =0; i!= m_total_gridpoints; ++i){
	if(i % 1000 == 0) std::cout << otag << "On Point " << i << "/" << m_total_gridpoints << std::endl;

	//calculate the covariance matrix at grid point i
	m_scaled_spec_grid[i].CollapseVector();
	m_scaled_spec_grid[i].CalcErrorVector();	
	collapsed_full_systematic_matrix = m_chi->FillSystMatrix(*m_full_fractional_covariance_matrix, m_scaled_spec_grid[i].full_vector, m_scaled_spec_grid[i].full_err_vector, true);	
        //total_covariance_matrix = m_chi->AddStatMatrix(&collapsed_full_systematic_matrix, m_scaled_spec_grid[i].collapsed_vector);	
        //total_covariance_matrix = m_chi->AddStatMatrix(&collapsed_full_systematic_matrix, m_data_spectrum->collapsed_vector);	
        total_covariance_matrix = m_chi->AddStatMatrixCNP(&collapsed_full_systematic_matrix, m_scaled_spec_grid[i].collapsed_vector, m_data_spectrum->collapsed_vector);
	inversed_total_covariance_matrix= m_chi->InvertMatrix(total_covariance_matrix);


	//save covariance matrix & spectra comparison at CV
	if(i == m_cv_spec_index){
            fout->cd();
            collapsed_full_systematic_matrix.Write("syst_collapsed_matrix_CV");
            total_covariance_matrix.Write("total_collapsed_matrix_CV");
            inversed_total_covariance_matrix.Write("inversed_matrix_CV");

	    //print out comparison between data and CV
	    if(is_verbose && m_bool_data_spectrum_loaded)
		m_chi->DrawComparisonIndividual(m_scaled_spec_grid[i], *m_data_spectrum, collapsed_full_systematic_matrix, tag+"_CVvsData", true);

	    //print collapsed fractional covaraicne matrix at CV
	    PrintCollapsedFractionalMatrix(*m_full_fractional_covariance_matrix, m_scaled_spec_grid[i]);

	    //print collapsed vector of CV and data
	    std::cout << "Collapsed Vector of GENIE CV: ";
	    m_scaled_spec_grid[i].PrintCollapsedVector();
	    std::cout << "Collapsed Vector of Data    : ";
	    m_data_spectrum->PrintCollapsedVector();
	}
	
	//no determinant
        vec_chi.push_back(m_chi->CalcChi(inversed_total_covariance_matrix, m_scaled_spec_grid[i].collapsed_vector, m_data_spectrum->collapsed_vector));

	//add log of matrix determinant
        //vec_chi.push_back(m_chi->CalcChi(inversed_total_covariance_matrix, m_scaled_spec_grid[i].collapsed_vector, m_data_spectrum->collapsed_vector) + log(abs(total_covariance_matrix.Determinant())));

        if(vec_chi.back() < best_chi){
            best_chi = vec_chi.back();
            best_point = i;
        }

    } //loop over the grid
    fout->Close();

    std::cout << otag << "chi2 minimum value " << best_chi << " is found at point " << best_point << "/" << m_total_gridpoints << std::endl;
    //best-fit vs data comparison	
    if(is_verbose && m_bool_data_spectrum_loaded){
	m_chi->DrawComparisonIndividualFracMatrix(m_scaled_spec_grid[best_point], *m_data_spectrum, *m_full_fractional_covariance_matrix, tag+"_BFvsData", true);
    }


    m_map={{best_point, vec_chi}};
    this->PrintOutFitInfo(m_map,  otag + "BF - "+tag, true);
    this->WriteOutInfo(m_map);
    this->CloseFiles();
    return 0;
}


int SBNsinglephoton::LoadCV(){
    if(is_verbose) std::cout << "SBNsinglephoton::LoadCV\t|| Setup CV spectrum" << std::endl;
    m_cv_spectrum = new SBNspec(tag+"_CV.SBNspec.root", xmlname, false);
    m_cv_spectrum->CollapseVector();
    m_cv_spectrum->CalcErrorVector();
    m_bool_cv_spectrum_loaded = true;

    return 0;
}

int SBNsinglephoton::LoadData(std::string filename){
    if(is_verbose) std::cout << "SBNsinglephoton::LoadData\t|| Load data spectrum from file: " << filename << std::endl;
    m_data_spectrum = new SBNspec(filename, xmlname, false);
    m_bool_data_spectrum_loaded = true;
    m_data_spectrum->CollapseVector(); 
    return 0;
}


int SBNsinglephoton::SetFullFractionalCovarianceMatrix(std::string filename, std::string matrix_name){
	otag = "SBNsinglephoton::SetFullFractionalCovarianceMatrix\t|| ";
	if(is_verbose) std::cout << otag << "Open total frac covariance matrix file: " << filename << std::endl;
	TFile* f_syst = new TFile(filename.c_str(), "read");
	m_full_fractional_covariance_matrix = (TMatrixT<double>*)f_syst->Get(matrix_name.c_str());

	if(m_full_fractional_covariance_matrix->GetNcols() != m_full_fractional_covariance_matrix->GetNrows()){
		std::cout << otag<<  "Matrix provided is not sysmetric" << std::endl;
		exit(EXIT_FAILURE);		
	}


	this->RemoveNan(m_full_fractional_covariance_matrix);

	if(is_verbose) std::cout << otag << "matrix size: " << m_full_fractional_covariance_matrix->GetNcols()<< std::endl;;

/*
	TMatrixT<double>* temp_genie_XS = (TMatrixT<double>*)f_syst->Get("individualDir/All_UBGenie_frac_covariance");
        *temp_genie_XS += *((TMatrixT<double>*)f_syst->Get("individualDir/AxFFCCQEshape_UBGenie_frac_covariance"));
        *temp_genie_XS += *((TMatrixT<double>*)f_syst->Get("individualDir/DecayAngMEC_UBGenie_frac_covariance"));
        *temp_genie_XS += *((TMatrixT<double>*)f_syst->Get("individualDir/NormCCCOH_UBGenie_frac_covariance"));
        *temp_genie_XS += *((TMatrixT<double>*)f_syst->Get("individualDir/NormNCCOH_UBGenie_frac_covariance"));
        *temp_genie_XS += *((TMatrixT<double>*)f_syst->Get("individualDir/RPA_CCQE_UBGenie_frac_covariance"));
        *temp_genie_XS += *((TMatrixT<double>*)f_syst->Get("individualDir/Theta_Delta2Npi_UBGenie_frac_covariance"));
        *temp_genie_XS += *((TMatrixT<double>*)f_syst->Get("individualDir/VecFFCCQEshape_UBGenie_frac_covariance"));
        *temp_genie_XS += *((TMatrixT<double>*)f_syst->Get("individualDir/XSecShape_CCMEC_UBGenie_frac_covariance"));
        this->RemoveNan(temp_genie_XS);

	*m_full_fractional_covariance_matrix = *temp_genie_XS;

	for(int i=0; i<m_full_fractional_covariance_matrix->GetNcols(); i++){
	    for(int j=0; j<m_full_fractional_covariance_matrix->GetNcols(); j++){
		//if(std::isnan((*m_full_fractional_covariance_matrix)(i,j))) (*m_full_fractional_covariance_matrix)(i,j) =0.0;
		if( i>89 && j <= 89) (*m_full_fractional_covariance_matrix)(i,j) =0.0;
		if( i<=89 && j > 89) (*m_full_fractional_covariance_matrix)(i,j) =0.0;
		//if( (i== 2 || i==8) && j == 329) (*m_full_fractional_covariance_matrix)(i,j) =0.0;
		//if( (j== 2 || j==8) && i == 329) (*m_full_fractional_covariance_matrix)(i,j) =0.0;
	    }	
	}
*/
	f_syst->Close();
	return 0;
}

int SBNsinglephoton::SetGenieFractionalCovarianceMatrix(std::string filename){
 	if(is_verbose) std::cout << "SBNsinglephoton::SetGenieFractionalCovarianceMatrix\t||\t Opening Flux+XS covariance matrix file: " << filename<< std::endl;

	//check with Gray, are matrices added directly before checking the nan's?
	TFile* f_syst = new TFile(filename.c_str(), "read");

//        m_genie_fractional_covariance_matrix  = (TMatrixT<double>*)f_syst->Get("frac_covariance");

	//m_genie_fractional_covariance_matrix = new TMatrixT<double>(num_bins_total, num_bins_total);
        m_genie_fractional_covariance_matrix  = (TMatrixT<double>*)f_syst->Get("individualDir/All_UBGenie_frac_covariance");
        *m_genie_fractional_covariance_matrix += *((TMatrixT<double>*)f_syst->Get("individualDir/AxFFCCQEshape_UBGenie_frac_covariance"));
        *m_genie_fractional_covariance_matrix += *((TMatrixT<double>*)f_syst->Get("individualDir/DecayAngMEC_UBGenie_frac_covariance"));
        *m_genie_fractional_covariance_matrix += *((TMatrixT<double>*)f_syst->Get("individualDir/NormCCCOH_UBGenie_frac_covariance"));
        *m_genie_fractional_covariance_matrix += *((TMatrixT<double>*)f_syst->Get("individualDir/NormNCCOH_UBGenie_frac_covariance"));
        *m_genie_fractional_covariance_matrix += *((TMatrixT<double>*)f_syst->Get("individualDir/RPA_CCQE_UBGenie_frac_covariance"));
        *m_genie_fractional_covariance_matrix += *((TMatrixT<double>*)f_syst->Get("individualDir/Theta_Delta2Npi_UBGenie_frac_covariance"));
        *m_genie_fractional_covariance_matrix += *((TMatrixT<double>*)f_syst->Get("individualDir/VecFFCCQEshape_UBGenie_frac_covariance"));
	//XSecShape_CCMEC_UBGenie_frac_covariance is bugged, so not using it now
        *m_genie_fractional_covariance_matrix += *((TMatrixT<double>*)f_syst->Get("individualDir/XSecShape_CCMEC_UBGenie_frac_covariance"));
	this->RemoveNan(m_genie_fractional_covariance_matrix);


/*
	SBNspec temp_cv = *m_cv_spectrum;
	temp_cv.Keep("NCPi0Coh", 1.0);
	temp_cv.CalcFullVector();
	std::vector<double> temp_full = temp_cv.full_vector;


	TMatrixT<double> M_cohOnly = *m_genie_fractional_covariance_matrix;
	for(int i=0; i<num_bins_total; i++){
	   for(int j=0; j< num_bins_total; j++){
	     if(temp_full[i]==0 && temp_full[j]==0 ) M_cohOnly(i,j) = 0;
	   }
	}

	TMatrixT<double>* M_updatedCohOnly = (TMatrixT<double>*)f_syst->Get("individualDir/updated_allUBGenie_frac_covariance");
	for(int i=0; i<num_bins_total; i++){
           for(int j=0; j< num_bins_total; j++){
             if(temp_full[i]==0 && temp_full[j]==0 ) (*M_updatedCohOnly)(i,j) = 0;
           }
        }


	TCanvas c("c", "c",100, 100);
	gStyle->SetOptStat(0);
	M_cohOnly.Draw("colz");
	c.Update();
	c.SaveAs("Original_Genie_Matrix_COHonly.pdf", "pdf");
	c.Clear();
	c.cd();
	M_updatedCohOnly->Draw("colz");
	c.Update();
	c.SaveAs("Updated_Genie_Matrix_COHonly.pdf", "pdf");
*/

/*	//update the genie uncertainty of COH pi0's to be 1.
	TMatrixT<double> temp_matrix = *m_genie_fractional_covariance_matrix;
	SBNspec temp_cv = *m_cv_spectrum;
	temp_cv.Keep("NCPi0Coh", 1.0);
	temp_cv.CalcFullVector();
	std::vector<double> temp_full = temp_cv.full_vector;

	for(int i=0; i<temp_full.size(); i++){
	   for(int j=0;j<temp_full.size(); j++){
		if(temp_full[i]!=0 && temp_full[j]!=0) (*m_genie_fractional_covariance_matrix)(i,j)=1.0;
	   }
	}
	std::cout << "check " << __LINE__ << std::endl;

	TMatrixT<double>* m_fluxxs = (TMatrixT<double>*)f_syst->Get("frac_covariance");
	this->RemoveNan(m_fluxxs);
	*m_fluxxs +=  *m_genie_fractional_covariance_matrix - temp_matrix;
	f_syst->cd();
	m_fluxxs->Write("updated_frac_covariance", TObject::kWriteDelete);
	TDirectory *individualDir = f_syst->GetDirectory("individualDir");
	individualDir->cd();
	m_genie_fractional_covariance_matrix->Write("updated_allUBGenie_frac_covariance", TObject::kWriteDelete);
	std::cout << "check " << __LINE__ << std::endl;

	for(int i=0; i<num_bins_total; i++){
	    for(int j=0; j< num_bins_total; j++){
		if(temp_matrix(i,j) != (*m_genie_fractional_covariance_matrix)(i,j)){
		    if(temp_full[i]==0 || temp_full[j] ==0)
			std::cout << "We got trouble" << std::endl;
		}
	    }
	}


	 TFile* f_all = new TFile("/uboone/app/users/gge/singlephoton/whipping_star/working_directory/SinglePhoton_test/FakeData/Poisson_10/updated_1g1p_1g0p_2g1p_2g0p_Combined_FluxXSDet_v11_Momentum_fracfixed.SBNcovar.root", "UPDATE");
         TMatrixT<double>* mall = (TMatrixT<double>*)f_all->Get("frac_covariance");
	this->RemoveNan(mall);
	TMatrixT<double> mall_temp  = *mall;
	mall_temp += *m_genie_fractional_covariance_matrix - temp_matrix;
	f_all->cd();
	mall_temp.Write("updated_frac_covariance",TObject::kWriteDelete);
	f_all->Close();
	std::cout << "check " << __LINE__ << std::endl;
*/


/*	//set genie matrix to be the same as total covariance matrix
	m_genie_fractional_covariance_matrix = new TMatrixT<double>(num_bins_total, num_bins_total);
	*m_genie_fractional_covariance_matrix = *m_full_fractional_covariance_matrix;
*/


    //get a sub of matrix
    /*	TMatrixT<double> temp_matrix = *m_genie_fractional_covariance_matrix;
        m_genie_fractional_covariance_matrix->ResizeTo(num_bins_total, num_bins_total);
     *m_genie_fractional_covariance_matrix  = temp_matrix.GetSub(90, 329, 90, 329, "S");
     */
    if(is_verbose) std::cout << "total bins " << num_bins_total << ", matrix size " << m_genie_fractional_covariance_matrix->GetNrows() << std::endl;


/*
	TMatrixT<double>* temp = (TMatrixT<double>*)f_syst->Get("frac_covariance");
	this->RemoveNan(temp);
	*m_full_fractional_covariance_matrix += *temp;	

	delete temp;
*/
	f_syst->Close();
	
	return 0;
}


int SBNsinglephoton::CalcFullButGenieFractionalCovarMatrix(){
    
    otag="SBNsinglephoton::CalcFullButGenieFractionalCovarMatrix\t|| ";

    if(m_full_fractional_covariance_matrix==nullptr || m_genie_fractional_covariance_matrix==nullptr){
        std::cout<< otag<< "Either full fractional covar matrix or genie fractional covariance matrix has NOT been setup yet."<< std::endl;
        exit(EXIT_FAILURE); 
    }

    if(is_verbose) std::cout << otag<< "as the name says.." << std::endl;
    m_full_but_genie_fractional_covariance_matrix = new TMatrixT<double>(num_bins_total, num_bins_total);
    *m_full_but_genie_fractional_covariance_matrix  = (*m_full_fractional_covariance_matrix) - (*m_genie_fractional_covariance_matrix);


/*
   // for agnoastic search only, zero out ANY correlation between "signal-like" and other subchannels
        SBNspec temp_cv = *m_cv_spectrum;
        temp_cv.Keep("NCDelta", 1.0);
        //temp_cv.Keep("NCDeltaLEE", 1.0);
        temp_cv.CalcFullVector();
        std::vector<double> temp_full = temp_cv.full_vector;
	temp_cv = *m_cv_spectrum;
	temp_cv.Scale("NCDelta", 0.0);
	//temp_cv.Scale("NCDeltaLEE", 0.0);
	std::vector<double> temp_others = temp_cv.full_vector;

        for(int i=0; i<temp_full.size(); i++){
           for(int j=0;j<temp_full.size(); j++){
                if(temp_full[i]*temp_full[j]==0 && temp_others[i]*temp_others[j]==0) (*m_full_but_genie_fractional_covariance_matrix)(i,j)=0.0;
           }
        }
*/

    return 0;
}

int SBNsinglephoton::AddCovarianceMatrix(const std::string &filename, const std::string &covar_name){

    if(is_verbose)
	std::cout << "SBNsinglephoton::AddCovarianceMatrix\t|| Add in another covariance matrix: " << covar_name << ", from file: " << filename << std::endl;

    TFile* ftemp = new TFile(filename.c_str(), "read");
    TMatrixT<double>* temp_matrix = (TMatrixT<double>*)ftemp->Get(covar_name.c_str());
    RemoveNan(temp_matrix);
    *m_full_fractional_covariance_matrix += *temp_matrix;

    ftemp->Close(); 
    return 0;

}

void SBNsinglephoton::ZeroOutCorrelation(TMatrixT<double> *pmatrix, const std::string &fname){
    otag="SBNsinglephoton::ZeroOutCorrelation\t|| ";
    if(is_verbose)
	std::cout << otag << "Zero out correlation between " << fname << " and other subchannels" << std::endl;

    // check size of the matrix
    if( (pmatrix->GetNrows() != pmatrix->GetNcols()) || (pmatrix->GetNrows() != num_bins_total )){
	std::cerr << otag << "Matrix size is wrong, Nrow: " << pmatrix->GetNrows() << ", Ncol: " << pmatrix->GetNcols() << ", num_bins_total: " << num_bins_total << std::endl;
 	exit(EXIT_FAILURE);
    }


    SBNspec temp_spec = *m_cv_spectrum;
    temp_spec.CalcFullVector();
    temp_spec.Scale(fname,0.0);
    std::vector<double> other_vec = temp_spec.full_vector;
    temp_spec = *m_cv_spectrum;
    temp_spec.Keep(fname, 1.0);
    std::vector<double> sig_vec = temp_spec.full_vector;

    for(int i=0; i<other_vec.size(); i++){
           for(int j=0;j<other_vec.size(); j++){
                if(other_vec[i]*other_vec[j]==0 && sig_vec[i]*sig_vec[j]==0) (*pmatrix)(i,j)=0.0;
           }
    }
  
}

void SBNsinglephoton::ZeroOutOffDiagonal(){
   int num_rows = m_full_fractional_covariance_matrix->GetNrows();
   for(int i = 0; i < num_rows; ++i)
       for( int j = 0; j < num_rows; ++j)
	   if(i != j ) (*m_full_fractional_covariance_matrix)(i,j) = 0.0;
}

void SBNsinglephoton::ZeroOutGenieCorrelation(const std::string &fname){

    if(m_full_fractional_covariance_matrix==nullptr || m_genie_fractional_covariance_matrix==nullptr){
        std::cout<< "SBNsinglephoton::ZeroOutGenieCorrelation\t|| Either full fractional covar matrix or genie fractional covariance matrix has NOT been setup yet."<< std::endl;
        exit(EXIT_FAILURE);
    }

    TMatrixT<double> temp_diff_matrix = (*m_full_fractional_covariance_matrix) - (*m_genie_fractional_covariance_matrix);

    //zero out correlation in genie matrix
    ZeroOutCorrelation(m_genie_fractional_covariance_matrix, fname);

    *m_full_fractional_covariance_matrix = temp_diff_matrix + (*m_genie_fractional_covariance_matrix);

}

void SBNsinglephoton::ZeroOutFullCorrelation(const std::string &fname){

    if(m_full_fractional_covariance_matrix==nullptr){
	std::cerr << "SBNsinglephoton::ZeroOutFullCorrelation\t|| Full covariance matrix hasn't been setup yet" << std::endl;
	exit(EXIT_FAILURE);
    }

    ZeroOutCorrelation(m_full_fractional_covariance_matrix, fname);
    
}

double SBNsinglephoton::ScaleFactor(double E, double factor, std::vector<double>& vec){
    otag = "SBNsinglephoton::ScaleFactor\t|| ";
    double scale = factor;
    switch(vec.size()){
        case 0:
            //std::cout<< otag << "No energy/momentum dependent scaling applied" << std::endl;
            break;
        case 1:
            //std::cout << otag << "Applying Linear energy/momentum dependent scaling!";
            scale += vec[0]*E;
            break;
        case 2:
            //std::cout << otag << "Applying 2nd order polynomial energy/momentum dependent scaling!" << std::endl;
            scale += vec[0]*E+vec[1]*E*E;
            break;
        default: scale += 0;
    }

    //if(scale < 0) scale = 1.0;  //don't allow negative weight
    //if(scale < 0) scale = 0;  //don't allow negative weight
    return scale;

}


int SBNsinglephoton::RemoveNan(TMatrixT<double>* M){

    int row_size = M->GetNrows();
    int col_size = M->GetNcols();

    if(is_verbose) std::cout<< "SBNsinglephoton::RemoveNan\t|| Remove the nan's from matrix " << M->GetName() << std::endl;
    for(int i=0; i< row_size; i++){
        for(int j=0; j< col_size; j++){
            if(std::isnan((*M)(i,j))) (*M)(i,j) =0.0;
        }
    }

    return 0;	
}


int SBNsinglephoton::GrabFitMap(){
	if(is_verbose) std::cout << "SBNsinglephoton::GrabFitMap\t||\tOpen file: "<< tag<<"_fit_output.root, and grab chi and BF index info"<< std::endl;
	fin = new TFile(Form("%s_fit_output.root", tag.c_str()), "UPDATE");
	TVectorD* v = (TVectorD*)fin->Get("bf_index");	
	int bf_index = (int)(*v)[0];
	std::vector<double>* p_vec_chi;
	fin->GetObject("vector_chi", p_vec_chi);

	m_map={{bf_index, *p_vec_chi}};		

	delete v;
	delete p_vec_chi;
	return 0;
}


int SBNsinglephoton::SaveHistogram(bool sensitivity){	
        return this->SaveHistogram(m_map, sensitivity);
}


int SBNsinglephoton::SaveHistogram(std::map<int, std::vector<double>>& inmap, bool sensitivity){
    otag="SBNsinglephoton::SaveHistogram\t|| ";

    if(inmap.empty()){
	std::cout << otag << "map is empty!!" << std::endl;
	exit(EXIT_FAILURE);	
    }else{
	std::cout << otag << std::endl;
        std::map<int, std::vector<double>>::const_iterator itmap = inmap.cbegin();
	int best_point = itmap->first;
	std::vector<double> vec_chi = itmap->second;
	double chi_min = *std::min_element(vec_chi.begin(), vec_chi.end());
	std::for_each(vec_chi.begin(), vec_chi.end(), [&chi_min](double& d){d -= chi_min;});  //get the delta_chi vector


	std::map<std::string, std::string> title_map={{"NCDelta", "Enhancement for NC #Delta->N#gamma BR"},
						      {"NCDeltaLEE", "Excess Scale Factor (x GENIE SM NC #Delta->N#gamma prediction)"},
						      {"NCPi0Coh", "Flat Scaling Factor for NC 1#pi^{0} Coherent"},
						      {"NCPi0NotCoh", "Flat Scaling Factor for NC 1#pi^{0} Non-Coherent"},
						      {"NCPi0NotCohMom", "Momentum-dependent Term Coefficient for NC 1#pi^{0} Non-Coherent"},
						      {"All", "Flat normalization factor for all"}};

        int m_poly_best_index = best_point/m_flat_total_gridpoints;
        int m_best_index = best_point%m_flat_total_gridpoints;

	  
	if(m_grid.f_num_dimensions == 2){
	
	   if(!m_bool_poly_grid){
	    	if(is_verbose) std::cout<< otag<< "Case: NCpi0 normalization fit, no energy/momentum dependent scaling!" << std::endl;
		auto grid_x = m_grid.f_dimensions.at(0);
		auto grid_y = m_grid.f_dimensions.at(1);

		// fill original 2D chi surface
		TH2D* hr = new TH2D("h_chi_surface", Form("#Delta#chi^{2}_{CNP} surface; %s;%s", title_map[grid_y.GetName()].c_str(), title_map[grid_x.GetName()].c_str()), grid_y.GetNPoints(), &grid_y.GetEdges()[0], grid_x.GetNPoints(), &grid_x.GetEdges()[0]);
		std::vector<double> mchi_grid_x(grid_x.GetNPoints(), DBL_MAX);
		std::vector<double> mchi_grid_y(grid_y.GetNPoints(), DBL_MAX);
		for(size_t i = 0; i != vec_chi.size() ; ++i){
		     size_t grid_x_index = i/grid_y.GetNPoints();
		     size_t grid_y_index = i%grid_y.GetNPoints();
		     mchi_grid_x[grid_x_index] = std::min(mchi_grid_x[grid_x_index], vec_chi[i]);
		     mchi_grid_y[grid_y_index] = std::min(mchi_grid_y[grid_y_index], vec_chi[i]);
		     hr->SetBinContent(grid_y_index+1, grid_x_index+1, vec_chi[i]);
		}
		hr->Write();

		//generate interpolated 2D chi surface
		TH2D h_chi_surface = this->Do2DInterpolation(m_interpolation_number, grid_y.GetPoints(), grid_x.GetPoints(),vec_chi, tag);
		h_chi_surface.SetName("h_chi_interpolated_surface");
		h_chi_surface.SetTitle(Form("#Delta#chi^{2}_{CNP} surface; %s;%s", title_map[grid_y.GetName()].c_str(), title_map[grid_x.GetName()].c_str()));
		h_chi_surface.Write();

		hr->GetXaxis()->SetRangeUser(h_chi_surface.GetXaxis()->GetXmin(), h_chi_surface.GetXaxis()->GetXmax());
                hr->GetYaxis()->SetRangeUser(h_chi_surface.GetYaxis()->GetXmin(), h_chi_surface.GetYaxis()->GetXmax());
		//std::cout << "Interpolated x range: " << h_chi_surface.GetXaxis()->GetXmin() << ", " << h_chi_surface.GetXaxis()->GetXmax() << std::endl;
		//std::cout << "Trimmed x range:" << hr->GetXaxis()->GetXmin() << ", " << hr->GetXaxis()->GetXmax()  << std::endl;
		//draw contours
		std::vector<TGraph> contour_graph = this->FindContour(h_chi_surface, 3, tag);  //want 3 contour
	        //TCanvas* c_canvas=(TCanvas*)this->DrawContour(hr, contour_graph, m_vec_grid[m_best_index]).Clone();	
	        this->DrawContour(hr, contour_graph, m_vec_grid[m_best_index]);	

		std::vector<TH1D*> vec_mchi_hist;
		vec_mchi_hist.push_back(new TH1D(Form("h_mchi_%s", (grid_x.GetName()).c_str()), Form("Marginalized #Delta#chi^{2}_{CNP} distribution; %s; #Delta#chi^{2}_{CNP}", title_map[grid_x.GetName()].c_str()), grid_x.GetNPoints(), &grid_x.GetEdges()[0]));
	 	vec_mchi_hist.push_back(new TH1D(Form("h_mchi_%s", (grid_y.GetName()).c_str()), Form("Marginalized #Delta#chi^{2}_{CNP} distribution; %s; #Delta#chi^{2}_{CNP}", title_map[grid_y.GetName()].c_str()), grid_y.GetNPoints(), &grid_y.GetEdges()[0]));

		for(size_t i = 0; i != grid_x.GetNPoints(); ++i)
		    vec_mchi_hist[0]->SetBinContent(i+1, mchi_grid_x[i]);

		for(size_t i =0; i!= grid_y.GetNPoints(); ++i)
		    vec_mchi_hist[1]->SetBinContent(i+1, mchi_grid_y[i]);

		DrawMarginalizedChi(vec_mchi_hist, tag, sensitivity);	
	   }else{
		if(is_verbose) std::cout<< otag<<"Case: NCpi0 normalization fit, with energy/momentum dependent scaling to " << m_poly_grid.f_num_dimensions << "nd order!" << std::endl;
		 auto grid_1order = m_poly_grid.f_dimensions.at(0);

		 //grab the flat grid associated with poly_grid
		 NGridDimension grid_flat = grid_1order;
		 NGridDimension grid_fother = grid_1order;
		 int period = 0;  // very important variable!! 
		 for(auto &grid:m_grid.f_dimensions){
			period *= grid.GetNPoints();
			if(grid.GetName() == grid_1order.GetName()){
			    grid_flat = grid;
			    period = 1;
			}else{
			    grid_fother = grid;
			}
		 }		

		 //2D marginalized chi^2 surface
		 TH2D h_mchi_poly("h_mchi_poly", Form("marginalized #Delta#chi^{2}_{CNP}; %s; %s", title_map[grid_flat.GetName()].c_str(), title_map["NCPi0NotCohMom"].c_str()), grid_flat.GetNPoints(), &grid_flat.GetEdges()[0], grid_1order.GetNPoints(), &grid_1order.GetEdges()[0]);
		 TH2D h_mchi_flat("h_mchi_flat", Form("marginalized #Delta#chi^{2}_{CNP}; %s; %s", title_map[grid_fother.GetName()].c_str(), title_map[grid_flat.GetName()].c_str()), grid_fother.GetNPoints(), &grid_fother.GetEdges()[0], grid_flat.GetNPoints(), &grid_flat.GetEdges()[0]);
		 TH2D h_mchi_mix("h_mchi_mix", Form("marginalized #Delta#chi^{2}_{CNP}; %s; %s", title_map[grid_fother.GetName()].c_str(), title_map["NCPi0NotCohMom"].c_str()), grid_fother.GetNPoints(), &grid_fother.GetEdges()[0], grid_1order.GetNPoints(), &grid_1order.GetEdges()[0]);

		 // 1D marginalzied chi^2 curve
		 std::vector<TH1D*> vh_mchi(3);
		 vh_mchi[0] = new TH1D(Form("h_mchi_flat_%s", (grid_flat.GetName()).c_str()), Form("marginalized #Delta#chi^{2}_{CNP}; %s;marginalized #Delta#chi^{2}_{CNP}",title_map[grid_flat.GetName()].c_str()), grid_flat.GetNPoints(), &grid_flat.GetEdges()[0]);
		 vh_mchi[1] = new TH1D(Form("h_mchi_poly_%s", (grid_1order.GetName()).c_str()), Form("marginalized #Delta#chi^{2}_{CNP}; %s;marginalized #Delta#chi^{2}_{CNP}",title_map["NCPi0NotCohMom"].c_str()), grid_1order.GetNPoints(), &grid_1order.GetEdges()[0]);
		 vh_mchi[2] = new TH1D(Form("h_mchi_%s", (grid_fother.GetName()).c_str()), Form("marginalized #Delta#chi^{2}_{CNP}; %s; marginalized #Delta#chi^{2}_{CNP}",title_map[grid_fother.GetName()].c_str()),grid_fother.GetNPoints(), &grid_fother.GetEdges()[0]);

		 // initialize histogram content
		 for(int i=0;i<grid_flat.GetNPoints(); i++){
			vh_mchi[0]->SetBinContent(i+1, DBL_MAX);
			for(int j=0;j<grid_1order.GetNPoints(); j++)
			    h_mchi_poly.SetBinContent(i+1, j+1, DBL_MAX);
			for(int j=0;j<grid_fother.GetNPoints(); j++)
			    h_mchi_flat.SetBinContent(j+1, i+1, DBL_MAX);
		 } 	
		 for(int i=0;i<grid_fother.GetNPoints(); i++){
			vh_mchi[2]->SetBinContent(i+1,DBL_MAX);
			for(int j=0;j<grid_1order.GetNPoints(); j++)
			    h_mchi_mix.SetBinContent(i+1, j+1, DBL_MAX);
		 } 
		for(int i=0; i<grid_1order.GetNPoints(); i++) vh_mchi[1]->SetBinContent(i+1, DBL_MAX);
	
		 //start filling chi^2 surface/curve histograms
		 for(int i=0; i<vec_chi.size(); i++){
		     int temp_poly_index = i/m_flat_total_gridpoints;
        	     int temp_mgrid_index = i%m_flat_total_gridpoints; 
		     int flat_index = (temp_mgrid_index%(period* grid_flat.GetNPoints()))/period; //get its index in flat dimension for current point
			// needs more work
		     int fother_index = temp_mgrid_index%grid_fother.GetNPoints(); //index for another flat dimension
		     if(period == 1) fother_index = temp_mgrid_index/grid_flat.GetNPoints();
		     if(vec_chi[i] < h_mchi_poly.GetBinContent(flat_index+1, temp_poly_index+1)) h_mchi_poly.SetBinContent(flat_index+1, temp_poly_index+1, vec_chi[i]);
		     if(vec_chi[i] < h_mchi_flat.GetBinContent(fother_index+1, flat_index+1)) h_mchi_flat.SetBinContent(fother_index+1, flat_index+1, vec_chi[i]);
		     if(vec_chi[i] < h_mchi_mix.GetBinContent(fother_index+1, temp_poly_index+1)) h_mchi_mix.SetBinContent(fother_index+1, temp_poly_index+1, vec_chi[i]);
		     if(vec_chi[i] < vh_mchi[0]->GetBinContent(flat_index+1)) vh_mchi[0]->SetBinContent(flat_index+1, vec_chi[i]);
		     if(vec_chi[i] < vh_mchi[1]->GetBinContent(temp_poly_index+1)) vh_mchi[1]->SetBinContent(temp_poly_index+1, vec_chi[i]);
		     if(vec_chi[i] < vh_mchi[2]->GetBinContent(fother_index+1)) vh_mchi[2]->SetBinContent(fother_index+1, vec_chi[i]);
		 }	 
		fin->cd(); 
		h_mchi_poly.Write();
		h_mchi_flat.Write();
		h_mchi_mix.Write();

		DrawMarginalizedChi(vh_mchi, tag, sensitivity);	


		//save marginalized chi
		std::vector<double> marginalized_chi_poly, marginalized_chi_mix;
		for(int i=0; i< grid_1order.GetNPoints();i++){
		    for(int j=0; j< grid_flat.GetNPoints(); j++)
			marginalized_chi_poly.push_back(h_mchi_poly.GetBinContent(j+1, i+1));
		    for(int j=0; j< grid_fother.GetNPoints(); j++)
			marginalized_chi_mix.push_back(h_mchi_mix.GetBinContent(j+1, i+1));
		}


		//save marginalized chi
		std::vector<double> marginalized_chi_flat;
		for(int i=0; i< grid_flat.GetNPoints();i++)
		    for(int j=0; j< grid_fother.GetNPoints(); j++)
			marginalized_chi_flat.push_back(h_mchi_flat.GetBinContent(j+1, i+1));

		//draw the contour
		TH2D h_mchi_flatinter = this->Do2DInterpolation(m_interpolation_number, grid_fother.GetPoints(), grid_flat.GetPoints(), marginalized_chi_flat, tag+"_Flat");
		std::vector<TGraph> vg_mchi_flat_contour = this->FindContour(h_mchi_flatinter, 2, tag+"_Flat");
		h_mchi_flat.GetXaxis()->SetRangeUser(h_mchi_flatinter.GetXaxis()->GetXmin(), h_mchi_flatinter.GetXaxis()->GetXmax());
                h_mchi_flat.GetYaxis()->SetRangeUser(h_mchi_flatinter.GetYaxis()->GetXmin(), h_mchi_flatinter.GetYaxis()->GetXmax());
		this->DrawContour(&h_mchi_flat, vg_mchi_flat_contour, tag+"_Flat",std::vector<double>{});

		//draw the contour
		TH2D h_mchi_polyinter = this->Do2DInterpolation(m_interpolation_number, grid_flat.GetPoints(), grid_1order.GetPoints(), marginalized_chi_poly, tag+"_Poly");
		std::vector<TGraph> vg_mchi_contour = this->FindContour(h_mchi_polyinter, 3, tag+"_Poly");
		h_mchi_poly.GetXaxis()->SetRangeUser(h_mchi_polyinter.GetXaxis()->GetXmin(), h_mchi_polyinter.GetXaxis()->GetXmax());
                h_mchi_poly.GetYaxis()->SetRangeUser(h_mchi_polyinter.GetYaxis()->GetXmin(), h_mchi_polyinter.GetYaxis()->GetXmax());
		this->DrawContour(&h_mchi_poly, vg_mchi_contour,tag+"_Poly", std::vector<double>{});

		TH2D h_mchi_mixinter = this->Do2DInterpolation(m_interpolation_number, grid_fother.GetPoints(), grid_1order.GetPoints(), marginalized_chi_mix, tag+"_Mix");
		std::vector<TGraph> vg_mchi_mix_contour = this->FindContour(h_mchi_mixinter, 3, tag+"_Mix");
		h_mchi_mix.GetXaxis()->SetRangeUser(h_mchi_mixinter.GetXaxis()->GetXmin(), h_mchi_mixinter.GetXaxis()->GetXmax());
                h_mchi_mix.GetYaxis()->SetRangeUser(h_mchi_mixinter.GetYaxis()->GetXmin(), h_mchi_mixinter.GetYaxis()->GetXmax());
		this->DrawContour(&h_mchi_mix, vg_mchi_mix_contour,tag+"_Mix", std::vector<double>{});


	   } //end of if_poly_grid loop
	} //end of 2 dimension case
	else if(m_grid.f_num_dimensions == 1){
	   TH1D* h_dchi=nullptr;
	   NGridDimension xgrid = m_grid.f_dimensions.at(0);
	   if(xgrid.GetName() == "NCDelta") h_dchi = new TH1D("h_delta_chi", Form("#Delta#chi^{2}_{CNP} distribution;%s; #Delta#chi^{2}_{CNP} ",title_map[xgrid.GetName()].c_str()), m_grid.f_num_total_points, &xgrid.GetEdges()[0]);
	   else if(xgrid.GetName() == "NCDeltaLEE" ) h_dchi = new TH1D("h_delta_chi", Form("#Delta#chi^{2}_{CNP} distribution;%s; #Delta#chi^{2}_{CNP} ",title_map[xgrid.GetName()].c_str()), m_grid.f_num_total_points, &xgrid.GetEdges()[0]);
	   //else if(xgrid.GetName() == "NCDeltaLEE" ) h_dchi = new TH1D("h_delta_chi", Form("#Delta#chi^{2}_{CNP} distribution;%s; #Delta#chi^{2}_{CNP} ",title_map[xgrid.GetName()].c_str()), m_grid.f_num_total_points, (xgrid.GetMin())*2, (xgrid.GetMax())*2);
	   for(int i=0 ;i< vec_chi.size(); i++){
		std::vector<double> ipoint = m_vec_grid[i];
		//h_dchi->Fill(ipoint[0], vec_chi[i]);
		h_dchi->SetBinContent(i+1, vec_chi[i]);
	   }

	   DrawMarginalizedChi(std::vector<TH1D*>{h_dchi}, tag,sensitivity, std::vector<double>{m_vec_grid[m_best_index][0], chi_min, vec_chi.at(m_cv_spec_index)});
	}//end of 1 dimension case
	else{

	   if(!m_bool_poly_grid){
		if(is_verbose) std::cout<< otag << "Case: NC delta and pi0 combined fit, no energy/momentum dependent scaling!" << std::endl;
		auto grid_x = m_grid.f_dimensions.at(0);
		auto grid_y = m_grid.f_dimensions.at(1);
		auto grid_z = m_grid.f_dimensions.at(2);   //assume grid_Z is the grid for NCDeltaLEE here.
		std::vector<double> temp_best_point = m_vec_grid[m_best_index];

		//marginalize over 1 parameter	
		//TH2D* h_mchi2_xy = new TH2D("h_mchi2_xy", Form("marginalized #Delta#chi^{2}_{CNP} surface; %s;%s", title_map[grid_x.GetName()].c_str(), title_map[grid_y.GetName()].c_str()), grid_x.GetNPoints(), grid_x.GetMin(), grid_x.GetMax(), grid_y.GetNPoints(), grid_y.GetMin(), grid_y.GetMax());
		TH2D* h_mchi2_xy = new TH2D("h_mchi2_xy", Form("marginalized #Delta#chi^{2}_{CNP} surface; %s;%s", title_map[grid_y.GetName()].c_str(), title_map[grid_x.GetName()].c_str()), grid_y.GetNPoints(), &grid_y.GetEdges()[0], grid_x.GetNPoints(), &grid_x.GetEdges()[0]);
		TH2D* h_mchi2_yz = new TH2D("h_mchi2_yz", Form("marginalized #Delta#chi^{2}_{CNP} surface; %s;%s", title_map[grid_y.GetName()].c_str(), title_map[grid_z.GetName()].c_str()), grid_y.GetNPoints(), &grid_y.GetEdges()[0], grid_z.GetNPoints(), (grid_z.GetMin())*2+1, (grid_z.GetMax())*2+1);
		TH2D* h_mchi2_xz = new TH2D("h_mchi2_xz", Form("marginalized #Delta#chi^{2}_{CNP} surface; %s;%s", title_map[grid_x.GetName()].c_str(), title_map[grid_z.GetName()].c_str()), grid_x.GetNPoints(), &grid_x.GetEdges()[0], grid_z.GetNPoints(), (grid_z.GetMin())*2+1, (grid_z.GetMax())*2+1);
		//global minimum
		TH2D* h_gchi2_xy = new TH2D("h_gchi2_xy", Form("h_gchi2_xy; %s;%s", title_map[grid_x.GetName()].c_str(), title_map[grid_y.GetName()].c_str()), grid_x.GetNPoints(), &grid_x.GetEdges()[0], grid_y.GetNPoints(), &grid_y.GetEdges()[0]);
		TH2D* h_gchi2_yz = new TH2D("h_gchi2_yz", Form("h_gchi2_yz; %s;%s", title_map[grid_y.GetName()].c_str(), title_map[grid_z.GetName()].c_str()), grid_y.GetNPoints(), &grid_y.GetEdges()[0], grid_z.GetNPoints(), (grid_z.GetMin())*2+1, (grid_z.GetMax())*2+1);
		TH2D* h_gchi2_xz = new TH2D("h_gchi2_xz", Form("h_gchi2_xz; %s;%s", title_map[grid_x.GetName()].c_str(), title_map[grid_z.GetName()].c_str()), grid_x.GetNPoints(), &grid_x.GetEdges()[0], grid_z.GetNPoints(), (grid_z.GetMin())*2+1, (grid_z.GetMax())*2+1);

		//minimize over two parameters
		TH1D* h_chi_delta = new TH1D("h_chi_delta", Form("h_chi_delta; %s;#Delta#chi^{2}_{CNP}",title_map[grid_z.GetName()].c_str()), grid_z.GetNPoints(), (grid_z.GetMin())*2+1, (grid_z.GetMax())*2+1);

		//vectors that stores marginalized chi
		std::vector<double> mchi_xy, mchi_yz, mchi_xz, gchi_xy, gchi_yz, gchi_xz;

		for(int ix=1;ix <= grid_x.GetNPoints(); ix++){
		        //for(int iy=1; iy <= grid_y.GetNPoints(); iy++) h_mchi2_xy->SetBinContent(ix, iy, DBL_MAX);
		        for(int iy=1; iy <= grid_y.GetNPoints(); iy++) h_mchi2_xy->SetBinContent(iy, ix, DBL_MAX);
        		for(int iz=1; iz <= grid_z.GetNPoints(); iz++) h_mchi2_xz->SetBinContent(ix, iz, DBL_MAX);
   		}

   		for(int iz=1; iz <= grid_z.GetNPoints(); iz++){
        		for(int iy=1; iy <= grid_y.GetNPoints(); iy++)  h_mchi2_yz->SetBinContent(iy, iz, DBL_MAX);
        		h_chi_delta->SetBinContent(iz, DBL_MAX);
   		}


		for(int ix=0; ix < grid_x.GetNPoints(); ix++){
       		    for(int iy=0; iy < grid_y.GetNPoints(); iy++){
           		for(int iz=0 ; iz< grid_z.GetNPoints(); iz++){
                	    int ip = ix*grid_y.GetNPoints()*grid_z.GetNPoints() + iy*grid_z.GetNPoints() + iz; // index of grid point
		            std::vector<double> point = m_vec_grid[ip];


                            //marginalized minimum
                            //conditional operator, saver the smaller chi.
                	 /*   if(vec_chi[ip]< h_mchi2_xy->GetBinContent(ix+1, iy+1)){
                        	 h_mchi2_xy->SetBinContent(ix+1, iy+1, vec_chi[ip]);
                         	//std::cout << "chi2 value: " << chi[ip] << std::endl;

		 	    }
                	   */ if(vec_chi[ip]< h_mchi2_xy->GetBinContent(iy+1, ix+1)){
                        	 h_mchi2_xy->SetBinContent(iy+1, ix+1, vec_chi[ip]);
                         	//std::cout << "chi2 value: " << chi[ip] << std::endl;

		 	    }
                	    if(vec_chi[ip]< h_mchi2_xz->GetBinContent(ix+1, iz+1)) h_mchi2_xz->SetBinContent(ix+1, iz+1, vec_chi[ip]);
                	    if(vec_chi[ip]< h_mchi2_yz->GetBinContent(iy+1, iz+1)) h_mchi2_yz->SetBinContent(iy+1, iz+1, vec_chi[ip]);


               		    //global minimum
                	    if(point[2] == temp_best_point[2]){
				 h_gchi2_xy->Fill(point[0], point[1], vec_chi[ip]);
				 gchi_xy.push_back(vec_chi[ip]);}
                	    if(point[1] == temp_best_point[1]){
				 h_gchi2_xz->Fill(point[0], point[2]*2+1, vec_chi[ip]);
				 gchi_xz.push_back(vec_chi[ip]);}
               		    if(point[0] == temp_best_point[0]){
				 h_gchi2_yz->Fill(point[1], point[2]*2+1, vec_chi[ip]);
				 gchi_yz.push_back(vec_chi[ip]);}

                 	    //marginalize two parameters
                	    if(vec_chi[ip] < h_chi_delta->GetBinContent(iz+1)) h_chi_delta->SetBinContent(iz+1, vec_chi[ip]);
           		}
        	    }
   		}

	/*	for(int iy=1;iy <= grid_y.GetNPoints(); iy++){
		        for(int ix=1; ix <= grid_x.GetNPoints(); ix++) mchi_xy.push_back(h_mchi2_xy->GetBinContent(ix, iy));
   		}
	*/	for(int ix=1;ix <= grid_x.GetNPoints(); ix++){
		        for(int iy=1; iy <= grid_y.GetNPoints(); iy++) mchi_xy.push_back(h_mchi2_xy->GetBinContent(iy, ix));
   		}

   		for(int iz=1; iz <= grid_z.GetNPoints(); iz++){
        		for(int iy=1; iy <= grid_y.GetNPoints(); iy++)  mchi_yz.push_back(h_mchi2_yz->GetBinContent(iy, iz));
        		for(int ix=1; ix <= grid_x.GetNPoints(); ix++)  mchi_xz.push_back(h_mchi2_xz->GetBinContent(ix, iz));
   		}
		
		h_mchi2_xy->Write(); h_mchi2_xz->Write(); h_mchi2_yz->Write();
		h_gchi2_xy->Write(); h_gchi2_xz->Write(); h_gchi2_yz->Write();
	
		DrawMarginalizedChi(std::vector<TH1D*>{h_chi_delta}, tag, sensitivity);


		if(grid_x.GetNPoints() >=5 && grid_y.GetNPoints()>=5){
			//draw 1,2,3 sigma contours for h_mchi2_zy
			TH2D h_mchi2_xy_inter =this->Do2DInterpolation(m_interpolation_number,grid_y.GetPoints(), grid_x.GetPoints(), mchi_xy, tag+"_XY"); 
			std::vector<TGraph> g_mchi_contour = this->FindContour(h_mchi2_xy_inter, 3, tag+"_XY");

			h_mchi2_xy->GetXaxis()->SetRangeUser(h_mchi2_xy_inter.GetXaxis()->GetXmin(), h_mchi2_xy_inter.GetXaxis()->GetXmax());
			h_mchi2_xy->GetYaxis()->SetRangeUser(h_mchi2_xy_inter.GetYaxis()->GetXmin(), h_mchi2_xy_inter.GetYaxis()->GetXmax());
			std::cout << "x max: " << h_mchi2_xy->GetXaxis()->GetXmax() << " "<< h_mchi2_xy_inter.GetXaxis()->GetXmax() << std::endl;
			this->DrawContour(h_mchi2_xy, g_mchi_contour, std::vector<double>{});
		}
	   } //end of the m_bool_poly_grid=false loop
	}//end of 3 dimension case
    fin->Close();
    } //end of map check

    return 0;
}



//2D interpolation
//x, y vector should only be in increasing order
//NOTE::elements of the vector should first change with vector 'x', then change with vector 'y'!! Order is important here!
TH2D SBNsinglephoton::Do2DInterpolation(int inter_number, const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& value, std::string intag){
	double x_min = x[0];
	double x_max = x[x.size()-1];
	double y_min = y[0];
	double y_max = y[y.size() -1];
	//interpolation step size
	double x_step = fabs(x_max - x_min)/double(inter_number -1);
	double y_step = fabs(y_max - y_min)/double(inter_number -1);

	std::cout << "SBNsinglephoton::Do2DInterpolation\t||\tCreate interpolated 2D plot with tag: "<< intag << ", interpolation number " << inter_number<< std::endl;
	//delete gROOT->FindObject("h_inter");
	TH2D h_inter(("h_inter_"+intag).c_str(), ("h_inter_"+intag).c_str(), inter_number, x_min, x_max, inter_number, y_min, y_max);
	
	const gsl_interp2d_type *T = gsl_interp2d_bicubic;  //bicubic interpolation
	gsl_spline2d *spline = gsl_spline2d_alloc(T, x.size(), y.size());  
	gsl_interp_accel *xacc = gsl_interp_accel_alloc();
        gsl_interp_accel *yacc = gsl_interp_accel_alloc();

        /* initialize interpolation */
	gsl_spline2d_init(spline,  &x[0], &y[0], &value[0] , x.size(), y.size());
	for (int i = 0; i < inter_number; i++){

              double xi;

              if(i == (inter_number -1)) xi = x_min + (i-0.5)*x_step;  //to fill the last bin
              else xi = x_min + i*x_step;
	    
              for (int j = 0; j < inter_number; j++){

                  double yj;
                  if(j == (inter_number -1)) yj = y_min + (j-0.5)*y_step;
                  else yj = y_min + j*y_step;

		  double zij = gsl_spline2d_eval(spline, xi, yj, xacc, yacc);  
		  h_inter.Fill(xi, yj, zij);
	      }
	}

	//freee up the pointers
	gsl_spline2d_free(spline);
        gsl_interp_accel_free(xacc);
        gsl_interp_accel_free(yacc);
	

	return h_inter;	
}


//return 1/2/3/4/5 sigma contours as vector of TGraph
std::vector<TGraph> SBNsinglephoton::FindContour(TH2D &hin, int n, std::string intag){
   if(is_verbose) std::cout<< "SBNsinglephoton::FindContour\t||\tTrying to find " << n << "contours" << std::endl;
   
   std::vector<double> full_contour{2.30, 6.18,11.83, 19.35, 28.23}; // chi2 value for 2 dof with 1, 2, 3, 4, 5 sigma confidence level
   std::vector<std::string> full_ctour_string{"1sigma","2sigma", "3sigma", "4sigma", "5sigma"};
   std::vector<std::string> full_name_string{"1#sigma", "2#sigma", "3#sigma", "4#sigma", "5#sigma"};
   std::vector<int> full_color{kGreen-7, kCyan-7, kMagenta-7, kRed-7, kBlue-7};//color of the contour: kBlue, kMagenta, kRed
   std::vector<int> full_style{kSolid, kDashed, kDotted, kDashDotted, 5};
   std::vector<int> full_marker_style{24,25,26,27,28};

   std::vector<double> contour(full_contour.begin(), full_contour.begin()+n); // chi2 value for 2 dof with 1 sigma, 90%, 2 sigma and 3 sigma confidence level
   std::vector<std::string> CL_string(full_ctour_string.begin(), full_ctour_string.begin()+n);
   std::vector<std::string> Name_string(full_name_string.begin(), full_name_string.begin()+n);
   std::vector<int> color_vec(full_color.begin(), full_color.begin()+n); 
   std::vector<int> style_vec(full_style.begin(), full_style.begin()+n);
   std::vector<int> mstyle_vec(full_marker_style.begin(), full_marker_style.begin()+n); 

   std::vector<TGraph> out_graph(n);
   //draw contour
   TCanvas c((intag+"c").c_str(), (intag+"c").c_str());
   c.cd();
   hin.SetContour((int)contour.size(), &contour[0]);
   hin.Draw("CONT Z LIST");//"LIST" generates a list of TGraph for each contour
   //grab contour object
   c.Update();
   gROOT->GetListOfSpecials()->Print();
   TObjArray *conts = (TObjArray*)gROOT->GetListOfSpecials()->FindObject("contours");
   
   //loop over contours
   for(int i=0; i<conts->GetSize(); i++){
	//create a new TList pointer inside a loop, it will get deleted everytime during the loop --> no memory leak here.
        TList* con_list = (TList*)conts->At(i);   //a list of TGraph for i'th contour.
        double x_min = DBL_MAX;
        double x_max = DBL_MIN;
        double y_min = DBL_MAX;
        double y_max = DBL_MIN;

        //out_graph[i] = new TGraph();   // one TGraph for one contour
        out_graph[i].SetName(Form("graph_%s", CL_string[i].c_str()));
        out_graph[i].SetTitle(Form("%s", Name_string[i].c_str()));
        out_graph[i].SetLineStyle(style_vec[i]);
        out_graph[i].SetLineWidth(2);
        out_graph[i].SetMarkerColor(color_vec[i]);
        out_graph[i].SetMarkerStyle(mstyle_vec[i]);
	out_graph[i].SetMarkerSize(0.1);
	TGraph * c_graph = (TGraph*)con_list->First();  //grab the TGraph

        for(int j=0; j< con_list->GetSize() ; j++){
            TGraph * temp_graph= (TGraph*)c_graph->Clone();
            double x,y;
            for(int k =0; k< temp_graph->GetN(); k++){
                temp_graph->GetPoint(k, x, y);
                if(x < x_min) x_min = x;
                if(x > x_max) x_max = x;

                if(y < y_min) y_min =y;
                if(y > y_max) y_max = y;

                out_graph[i].SetPoint(out_graph[i].GetN(), x, y);
            }
            c_graph=(TGraph*)con_list->After(c_graph);
        }

        std::cout << "Contour " << CL_string[i] << ": " << contour[i] << std::endl;
        std::cout << "(2D projected) range for x : " << std::setprecision(3) << x_min << "~" << x_max << std::endl;
        std::cout << "(2D projected) range for y : " << std::setprecision(3) << y_min << "~" << y_max << std::endl;
    }

   //conts->SetOwner(kTRUE);
   //conts->Clear();
   return out_graph;
}



void SBNsinglephoton::DrawContour(TH2D* h, std::vector<TGraph>& v_graph, std::vector<double> bf_point){
	this->DrawContour(h, v_graph, tag, bf_point);
	return;
}

void SBNsinglephoton::DrawContour(TH2D* h, std::vector<TGraph>& v_graph, std::string intag,  std::vector<double> bf_point){

   if(is_verbose) std::cout << "SBNsinglephoton::DrawContour||\tCreating contour pdf: "<<intag<<"_pretty_chi_contours.pdf"<< std::endl;
   //delete gROOT->FindObject((intag+"_pretty_contour").c_str());

   TCanvas *c_canvas = new TCanvas((intag+"_pretty_contour").c_str(), (intag+"_pretty_contour").c_str());
   gStyle->SetOptStat(0);
   c_canvas->cd();
   TLegend *leg = new TLegend(0.68, 0.7, 0.9,0.9);
   leg->SetFillStyle(0); 
   leg->SetBorderSize(0);

   h->Draw("colz");
   TMarker* marker=nullptr;
  
   if(bf_point.size() == 2){
	marker=new TMarker(bf_point.at(1), bf_point.at(0), 29);
	marker->SetMarkerColor(kWhite);
	marker->SetMarkerSize(3);
	marker->Draw();
	leg->AddEntry(marker, "Best-fit Point", "P");
   }

   for(auto& g:v_graph){
        g.SetLineColor(kWhite);
        g.SetMarkerColor(kWhite);
        g.Draw("same l");
        leg->AddEntry(&g, g.GetTitle(), "LP");
  }
  leg->SetTextSize(0.043);
  leg->Draw();
  c_canvas->SaveAs((intag+"_pretty_chi_contours.pdf").c_str(),"pdf");

  //delete marker;
  return;
}


int SBNsinglephoton::PrintOutFitInfo(const std::map<int, std::vector<double>>& inmap, std::string intag, bool print_all){
    std::map<int, std::vector<double>>::const_iterator itmap = inmap.begin();

    int best_point = itmap->first;
    std::vector<double> vec_chi = itmap->second;
    double chi_min = *std::min_element(vec_chi.begin(), vec_chi.end());
    int ndof = num_bins_total_compressed - m_fit_dimension;  //degree of freedom of the chi minimum
    if(intag.find("CV") != std::string::npos){
    	std::cout << intag<<": Chi2 minimum/ndof is " << chi_min << "/"<< ndof << ", Pvalue: " << TMath::Prob(chi_min, ndof) << std::endl;
    }else{
    	std::cout << intag<<": BestFit Chi2 minimum/ndof is " << chi_min << "/"<< ndof << ", Pvalue: " << TMath::Prob(chi_min, ndof) << std::endl;
    }

    //----- Print out BF information -----------------
    //best-fit point index at two grid
    int m_poly_grid_index = best_point/m_flat_total_gridpoints;
    int m_grid_index = best_point%m_flat_total_gridpoints;

    //print out info
    std::vector<double> m_point = m_vec_grid[m_grid_index];
    std::cout << intag<<": At flat point";
    for(int i=0; i<m_point.size(); i++)
        std::cout << ", " << m_grid.f_dimensions[i].GetName() << ": "<< m_point[i] ;
    if(m_bool_poly_grid){
        std::vector<double> m_poly_point = m_vec_poly_grid[m_poly_grid_index];

        std::cout << ", and energy/momentum dependent shift for NCpi0 non-coherent component with";
        for(int i=0 ; i< m_poly_point.size();i++)
            std::cout << ", "<< i+1 <<"nd order factor: "<< m_poly_point[i];
    }
    std::cout << std::endl;

/*    //------- Print out CV chi^2 -------------------
    if(intag.find("CV") != std::string::npos){
	std::cout << intag << ": GENIE CV Chi2/ndof evaluated at CentralValue: " << vec_chi[m_cv_spec_index]<< "/" << num_bins_total_compressed << ", Pvalue: " << TMath::Prob(vec_chi[m_cv_spec_index], num_bins_total_compressed) << std::endl;
    }else{
	std::cout << intag << ": GENIE CV Chi2 evaluated at BestFit: " << vec_chi[m_cv_spec_index] << ", Dchi: " << vec_chi[m_cv_spec_index] - chi_min << ", Fitting variable dimensions: " << m_fit_dimension << std::endl;
    }
*/

    //------- Print out whole grid --------------
    if(print_all){
        std::cout << intag<< "========================Below are detailed chi2 and coordinate values==========================" << std::endl;
        for(int i=0;i<vec_chi.size(); i++){
            int temp_poly_grid_index = i/m_flat_total_gridpoints;
            int temp_grid_index = i%m_flat_total_gridpoints;

            std::vector<double> temp_point = m_vec_grid[temp_grid_index];
            std::cout << "Point " << i << " Chi:" << vec_chi[i] << ", Coordinate: ";
            for(int j=0 ;j<temp_point.size();j++) std::cout << temp_point[j] << ", ";
            if(m_bool_poly_grid){
                std::vector<double> temp_poly_point = m_vec_poly_grid[temp_poly_grid_index];
                std::cout << "PolyFactor: ";
                for(int j=0; j<temp_poly_point.size();j++) std::cout << j+1<<"nd order: "<< temp_poly_point[j]<< ", ";
            }
            std::cout<< std::endl;
        }
        std::cout <<intag << "=========================End of detailed chi info=============================================="<< std::endl;
    }
    return 0;
}

int SBNsinglephoton::WriteOutInfo(std::map<int, std::vector<double>>& inmap){
	std::map<int, std::vector<double>>::iterator itmap = inmap.begin();

	TFile* fout = new TFile(Form("%s_fit_output.root", tag.c_str()), "UPDATE");
	TVectorD v(1);
	v[0] = itmap->first;
	std::vector<double> vchi=itmap->second;

	fout->cd();
	v.Write("bf_index");
	fout->WriteObject(&vchi, "vector_chi");
	fout->Close();
	return 0;
}


int SBNsinglephoton::SetFlatFullFracCovarianceMatrix(double flat){

    m_full_fractional_covariance_matrix = new TMatrixT<double>(num_bins_total, num_bins_total);
    for(int i=0;i<num_bins_total; i++) (*m_full_fractional_covariance_matrix)(i,i) = flat*flat;

    return 0;	
}

int SBNsinglephoton::SetStatOnly(){
    if(is_verbose) std::cout<<"SBNsinglephoton::SetStatOnly\t||\t Setting covariance matrices to be ZERO!" << std::endl;
    m_chi->is_stat_only = true;
    m_full_fractional_covariance_matrix = new TMatrixT<double>(num_bins_total, num_bins_total);
    m_full_but_genie_fractional_covariance_matrix = new TMatrixT<double>(num_bins_total, num_bins_total);
    m_genie_fractional_covariance_matrix = new TMatrixT<double>(num_bins_total, num_bins_total);

    m_full_fractional_covariance_matrix->Zero();
    m_full_but_genie_fractional_covariance_matrix->Zero();
    m_genie_fractional_covariance_matrix->Zero();

    return 0;
}

int SBNsinglephoton::ModifyCV(double infactor){

    std::vector<double> temp_constrain_param;
    for(auto const &dim:m_grid.f_dimensions){
        if(dim.f_has_constrain){
            temp_constrain_param.push_back(dim.f_constrain_value);
        }
    }

    if(m_bool_poly_grid){

        for(auto const &dim:m_poly_grid.f_dimensions){
            if(dim.f_has_constrain){
                temp_constrain_param.push_back(dim.f_constrain_value);
            }
        }
    }

    return this->ModifyCV(infactor, temp_constrain_param);
}

int SBNsinglephoton::ModifyCV(double infactor, std::vector<double> param){


    CheckCVLoad();

    std::cout<< "SBNsinglephoton::ModifyCV\t|| Start modifying CV spectrum" <<std::endl;
    int index = 0;
    double non_coh_factor=0;

    //apply flat scaling
    for(auto const &dim:m_grid.f_dimensions){
        if(dim.f_has_constrain  && (index < param.size()) ){
	    if(dim.GetName() == "NCPi0NotCoh") non_coh_factor = param[index];  //remember the scaling for NCpi0NotCoh, for the poly grid.
            m_cv_spectrum->Scale(dim.GetName(), param[index]);
            ++index;
        }
    }

    //energy/momentum dependent scaling
    //now that we reset negative weight to 1, this needs to be modified, tho it's probably never gonna be used
    if(m_bool_poly_grid){

   	std::ostringstream prescale_tag;
   	std::vector<double> temp_scale_parameter;
   	for(auto const &dim:m_poly_grid.f_dimensions){
      	     if(dim.f_has_constrain && (index < param.size())  ){   
	      	prescale_tag << "_" << std::fixed<< std::setprecision(3) << param[index];
	     	temp_scale_parameter.push_back(param[index]);
	     	++index;
             }
        }


   	if(temp_scale_parameter.size() != 0){
             m_cv_spectrum->Scale("NCPi0NotCoh", 0.0);

	     if(!m_file_open) this->OpenFiles();
	     if(!m_bool_cv_spectrum_generated) m_bool_cv_spectrum_generated=true;
             SBNspec spec_polygrid_part(xmlname, -1, false);
             this->ScaleSpectrum(&spec_polygrid_part, non_coh_factor, temp_scale_parameter);
             m_cv_spectrum->Add(&spec_polygrid_part); 
             //check if a file exits
             //if(gSystem->AccessPathName(temp_filename.c_str()))
   	}
    }

    

    if(m_bool_data_spectrum_loaded){
        m_cv_spectrum->Scale("NCDeltaLEE", (infactor-1)*0.5);
        m_chi->DrawComparisonIndividual(*m_cv_spectrum, *m_data_spectrum, "ModifiedCV_Data");
    }

    m_bool_modify_cv = true;
    m_cv_delta_scaling = infactor;

    //m_cv_spectrum->WriteOut("ModifiedCV");

    return 0;
}

SBNspec SBNsinglephoton::GeneratePointSpectra(int np){
    int m_poly_grid_index = np/m_flat_total_gridpoints;
    int m_grid_index = np%m_flat_total_gridpoints;

    //SBNspec spec_rgrid_part(tag+"_CV.SBNspec.root", xmlname, false);
    SBNspec spec_rgrid_part = *m_cv_spectrum;

    std::vector<double> temp_p_grid = m_vec_grid[m_grid_index];
    for(int i=0;i<m_grid.f_num_dimensions; i++){
        if(m_grid.f_dimensions[i].GetName() == "NCDeltaLEE" && temp_p_grid[i]<0 ){
            spec_rgrid_part.Scale("NCDeltaLEE", 0.0);
            spec_rgrid_part.Scale("NCDelta", 1+temp_p_grid[i]*2.0);
        }
        else
            spec_rgrid_part.Scale(m_grid.f_dimensions[i].GetName(), temp_p_grid[i]);
            //spec_rgrid_part.ScaleAll(temp_p_grid[i]);
    }

    if(m_bool_poly_grid){
        double non_coh_factor;   
        for(int i=0;i<m_grid.f_num_dimensions; i++){
            if(m_grid.f_dimensions[i].GetName() == "NCPi0NotCoh"){
                non_coh_factor=temp_p_grid[i];
                spec_rgrid_part.Scale("NCPi0NotCoh", 0.0);
            }
        }

	std::vector<double> temp_p_polygrid=m_vec_poly_grid[m_poly_grid_index];

	if(!m_file_open) this->OpenFiles();
	if(!m_bool_cv_spectrum_generated) m_bool_cv_spectrum_generated=true;
	SBNspec spec_polygrid_part(xmlname, -1, false);
	this->ScaleSpectrum(&spec_polygrid_part, non_coh_factor, temp_p_polygrid);

	spec_rgrid_part.Add(&spec_polygrid_part);
    }


    return spec_rgrid_part;
}


int SBNsinglephoton::PoissonFluctuation(SBNspec *inspec){
    std::cout << "SBNsinglephoton::PoissonFluctuation\t || Poisson fluctuation the spectrum!!" << std::endl;
    inspec->CollapseVector();

    inspec->ScalePoisson();
    inspec->CollapseVector();
    inspec->CalcErrorVector();

    return 0;
}


double SBNsinglephoton::CalcChi(bool use_cnp){

    m_cv_spectrum->CollapseVector();
    m_cv_spectrum->CalcErrorVector();
    m_data_spectrum->CollapseVector();
    //get the collapsed systematic covariance matrix
    TMatrixT<double> Mtemp = m_chi->FillSystMatrix(*m_full_fractional_covariance_matrix, m_cv_spectrum->full_vector, m_cv_spectrum->full_err_vector, true);

    TMatrixT<double> Mtotal, InvertM;
    if(!use_cnp){
	Mtotal = m_chi->AddStatMatrixCNP(&Mtemp, m_cv_spectrum->collapsed_vector, m_data_spectrum->collapsed_vector);
    }else{
	Mtotal = m_chi->AddStatMatrix(&Mtemp, m_cv_spectrum->collapsed_vector);
    }

    InvertM = m_chi->InvertMatrix(Mtotal);
    return m_chi->CalcChi(InvertM, m_cv_spectrum->collapsed_vector, m_data_spectrum->collapsed_vector, true);

}

void SBNsinglephoton::SetInterpolationNumber(int in){
	m_interpolation_number = in;
}

void SBNsinglephoton::MaskScaleFactor(double& factor){
    if(factor < 0) factor = 0.0;
}

void SBNsinglephoton::DrawMarginalizedChi(std::vector<TH1D*> vec_hist, std::string intag, bool sensitivity, std::vector<double> fit_result){

    fin->cd();

    std::vector<double> chi_contour{1.0, 2.71, 6.63};
    std::vector<std::string> contour_label{"1#sigma", "90%", "99%"};
    std::vector<int> line_style{ 2, 3,4}; //corresponds to Solid, Dashed, Dotted respectively
    int line_color = 593, bf_color = 625; //625=kRed -7, 593=kBlue-7

    //sensitivity text band
    TText* txt = new TText(0.5, 0.5, "SENSITIVITY");
    txt->SetTextAlign(11);
    txt->SetTextAngle(30);
    txt->SetTextSize(0.2);
    txt->SetTextColorAlpha(kRed -10, 0.25);


    //text of fitting information
    TLatex* lat = new TLatex();
    lat->SetTextSize(0.043);
    lat->SetTextFont(42);
    lat->SetTextAlign(13); //align at the top
    std::string plot_text;
    if(fit_result.size()){
	   std::ostringstream ssBF, ssBFchi, ssCVchi;;
	   ssBF << std::setprecision(2) << std::fixed << fit_result.at(0);	   
	   ssBFchi << std::setprecision(2) << std::fixed << fit_result.at(1);
	   //ssCVchi << std::setprecision(2) << std::fixed << fit_result.at(2);
	   ssCVchi << std::setprecision(2) << std::fixed << 100*(1-TMath::Prob(fit_result.at(2),1));
	   plot_text = "#splitline{#color[625]{#bf{BF point}}: "+ ssBF.str() + "}{#splitline{#chi^{2}/ndof (data|BF): " + ssBFchi.str() + "/"+ std::to_string(num_bins_total_compressed - m_fit_dimension) + "}{SM #Delta#chi^{2} is at " + ssCVchi.str() +"\% CL}}";
    }

    //an arrow indicating BF point
    TArrow* arr = new TArrow();
    arr->SetAngle(40);
    arr->SetLineWidth(3);
    arr->SetLineColor(bf_color);
    arr->SetFillColor(bf_color);

    for(auto &h:vec_hist){
	h->Write();
	TCanvas* c=new TCanvas((std::string("canvas_")+h->GetName()).c_str(), "c");
	TLegend* leg = new TLegend(0.15, 0.7, 0.6, 0.9);
	leg->SetFillStyle(0);
	leg->SetBorderSize(0);
	leg->SetTextSize(0.043);
	leg->SetHeader("C.L. lines (Wilks)", "");
	gStyle->SetOptStat(0);
        h->SetLineColor(line_color); h->SetLineWidth(3);
	h->Draw("hist");
        if(h->GetMaximum() > 60 ) h->SetMaximum(60.); 

	std::vector<TLine*> lines;
    	for(int i = 0; i != chi_contour.size(); ++i){
	    if( h->GetMaximum() > chi_contour[i]){
		lines.emplace_back(new TLine(h->GetXaxis()->GetXmin(), chi_contour[i], h->GetXaxis()->GetXmax(), chi_contour[i]));
		//lines.back()->SetLineColor(kGreen -3);
		lines.back()->SetLineColor(line_color);
		lines.back()->SetLineWidth(2);
		lines.back()->SetLineStyle(line_style[i]);
		lines.back()->Draw("same");
		leg->AddEntry(lines.back(), contour_label[i].c_str(),"L");	
	    }else 
		 break;
 	}	
	leg->Draw();
	if(sensitivity) txt->Draw();
	if(fit_result.size()){
		lat->DrawLatex(h->GetXaxis()->GetXmax()*0.5, h->GetMaximum()*0.98, plot_text.c_str());
		arr->DrawArrow(fit_result[0], h->GetMaximum()*0.3, fit_result[0], 0, 0.05, "|>");
	}
	c->Update();
	c->Write();
	c->SaveAs((intag+"_"+h->GetName()+".pdf").c_str(), "pdf");

   }

}


void SBNsinglephoton::LocateCVGlobalIndex(){

    otag="SBNsinglephoton::LocateCVGlobalIndex\t|| ";

    std::vector<double> central_value_flat(m_grid.f_num_dimensions, 1.0);  //flat grid point for GENIE central value
    for(int i=0; i != m_grid.f_num_dimensions; ++i){
        if(m_grid.f_dimensions[i].GetName().find("NCDeltaLEE") != std::string::npos)  //for LEE component, the scale factor for CV is 0
                central_value_flat[i]= 0.0;
    }

    for(int j=0; j<m_flat_total_gridpoints; ++j){

        std::vector<double> jpoint = m_vec_grid[j];
	if(jpoint == central_value_flat){
	    if(m_bool_poly_grid){
		std::vector<double> central_value_poly(m_poly_grid.f_num_dimensions, 0.0);  //poly grid point for GENIE central value
		for(int i=0; i<m_poly_total_gridpoints; ++i){
		    if(m_vec_poly_grid.at(i) == central_value_poly){
			m_cv_spec_index = i*m_flat_total_gridpoints + j;
    			std::cout << otag << "Locate the index of GENIE CV: " << m_cv_spec_index << std::endl;
			return;
		    }
		}

	    }else{
		m_cv_spec_index = j;
    		std::cout << otag << "Locate the index of GENIE CV: " << m_cv_spec_index << std::endl;
		return;
	    }
	}	
    } //loop over flat grid
    
}

void SBNsinglephoton::SetupAsimovDataset(const std::string &in_string){
    asimov_dataset = in_string;
}

void SBNsinglephoton::CheckCVLoad(){

    if(!m_bool_cv_spectrum_loaded){
        std::cout << "SBNsinglephoton::CheckCVLoad\t|| CV spec hasn't been loaded yet, load CV..." << std::endl;
        this->LoadCV();
    }
}

void SBNsinglephoton::CheckNumScaledSpectra(){
    if(m_scaled_spec_grid.size() != m_total_gridpoints){
        std::cerr << "SBNsinglephoton::CheckNumScaledSpectra\t|| ERROR!! # of scaled spectra: "<<m_scaled_spec_grid.size()<<" does NOT match grid size: " << m_total_gridpoints << " !!" <<std::endl;
        exit(EXIT_FAILURE);
    }

}

void SBNsinglephoton::CheckDataLoad(){

    otag="SBNsinglephoton::CheckDataLoad\t|| ";
    if(!m_bool_data_spectrum_loaded){
	std::cout << otag << "WARNING!! Data spec hasn't been loaded, will do SAFE-MODE, aka sensitivity study instead!" << std::endl;

	// interpret the asimov dataset configuration
	std::vector<std::string> asimov_config_channels;
	std::vector<double> asimov_config_scalings;
	std::string channel_delimiter = ";", value_delimiter=",";
	size_t pos_front=0, pos_end = asimov_dataset.find(channel_delimiter, pos_front);
	while(1){
	    
	    std::string individual_channel_info;
	    if(pos_end != std::string::npos) individual_channel_info = asimov_dataset.substr(pos_front, pos_end - pos_front);
	    else individual_channel_info = asimov_dataset.substr(pos_front);

	    size_t pos_middle = individual_channel_info.find(value_delimiter);
	    asimov_config_channels.push_back(individual_channel_info.substr(0, pos_middle));
	    asimov_config_scalings.push_back(std::stod(individual_channel_info.substr(pos_middle + value_delimiter.length())));

	    if(pos_end == std::string::npos) break;
	    pos_front = pos_end + channel_delimiter.length();
	    pos_end = asimov_dataset.find(channel_delimiter, pos_front);
	}


	//prints out
	std::ostringstream writeout_tag;
	std::cout << otag << "Asimov dataset is configured as such: "; 
	for(size_t i = 0; i != asimov_config_channels.size(); ++ i){
	    writeout_tag << "_" << asimov_config_channels[i] << "_" << std::fixed<< std::setprecision(1) << asimov_config_scalings[i];
	    std::cout << asimov_config_channels[i] << ": " << asimov_config_scalings[i] << " ";
	}
        std::cout << std::endl; 

	m_data_spectrum = new SBNspec(tag+"_CV.SBNspec.root", xmlname, false);
	for(size_t i = 0; i != asimov_config_channels.size(); ++ i)
	    m_data_spectrum->Scale(asimov_config_channels[i], asimov_config_scalings[i]);
	m_data_spectrum->CollapseVector();
	//this->PoissonFluctuation(m_data_spectrum);
	m_data_spectrum->WriteOut("ToyData"+writeout_tag.str());

    }else{
	std::cout << otag << "Data spec is loaded! Will be fitting to REAL data!" << std::endl;
    }    
}


void SBNsinglephoton::PrintCollapsedFractionalMatrix(const TMatrixT<double> &frac, const SBNspec &spec){

     otag="SBNsinglephoton::PrintCollapsedFractionalMatrix\t|| ";
     TMatrixT<double> collapse_matrix = m_chi->FillSystMatrix(frac, spec.full_vector, spec.full_err_vector, true);

     size_t Ndim = collapse_matrix.GetNcols();
     std::cout << otag << "GENIE CV" << std::endl;

     //3 significant digits
     std::cout.precision(3);    

     // set 4 decimal places
     //std::cout.precision(4);
     //std::cout.setf(std::ios::fixed, std::ios::floatfield);
     for(size_t i =0; i != Ndim; ++i){
	std::cout << otag;
	for(size_t j=0; j!=Ndim; ++j){
	    double element = collapse_matrix(i,j)/(spec.collapsed_vector[i]*spec.collapsed_vector[j]);
	    std::cout << element << " ";
	}
	std::cout << std::endl;
     }

     return;
}

void SBNsinglephoton::PublicDataPrintOut(const std::vector<std::string>& subch, const std::vector<double>& scale){
    CheckCVLoad();
    CheckDataLoad();

    auto current_spec = *m_cv_spectrum;
    std::ofstream txt_file;

    //scale the MC spectrum
    std::ostringstream file_suffix;
    file_suffix.precision(2);
    for(int i =0; i!= subch.size(); ++i){
	current_spec.Scale(subch[i], scale[i]);
	file_suffix << "_" << subch[i] << "_" << std::fixed << scale[i];	
    }

  
    //calculate the collapsed systematic covariance matrix
    current_spec.CollapseVector();
    current_spec.CalcErrorVector();
    std::vector<double> current_spec_full = current_spec.full_vector;
    std::vector<double> current_spec_collapse = current_spec.collapsed_vector;
    TMatrixT<double> systematic_matrix = m_chi->FillSystMatrix(*m_full_fractional_covariance_matrix, current_spec.full_vector, current_spec.full_err_vector, true);


    // map each channel to correponding normalization bin-width
    std::map<std::string, double> channel_bin_norm_map = {{"1g1p", 0.1}, {"1g0p", 0.1}, {"2g1p", 0.075}, {"2g0p", 0.075}};

    //now start to look at each channel, and print out observed data, MC prediction and overall uncertainty at each channel
    for(int im = 0; im <mode_names.size(); ++im){
        for(int id = 0; id <detector_names.size(); ++id){

	    int channel_full_bin_start = num_bins_detector_block * id + num_bins_mode_block * im;
	    int channel_collapse_bin_start = num_bins_detector_block_compressed * id + num_bins_mode_block_compressed * im;

            for(int ic = 0; ic <channel_names.size(); ++ic){


		//open a new file
                std::string file_name = mode_names.at(im)+"_"+detector_names.at(id)+"_"+channel_names.at(ic) + file_suffix.str() + ".txt";
 		std::cout << file_name << std::endl;

		txt_file.open(file_name, std::ofstream::out | std::ofstream::trunc);
    		txt_file.setf(std::ios::left);

		// print out helper info
		txt_file << "Event Counts Per " << channel_bin_norm_map[channel_names[ic]] << " GeV\n" << std::endl;
 		txt_file << std::setw(15) << "LowBinEdge" << std::setw(15) << "HighBinEdge";
		for(auto& subch_name : subchannel_names[ic])
	 	    txt_file << std::setw(15) << subch_name;
		txt_file << std::setw(15) << "TotalCount" << std::setw(15) << "Uncertainty" << std::setw(15) <<  "Data" << std::setw(15) << "GausError" << std::setw(20) << "PoissonErrorUp" << std::setw(20) << "PoissonErrorLow" << std::endl;
		txt_file << std::setw(15) << "----------" << std::setw(15) <<"----------";
                for(auto& subch_name : subchannel_names[ic])
                    txt_file << std::setw(15) <<"----------";
                txt_file << std::setw(15) <<"----------" << std::setw(15) <<"----------" << std::setw(15) <<"----------" << std::setw(15) <<"----------" << std::setw(20) <<"-------------" << std::setw(20) <<"-------------" << std::endl;


		// bin edge and width of current channel	
		const std::vector<double>& channel_bin_edges = bin_edges[ic];
		const std::vector<double>& channel_bin_width = bin_widths[ic];

		// use TH1D to store observation, for extraction of Poisson error
		TH1D* hchan_data = new TH1D("hchan_data", "hchan_data", num_bins[ic], &channel_bin_edges[0]);
		hchan_data->SetBinErrorOption(TH1::kPoisson); 

		//iterate through each bin
		for(int bin = 0; bin != num_bins[ic]; ++bin){
		    

		    double bin_norm_factor = channel_bin_width[bin]/channel_bin_norm_map[channel_names[ic]];
 
		    //print out MC prediction and uncertainties.
		    txt_file << std::setw(15) << channel_bin_edges[bin] << std::setw(15) << channel_bin_edges[bin+1];
		    for(int isub = 0; isub != num_subchannels[ic]; ++isub){
			txt_file << std::setw(15) << current_spec_full[channel_full_bin_start + bin + num_bins[ic]*isub]/bin_norm_factor;  
		    }

		    txt_file << std::setw(15) << current_spec_collapse[channel_collapse_bin_start + bin]/bin_norm_factor << std::setw(15) << sqrt(systematic_matrix(channel_collapse_bin_start + bin, channel_collapse_bin_start + bin))/bin_norm_factor;

		    //pint out data-related quantities
		    hchan_data->SetBinContent(bin+1, m_data_spectrum->collapsed_vector[channel_collapse_bin_start + bin]);
		    txt_file << std::setw(15) << hchan_data->GetBinContent(bin+1)/bin_norm_factor << std::setw(15) << sqrt(hchan_data->GetBinContent(bin+1))/bin_norm_factor << std::setw(20) << hchan_data->GetBinErrorUp(bin+1)/bin_norm_factor << std::setw(20) << hchan_data->GetBinErrorLow(bin+1)/bin_norm_factor << std::endl ;


		} // end of iterating through each bin

		txt_file.close();


		//update the channel bin start index
		channel_full_bin_start += num_bins[ic] * num_subchannels[ic];
		channel_collapse_bin_start += num_bins[ic]; 
            } //channel loop
        } //detector loop

    } //mode loop

    return;
}
