#include "SBNsinglephoton.h"

using namespace sbn;


SBNsinglephoton::SBNsinglephoton(std::string xmlname, std::string intag, NGrid ingrid):SBNsinglephoton(xmlname, intag, ingrid,ingrid, false){}

SBNsinglephoton::SBNsinglephoton(std::string xmlname, std::string intag, NGrid ingrid, NGrid in_polygrid, bool has_polygrid): SBNconfig(xmlname), tag(intag), m_grid(ingrid), m_bool_poly_grid(has_polygrid){

    if(is_verbose) std::cout << "SBNsinglephoton::SBNsinglephoton\t|| Setup grid" << std::endl;
    m_vec_grid = m_grid.GetGrid();
    m_num_total_gridpoints = m_grid.f_num_total_points;

    m_poly_total_gridpoints = 1;
    if(m_bool_poly_grid){
	 m_poly_grid=in_polygrid;
	 m_vec_poly_grid = m_poly_grid.GetGrid();
         m_poly_total_gridpoints = m_poly_grid.f_num_total_points;
    }

    m_total_gridpoints = m_num_total_gridpoints*m_poly_total_gridpoints;
    m_max_number_iterations = 20;
    m_chi_min_convergance_tolerance = 0.001;
    m_bool_modify_cv = false;
    m_bool_cv_spectrum_generated = false;
    m_bool_cv_spectrum_loaded = false;
    m_bool_data_spectrum_loaded = false;


    m_full_fractional_covariance_matrix = NULL;
    m_full_but_genie_fractional_covariance_matrix = NULL;
    m_genie_fractional_covariance_matrix = NULL;
}

int SBNsinglephoton::SetPolyGrid(NGrid in_polygrid){
    if(is_verbose) std::cout << "SBNsinglephoton::SetPolyGrid\t|| Setup polynomial grid " << std::endl;
    m_bool_poly_grid = true;
    m_poly_grid=in_polygrid;
    m_vec_poly_grid = m_poly_grid.GetGrid();
    m_poly_total_gridpoints = m_poly_grid.f_num_total_points;

    m_total_gridpoints = m_num_total_gridpoints*m_poly_total_gridpoints;
 
    return 0;
}

int SBNsinglephoton::OpenFiles(){

    num_files = montecarlo_file.size();
    montecarlo_additional_weight.resize(num_files,1.0);
    montecarlo_additional_weight_formulas.resize(num_files);


    if(is_verbose) std::cout<< "SBNsinglephoton::OpenFiles\t|| Open the files"<< std::endl;
    for(auto &fn: montecarlo_file){
        files.push_back(new TFile(fn.c_str()));
        if(files.back()->IsZombie() || !files.back()->IsOpen()){
            std::cout<<"SBNsinglephoton::OpenFiles\t|| ERROR! Failed to open the file "<<fn<<std::endl;
            exit(EXIT_FAILURE);
        }
    }


    for(int i=0; i<montecarlo_name.size(); i++){
        std::cout<<"Getting TTree "<<montecarlo_name[i]<<" from file "<<montecarlo_file[i]<<std::endl;
        trees.push_back((TTree*)files.at(i)->Get(montecarlo_name.at(i).c_str()) );
        std::cout<<"--TTree has "<<trees.back()->GetEntries()<<" entries. "<<std::endl;
    }

   if(is_verbose) std::cout << "SBNsinglephoton::OpenFiles\t||Setup any friendtrees" << std::endl;
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
                trees.at(i)->SetBranchAddress( branch_variable->true_param_name.c_str(), branch_variable->GetTrueValue() );
            }
        }

        if(montecarlo_additional_weight_bool[i]){
            std::cout<<"Setting Additional weight of : "<< montecarlo_additional_weight_names[i].c_str()<<std::endl;
            montecarlo_additional_weight_formulas[i] =  new TTreeFormula(("a_w"+std::to_string(i)).c_str(),montecarlo_additional_weight_names[i].c_str(),trees[i]);
        }


    }

    if(is_verbose) std::cout << "SBNsinglephoton::OpenFiles\t|| Finish opening files and setting up TTrees " << std::endl; 

    return 0;
}


int SBNsinglephoton::CloseFiles(){
    if(is_verbose) std::cout<< "SBNsinglephoton::CloseFiles\t|| Closing TFiles..."<< std::endl;
    for(auto f: files){
            std::cout <<" TFile::Close() file=" << f->GetName() << " @" << f << std::endl;
            f->Close();
    }
    return 0;
}

int SBNsinglephoton::PreScaleSpectrum(std::string xmlname, double flat_factor, std::vector<double>& param){
    SBNspec tm(xmlname,-1,false);
    SBNspec temp_cv_spectrum = tm;
    SBNspec spec_prescale  = tm;   //prescaled spectrum

    std::cout<<"SBNsinglephoton::PreScaleSpectrum\t|| -----------------------------------------------\n";
    std::cout<<"SBNsinglephoton::PreScaleSpectrum\t|| -----------------------------------------------\n";

    for(int j=0;j<num_files;j++){


        for(int i=0; i< std::min(  montecarlo_maxevents.at(j)  ,nentries.at(j)); i++){
            trees.at(j)->GetEntry(i);

            if(i%100==0) std::cout<<"SBNsinglephoton::PreScaleSpectrum\t|| On event: "<<i<<" of "<<nentries[j]<<" from File: "<<montecarlo_file[j]<<std::endl;

            double global_weight = 1.0;
            if( montecarlo_additional_weight_bool[j]){
                    montecarlo_additional_weight_formulas[j]->GetNdata();
                    global_weight = montecarlo_additional_weight_formulas[j]->EvalInstance();
            };//this will be 1.0 unless specified
            global_weight = global_weight*montecarlo_scale[j];



            if(std::isinf(global_weight) || global_weight != global_weight){
                std::cout<<"SBNsinglephoton::PreScaleSpectrum\t|| ERROR  error @ "<<i<<" in File "<<montecarlo_file.at(j)<<" as its either inf/nan: "<<global_weight<<std::endl;
                exit(EXIT_FAILURE);
            }



            for(int t=0; t<branch_variables[j].size();t++){
                    const auto branch_variable = branch_variables[j][t];
                    int ih = temp_cv_spectrum.map_hist.at(branch_variable->associated_hist);

		    //grab reco info and determine which reco bin it belongs to
                    branch_variable->GetFormula()->GetNdata();
                    double reco_var = branch_variable->GetFormula()->EvalInstance();
                    int reco_bin = temp_cv_spectrum.GetGlobalBinNumber(reco_var,ih);

                    if(branch_variables[j][t]->GetOscillate()){
			//truth info
  			double true_var = *(static_cast<double*>(branch_variables[j][t]->GetTrueValue()));

                        double prescale_factor = this->ScaleFactor(true_var, flat_factor, param);

                        spec_prescale.hist[ih].Fill(reco_var, global_weight*prescale_factor);
         		if(!m_bool_cv_spectrum_generated)  temp_cv_spectrum.hist[ih].Fill(reco_var,global_weight);
 		    }else{
                        if(!m_bool_cv_spectrum_generated) temp_cv_spectrum.hist[ih].Fill(reco_var,global_weight);
  		    }
            }
        } //end of entry loop
    } // end of nfile loop

 
    if(is_verbose) std::cout<< "SBNsinglephoton::PreScaleSpectrum\t||Write out spectra" << std::endl;
    std::ostringstream prescale_tag;
    if(flat_factor != 0) prescale_tag << "Flat_" << std::fixed<< std::setprecision(3) << flat_factor << "_PreScaled";
    else prescale_tag << "PreScaled";
    //generate tag for prescaled spectra
    for(int i=0; i< param.size(); i++){
	prescale_tag << "_" << std::fixed<< std::setprecision(3) << param[i]; 
    }   
    if(param.size() != 0) spec_prescale.WriteOut(tag+"_"+prescale_tag.str());
    //if(param.size() != 0) spec_prescale.WriteOut(tag+"_PreScaled"+prescale_tag.str());

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
    this->OpenFiles();
    if(!m_bool_poly_grid){
	std::cout << "SBNsinglephoton::GeneratePreScaledSpectra\t|| No grid for polynomial scaling present------------------\n"<< std::endl;
	std::cout << "SBNsinglephoton::GeneratePreScaledSpectra\t|| Gonna generate only CV, then exiting GeneratePreScaledSpectra ..." << std::endl;
	std::vector<double> empty{};
	this->PreScaleSpectrum(xmlname, empty);
    }else{
	
    	//m_vec_poly_grid = m_poly_grid.GetGrid();
    	//m_poly_total_gridpoints = m_poly_grid.f_num_total_points;

    	for(int i=0; i< m_poly_total_gridpoints; i++){
		std::vector<double> point = m_vec_poly_grid[i];
		this->PreScaleSpectrum(xmlname, point);	
    	}
    }

    this->CloseFiles();

    return 0;
}


int SBNsinglephoton::LoadSpectraApplyFullScaling(){
   //	if(!m_bool_cv_spectrum_generated){
   //		std::cout << "SBNsinglephoton::LoadSpectraApplyFullScaling\t|| CV spectrum hasn't been generated yet" << std::endl;
   //		exit(EXIT_FAILURE);
   //	}

	SBNspec temp_cv((tag+"_CV.SBNspec.root").c_str(), this->xmlname, false);
	temp_cv.Scale("BNBext", 0.0);
	temp_cv.Scale("Dirt", 0.0);

	//m_scaled_spec_grid.clear();
	m_scaled_spec_grid.resize(m_total_gridpoints); 


	int ip_processd = 0; 
	if(m_bool_poly_grid){
	    std::cout<<"SBNsinglephoton::LoadSpectraApplyFullScaling\t|| Grab pre-scaled spectra and build the final spectra!" << std::endl; 
 	    //loop over polynomial grid
	    for(int i=0;i< m_poly_total_gridpoints; i++){
		std::vector<double> ipoint = m_vec_poly_grid[i];

		//first initilize SBNspec with pre-scaled root file
		std::ostringstream prescale_tag;
        	for(int ip=0; ip< ipoint.size(); ip++){
                	prescale_tag << "_" << std::fixed<< std::setprecision(3) << ipoint[ip];
                }	
		std::string full_filename = this->tag+"_PreScaled"+prescale_tag.str() + ".SBNspec.root";
		SBNspec temp_prescaled(full_filename, this->xmlname);	

		//loop over regular grid
		for(int j=0; j< m_num_total_gridpoints;j++){
			if(is_verbose) std::cout << "On Point " << i*m_num_total_gridpoints+j <<"/"<<m_total_gridpoints << std::endl;

			std::vector<double> jpoint = m_vec_grid[j];
		
			//m_scaled_spec_grid.push_back(new SBNspec((tag+"_CV.SBNspec.root").c_str(), this->xmlname, false));//start with genie CV
			m_scaled_spec_grid[ip_processd] = temp_cv; //start with genie CV
			for(int jp=0; jp< jpoint.size(); jp++){
				//m_scaled_spec_grid.back()->Scale(m_grid.f_dimensions[jp].f_name, jpoint[jp]); // scale corresponding subchannel
				m_scaled_spec_grid[ip_processd].Scale(m_grid.f_dimensions[jp].f_name, jpoint[jp]); // scale corresponding subchannel
			}
			
			//m_scaled_spec_grid.back()->Add(&temp_prescaled); //here we get the scaled spectra at final stage!!!
			m_scaled_spec_grid[ip_processd].Add(&temp_prescaled); //here we get the scaled spectra at final stage!!!

			// do not want to keep LEE in it
			//m_scaled_spec_grid.back()->Scale("NCDeltaRadOverlayLEE", 0.0);
			//m_scaled_spec_grid[ip_processd].Scale("NCDeltaRadOverlayLEE", 0.0);
			ip_processd++;

		}//end loop over regular grid

	     }// end loop over poly grid
	}
	else{
	    std::cout<<"SBNsinglephoton::LoadSpectraApplyFullScaling\t|| Scale the spectra for the whole grid " << std::endl;
	     for(int j=0; j< m_num_total_gridpoints;j++){
                        if(is_verbose && (j%1000 ==0)) std::cout << "On Point " << j <<"/"<<m_total_gridpoints << std::endl;

			std::vector<double> jpoint = m_vec_grid[j];

                        //m_scaled_spec_grid.push_back(new SBNspec((tag+"_CV.SBNspec.root").c_str(), this->xmlname, false));//start with genie CV
			m_scaled_spec_grid[ip_processd] = temp_cv; //start with genie CV
			for(int jp=0; jp< jpoint.size(); jp++){
                                //m_scaled_spec_grid.back()->Scale(m_grid.f_dimensions[jp].f_name, jpoint[jp]); // scale corresponding subchannel
                                m_scaled_spec_grid[ip_processd].Scale(m_grid.f_dimensions[jp].f_name, jpoint[jp]); // scale corresponding subchannel
			}

	     		// do not want to keep LEE in it
			//m_scaled_spec_grid.back()->Scale("NCDeltaRadOverlayLEE", 0.0);
			//m_scaled_spec_grid[ip_processd].Scale("NCDeltaRadOverlayLEE", 0.0);
	
			ip_processd++;
	    }
	}

	return 0;
}



int SBNsinglephoton::CalcChiGridScanShapeOnlyFit(){
    if(!m_bool_cv_spectrum_loaded){
	std::cout << "SBNsinglephoton::CalcChiGridScanShapeOnlyFit\t|| CV spec hasn't been loaded yet, load CV..." << std::endl;
	this->LoadCV();
    }	

    if(!m_bool_data_spectrum_loaded){
	std::cout << "SBNsinglephoton::CalcChiGridScanShapeOnlyFit\t|| WARNING!! Data spec hasn't been loaded, will do a sensitivity study instead!" << std::endl;
	m_data_spectrum = new SBNspec();
	*m_data_spectrum = *m_cv_spectrum;
	//m_data_spectrum->Scale("NCDeltaRadOverlayLEE", 0.0);
	m_data_spectrum->CollapseVector();
    }    


    if(m_scaled_spec_grid.size() != m_total_gridpoints){
	std::cout << "SBNsinglephoton::CalcChiGridScanShapeOnlyFit\t|| ERROR!! # of scaled spectra: "<<m_scaled_spec_grid.size()<<" does NOT match grid size: " << m_total_gridpoints << " !!" <<std::endl;
	exit(EXIT_FAILURE);
    }

	std::cout << "SBNsinglephoton::CalcChiGridScanShapeOnlyFit\t" << __LINE__ << std::endl;

    SBNspec background_spectrum;
    SBNspec last_best_spectrum;
    double best_chi, last_best_chi;
    int best_point, last_best_point;
    std::vector<double> vec_chi, last_vec_chi;
    TFile* fout = new TFile("NCpi0_normalization_fit_output.root", "recreate");  //save matrix plot etc.

    m_chi = new SBNchi(this->xmlname);
    TMatrixT<double> full_systematic_covariance(num_bins_total, num_bins_total);
    TMatrixT<double> genie_systematic_matrix(num_bins_total, num_bins_total);
    TMatrixT<double> collapsed_full_systematic_matrix(num_bins_total_compressed, num_bins_total_compressed);
    TMatrixT<double> total_covariance_matrix(num_bins_total_compressed, num_bins_total_compressed);
    TMatrixT<double> inversed_total_covariance_matrix(num_bins_total_compressed, num_bins_total_compressed);

    for(int n_iter =0; n_iter < m_max_number_iterations; n_iter ++){
	std::cout << "SBNsinglephoton::CalcChiGridScanShapeOnlyFit\t|| On fit iteration "<< n_iter << std::endl;
	//reset everything at the biginning of each iteration
	best_chi = DBL_MAX;
	vec_chi.clear();


	if(n_iter == 0){
		last_best_spectrum = *m_cv_spectrum;
		//last_best_spectrum.Scale("NCDeltaRadOverlayLEE", 0.0);
	}else{
		last_best_spectrum = m_scaled_spec_grid[last_best_point];		
	}

	background_spectrum = last_best_spectrum;

	//============================Calculate covariance matrix and its invert====================================/
	//full systematic covariance matrix, except genie uncertainty
	full_systematic_covariance = m_chi->FillSystMatrix(*m_full_but_genie_fractional_covariance_matrix, last_best_spectrum.full_vector);
	//calculate the shape only covariance matrix for genie uncertainty, to get rid of normalization uncertainty
	for(int i=0; i<m_grid.f_num_dimensions; i++){
	   SBNspec temp_comp = last_best_spectrum;
	   temp_comp.Keep( m_grid.f_dimensions[i].f_name, 1.0 );
	   genie_systematic_matrix = m_chi->FillSystMatrix(*m_genie_fractional_covariance_matrix, temp_comp.full_vector);
	   fout->cd();
	   genie_systematic_matrix.Write(Form("full_genie_%s_%d", m_grid.f_dimensions[i].f_name.c_str(), n_iter));	

	   genie_systematic_matrix = m_chi->CalcShapeOnlyCovarianceMatrix(*m_genie_fractional_covariance_matrix, &temp_comp, &temp_comp);
	   genie_systematic_matrix.Write(Form("ShapeOnly_genie_%s_%d", m_grid.f_dimensions[i].f_name.c_str(), n_iter));	
	   full_systematic_covariance += genie_systematic_matrix;

	   background_spectrum.Scale(m_grid.f_dimensions[i].f_name, 0.0);   
	}
	//add genie uncertainties for other subchannels to total covariance matrix
	genie_systematic_matrix = m_chi->FillSystMatrix(*m_genie_fractional_covariance_matrix, background_spectrum.full_vector);
	full_systematic_covariance += genie_systematic_matrix;

	m_chi->CollapseModes(full_systematic_covariance, collapsed_full_systematic_matrix);

	//use CNP statistical covariance matrix
	total_covariance_matrix = m_chi->AddStatMatrixCNP(&collapsed_full_systematic_matrix, last_best_spectrum.collapsed_vector, m_data_spectrum->collapsed_vector);
	//std::cout << "SBNsinglephoton::CalcChiGridScanShapeOnlyFit\t||check " << __LINE__ << std::endl;

	inversed_total_covariance_matrix= m_chi->InvertMatrix(total_covariance_matrix);

	fout->cd();
	genie_systematic_matrix.Write(Form("full_genie_otherbkg_%d", n_iter));
	full_systematic_covariance.Write(Form("syst_uncollapsed_matrix_%d", n_iter));
	collapsed_full_systematic_matrix.Write(Form("syst_collapsed_matrix_%d", n_iter));
	total_covariance_matrix.Write(Form("total_collapsed_matrix_%d", n_iter));

	if(is_verbose && m_bool_data_spectrum_loaded) last_best_spectrum.CompareSBNspecs(collapsed_full_systematic_matrix, m_data_spectrum, tag+"_Iter_"+std::to_string(n_iter));	
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

	if(is_verbose) std::cout << "SBNsinglephoton::CalcChiGridScanShapeOnlyFit\t|| chi2 minimum :" << best_chi << "at point " << best_point << "/" << m_total_gridpoints << std::endl;

	if(n_iter != 0){
	   if(fabs(best_chi - last_best_chi) < m_chi_min_convergance_tolerance){
		std::cout << "SBNsinglephoton::CalcChiGridScanShapeOnlyFit\t|| chi2 has converged with best chi2 value " << best_chi << " at point " << best_point << "/" << m_total_gridpoints << std::endl;
		break;
	   }
	}

	last_best_chi = best_chi;
	last_best_point = best_point;
	last_vec_chi = vec_chi;
    }//end loop for iteration

    fout->Close();
	
    m_map={{best_point, vec_chi}};
    this->PrintOutFitInfo(m_map, "SBNsinglephoton::CalcChiGridScanShapeOnlyFit\t||"+tag, true);
    //std::cout << "SBNsinglephoton::CalcChiGridScanShapeOnlyFit\t||check " << __LINE__ << std::endl;

    return 0;
}



int SBNsinglephoton::CalcChiGridScan(){
    if(!m_bool_cv_spectrum_loaded){
	std::cout << "SBNsinglephoton::CalcChiGridScan\t|| CV spec hasn't been loaded yet, load CV..." << std::endl;
	this->LoadCV();
    }	

    if(!m_bool_data_spectrum_loaded){
	std::cout << "SBNsinglephoton::CalcChiGridScan\t|| WARNING!! Data spec hasn't been loaded, will do a sensitivity study instead!" << std::endl;
	m_data_spectrum = new SBNspec();
	*m_data_spectrum = *m_cv_spectrum;
	//m_data_spectrum->Scale("NCDeltaRadOverlayLEE", 0.0);
	if(m_bool_modify_cv) m_data_spectrum->Scale("NCDeltaRadOverlayLEE", (m_cv_delta_scaling-1)*0.5);
	m_data_spectrum->CollapseVector();
    }else{
	m_cv_spectrum->CompareSBNspecs(m_data_spectrum, tag+"_CVvsData_NoErrorBar");
    }    


    if(m_scaled_spec_grid.size() != m_total_gridpoints){
	std::cout << "SBNsinglephoton::CalcChiGridScan\t|| ERROR!! # of scaled spectra: "<<m_scaled_spec_grid.size()<<" does NOT match grid size: " << m_total_gridpoints << " !!" <<std::endl;
	exit(EXIT_FAILURE);
    }


    SBNspec last_best_spectrum;
    double best_chi, last_best_chi;
    int best_point, last_best_point;
    std::vector<double> vec_chi, last_vec_chi;
    TFile* fout = new TFile("NCdelta_fit_output.root", "recreate");  //save matrix plot etc.

    m_chi = new SBNchi(this->xmlname);
    TMatrixT<double> collapsed_full_systematic_matrix(num_bins_total_compressed, num_bins_total_compressed);
    TMatrixT<double> total_covariance_matrix(num_bins_total_compressed, num_bins_total_compressed);
    TMatrixT<double> inversed_total_covariance_matrix(num_bins_total_compressed, num_bins_total_compressed);

    for(int n_iter =0; n_iter < m_max_number_iterations; n_iter ++){
	std::cout << "SBNsinglephoton::CalcChiGridScan\t|| On fit iteration "<< n_iter << std::endl;
	//reset everything at the biginning of each iteration
	best_chi = DBL_MAX;
	vec_chi.clear();


	if(n_iter == 0){
		last_best_spectrum = *m_cv_spectrum;
		//last_best_spectrum.Scale("NCDeltaRadOverlayLEE", 0.0);
		last_best_spectrum.Scale("BNBext", 0.0);
                last_best_spectrum.Scale("Dirt", 0.0);
	}else{
		last_best_spectrum = m_scaled_spec_grid[last_best_point];		
	}


	//============================Calculate covariance matrix and its invert====================================/
	//full systematic covariance matrix, except genie uncertainty
	collapsed_full_systematic_matrix = m_chi->FillSystMatrix(*m_full_fractional_covariance_matrix, last_best_spectrum.full_vector, true);

	//use CNP statistical covariance matrix
	total_covariance_matrix = m_chi->AddStatMatrixCNP(&collapsed_full_systematic_matrix, last_best_spectrum.collapsed_vector, m_data_spectrum->collapsed_vector);

	inversed_total_covariance_matrix= m_chi->InvertMatrix(total_covariance_matrix);

	fout->cd();
	collapsed_full_systematic_matrix.Write(Form("syst_collapsed_matrix_%d", n_iter));
	total_covariance_matrix.Write(Form("total_collapsed_matrix_%d", n_iter));

	if(is_verbose && m_bool_data_spectrum_loaded) last_best_spectrum.CompareSBNspecs(collapsed_full_systematic_matrix, m_data_spectrum, tag+"_Iter_"+std::to_string(n_iter));	
	//============================Done calculating covariance matrix ============================================/


	for(int i=0 ;i <m_total_gridpoints; i++){
	   if(i%1000 ==0) std::cout<< "On Point " << i << "/" << m_total_gridpoints << std::endl;
	   double temp_chi = m_chi->CalcChi(inversed_total_covariance_matrix, m_scaled_spec_grid[i].collapsed_vector, m_data_spectrum->collapsed_vector);
	   vec_chi.push_back(temp_chi);

	   if(temp_chi < best_chi){
		best_chi = temp_chi;
	        best_point = i;
	   }
	}

	if(is_verbose) std::cout << "SBNsinglephoton::CalcChiGridScan\t|| chi2 minimum :" << best_chi << "at point " << best_point << "/" << m_total_gridpoints << std::endl;

	if(n_iter != 0){
	   if(fabs(best_chi - last_best_chi) < m_chi_min_convergance_tolerance){
		std::cout << "SBNsinglephoton::CalcChiGridScan\t|| chi2 has converged with best chi2 value " << best_chi << " at point " << best_point << "/" << m_total_gridpoints << std::endl;
		break;
	   }
	}

	last_best_chi = best_chi;
	last_best_point = best_point;
	last_vec_chi = vec_chi;
    }//end loop for iteration

    fout->Close();

    //best-fit vs data comparison	
    if(is_verbose && m_bool_data_spectrum_loaded){
        collapsed_full_systematic_matrix = m_chi->FillSystMatrix(*m_full_fractional_covariance_matrix, m_scaled_spec_grid[best_point].full_vector, true);
	SBNspec temp_best_spec = this->GeneratePointSpectra(best_point);
        temp_best_spec.CompareSBNspecs(collapsed_full_systematic_matrix, m_data_spectrum, tag+"_BFvsData");
    }

    m_map={{best_point, vec_chi}};
    this->PrintOutFitInfo(m_map, "SBNsinglephoton::CalcChiGridScan\t|| "+tag, true);

    return 0;
}


int SBNsinglephoton::LoadCV(){
    if(is_verbose) std::cout << "SBNsinglephoton::LoadCV\t|| Setup CV spectrum" << std::endl;
    m_cv_spectrum = new SBNspec((tag+"_CV.SBNspec.root").c_str(), xmlname);
    m_cv_spectrum->CollapseVector();
    m_bool_cv_spectrum_loaded = true;

    return 0;
}

int SBNsinglephoton::LoadData(std::string filename){
    if(is_verbose) std::cout << "SBNsinglephoton::LoadCV\t|| Load data spectrum from file: " << filename << std::endl;
    m_data_spectrum = new SBNspec(filename.c_str(), xmlname);
    m_bool_data_spectrum_loaded = true;
    m_data_spectrum->CollapseVector(); 
    return 0;
}


int SBNsinglephoton::SetFullFractionalCovarianceMatrix(std::string filename, std::string matrix_name){

	TFile* f_syst = new TFile(filename.c_str(), "read");
	m_full_fractional_covariance_matrix = (TMatrixT<double>*)f_syst->Get(matrix_name.c_str());

	if(m_full_fractional_covariance_matrix->GetNcols() != m_full_fractional_covariance_matrix->GetNrows()){
		std::cout << "SBNsinglephoton::SetFullFractionalCovarianceMatrix\t|| Matrix provided is not sysmetric" << std::endl;
		exit(EXIT_FAILURE);		
	}

	this->RemoveNan(m_full_fractional_covariance_matrix);

/*	for(int i=0; i<m_full_fractional_covariance_matrix->GetNcols(); i++){
	    for(int j=0; j<m_full_fractional_covariance_matrix->GetNcols(); j++){
		if(std::isnan((*m_full_fractional_covariance_matrix)(i,j))) (*m_full_fractional_covariance_matrix)(i,j) =0.0;
	    }	
	}
*/
	f_syst->Close();
	return 0;
}

int SBNsinglephoton::SetGenieFractionalCovarianceMatrix(std::string filename){

	//check with Gray, are matrices added directly before checking the nan's?
	TFile* f_syst = new TFile(filename.c_str(), "read");
	m_genie_fractional_covariance_matrix  = (TMatrixT<double>*)f_syst->Get("individualDir/All_UBGenie_frac_covariance");
        *m_genie_fractional_covariance_matrix += *((TMatrixT<double>*)f_syst->Get("individualDir/AxFFCCQEshape_UBGenie_frac_covariance"));
        *m_genie_fractional_covariance_matrix += *((TMatrixT<double>*)f_syst->Get("individualDir/DecayAngMEC_UBGenie_frac_covariance"));
        *m_genie_fractional_covariance_matrix += *((TMatrixT<double>*)f_syst->Get("individualDir/NormCCCOH_UBGenie_frac_covariance"));
        *m_genie_fractional_covariance_matrix += *((TMatrixT<double>*)f_syst->Get("individualDir/NormNCCOH_UBGenie_frac_covariance"));
        *m_genie_fractional_covariance_matrix += *((TMatrixT<double>*)f_syst->Get("individualDir/RPA_CCQE_UBGenie_frac_covariance"));
        *m_genie_fractional_covariance_matrix += *((TMatrixT<double>*)f_syst->Get("individualDir/Theta_Delta2Npi_UBGenie_frac_covariance"));
        *m_genie_fractional_covariance_matrix += *((TMatrixT<double>*)f_syst->Get("individualDir/VecFFCCQEshape_UBGenie_frac_covariance"));
        *m_genie_fractional_covariance_matrix += *((TMatrixT<double>*)f_syst->Get("individualDir/XSecShape_CCMEC_UBGenie_frac_covariance"));


	this->RemoveNan(m_genie_fractional_covariance_matrix);
/*	for(int i=0; i<m_genie_fractional_covariance_matrix->GetNcols(); i++){
            for(int j=0; j<m_genie_fractional_covariance_matrix->GetNcols(); j++){
                if(std::isnan((*m_genie_fractional_covariance_matrix)(i,j))) (*m_genie_fractional_covariance_matrix)(i,j) =0.0;
            }
        }
*/
	f_syst->Close();
	
	return 0;
}


int SBNsinglephoton::CalcFullButGenieFractionalCovarMatrix(){

	if(m_full_fractional_covariance_matrix==NULL || m_genie_fractional_covariance_matrix==NULL){
	   std::cout<< "SBNsinglephoton::CalcFullButGenieFractionalCovarMatrix\t|| Either full fractional covar matrix or genie fractional covariance matrix has NOT been setup yet."<< std::endl;
	   exit(EXIT_FAILURE); 
	}

	if(is_verbose) std::cout << "SBNsinglephoton::CalcFullButGenieFractionalCovarMatrix\t|| as the name says.." << std::endl;
	m_full_but_genie_fractional_covariance_matrix = new TMatrixT<double>(num_bins_total, num_bins_total);
	*m_full_but_genie_fractional_covariance_matrix  = (*m_full_fractional_covariance_matrix) - (*m_genie_fractional_covariance_matrix);

	return 0;
}


double SBNsinglephoton::ScaleFactor(double E, double factor, std::vector<double>& vec){
	double scale = factor;
	switch(vec.size()){
	   case 0:
		std::cout<< "SBNsinglephoton::ScaleFactor\t|| No energy/momentum dependent scaling applied" << std::endl;
		break;
	   case 1:
		std::cout << "SBNsinglephoton::ScaleFactor\t|| Applying Linear energy/momentum dependent scaling!"<< std::endl;
		scale += vec[0]*E;
		break;
	   case 2:
		std::cout << "SBNsinglephoton::ScaleFactor\t|| Applying 2nd order polynomial energy/momentum dependent scaling!" << std::endl;
		scale += vec[0]*E+vec[1]*E*E;
		break;
	   default: scale += 0;
	}

	return scale;

}

/*int SBNsinglephoton::WriteOutSpec(SBNspec* inspec, std::string tag){
	std::cout<<"SBNsinglephoton::WriteOut()\t\t||\t\tWriting out SBNspec with tag "<<tag<<std::endl
	inspec->WriteOut(tag);
	return 0;
}
*/

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

int SBNsinglephoton::SaveHistogram(){
        return this->SaveHistogram(m_map);
}
int SBNsinglephoton::SaveHistogram(std::map<int, std::vector<double>>& inmap){


    if(inmap.empty()){
	std::cout << "SBNsinglephoton::SaveHistogram\t|| map is empty!!" << std::endl;
	exit(EXIT_FAILURE);	
    }else{
        std::map<int, std::vector<double>> map_chi = inmap;
        std::map<int, std::vector<double>>::iterator itmap = map_chi.begin();
	int best_point = itmap->first;
	std::vector<double> vec_chi = itmap->second;
	double chi_min = *std::min_element(vec_chi.begin(), vec_chi.end());
	std::for_each(vec_chi.begin(), vec_chi.end(), [&chi_min](double& d){d -= chi_min;});  //get the delta_chi vector


	int interpolation_number = 1000;
	if(m_grid.f_num_dimensions == 2){
	   TFile* fout = new TFile("NCpi0_normalization_fit_output.root", "UPDATE");
	
	   if(!m_bool_poly_grid){
	    	if(is_verbose) std::cout<< "SBNsinglephoton::SaveHistogram\t|| Case: NCpi0 normalization fit, no energy/momentum dependent scaling!" << std::endl;
		auto grid_x = m_grid.f_dimensions.at(0);
		auto grid_y = m_grid.f_dimensions.at(1);

		TH2D* h_chi_surface = this->Do2DInterpolation(interpolation_number, grid_y.f_points, grid_x.f_points, vec_chi);
		h_chi_surface->SetName("h_chi_interpolated_surface");
		h_chi_surface->SetTitle(Form("#Delta#chi^{2} surface; %s;%s",grid_y.f_name.c_str(), grid_x.f_name.c_str()));
		h_chi_surface->Write();
	   }else{
		if(is_verbose) std::cout<< "SBNsinglephoton::SaveHistogram\t|| Case: NCpi0 normalization fit, with energy/momentum dependent scaling to " << m_poly_grid.f_num_dimensions << "nd order!" << std::endl;
		
	   }
	   fout->Close();
	} //end of 2 dimension case
	else if(m_grid.f_num_dimensions == 1){
	   TFile* fout = new TFile("NCdelta_fit_output.root", "UPDATE");
	   TH1D* h_dchi = new TH1D("h_delta_chi", Form("#Delta#chi^{2} distribution;%s; #Delta#chi^{2} ",m_grid.f_dimensions.at(0).f_name.c_str()), m_grid.f_num_total_points, m_grid.f_dimensions.at(0).f_min, m_grid.f_dimensions.at(0).f_max);
	   for(int i=0 ;i< vec_chi.size(); i++){
		std::vector<double> ipoint = m_vec_grid[i];
		//h_dchi->Fill(ipoint[0], vec_chi[i]);
		h_dchi->SetBinContent(i+1, vec_chi[i]);
	   }

	   h_dchi->Write();
	   fout->Close();
	}//end of 1 dimension case
	else{

	   TFile* fout = new TFile("NCdelta_fit_output.root", "UPDATE");
	   if(!m_bool_poly_grid){
		if(is_verbose) std::cout<< "SBNsinglephoton::SaveHistogram\t|| Case: NC delta and pi0 combined fit, no energy/momentum dependent scaling!" << std::endl;
		auto grid_x = m_grid.f_dimensions.at(0);
		auto grid_y = m_grid.f_dimensions.at(1);
		auto grid_z = m_grid.f_dimensions.at(2);
		std::vector<double> temp_best_point = m_vec_grid[best_point];

		//marginalize over 1 parameter	
		TH2D* h_mchi2_xy = new TH2D("h_mchi2_xy", Form("h_mchi2_xy; %s;%s", grid_x.f_name.c_str(), grid_y.f_name.c_str()), grid_x.f_N, grid_x.f_min, grid_x.f_max, grid_y.f_N, grid_y.f_min, grid_y.f_max);
		TH2D* h_mchi2_yz = new TH2D("h_mchi2_yz", Form("h_mchi2_yz; %s;%s", grid_y.f_name.c_str(), grid_z.f_name.c_str()), grid_y.f_N, grid_y.f_min, grid_y.f_max, grid_z.f_N, grid_z.f_min, grid_z.f_max);
		TH2D* h_mchi2_xz = new TH2D("h_mchi2_xz", Form("h_mchi2_xz; %s;%s", grid_x.f_name.c_str(), grid_z.f_name.c_str()), grid_x.f_N, grid_x.f_min, grid_x.f_max, grid_z.f_N, grid_z.f_min, grid_z.f_max);
		//global minimum
		TH2D* h_gchi2_xy = new TH2D("h_gchi2_xy", Form("h_gchi2_xy; %s;%s", grid_x.f_name.c_str(), grid_y.f_name.c_str()), grid_x.f_N, grid_x.f_min, grid_x.f_max, grid_y.f_N, grid_y.f_min, grid_y.f_max);
		TH2D* h_gchi2_yz = new TH2D("h_gchi2_yz", Form("h_gchi2_yz; %s;%s", grid_y.f_name.c_str(), grid_z.f_name.c_str()), grid_y.f_N, grid_y.f_min, grid_y.f_max, grid_z.f_N, grid_z.f_min, grid_z.f_max);
		TH2D* h_gchi2_xz = new TH2D("h_gchi2_xz", Form("h_gchi2_xz; %s;%s", grid_x.f_name.c_str(), grid_z.f_name.c_str()), grid_x.f_N, grid_x.f_min, grid_x.f_max, grid_z.f_N, grid_z.f_min, grid_z.f_max);

		//minimize over two parameters
		TH1D* h_chi_delta = new TH1D("h_chi_delta", Form("h_chi_delta; %s;#Delta#chi^{2}",grid_z.f_name.c_str()), grid_z.f_N, grid_z.f_min, grid_z.f_max);

		for(int ix=1;ix <= grid_x.f_N; ix++){
		        for(int iy=1; iy <= grid_y.f_N; iy++) h_mchi2_xy->SetBinContent(ix, iy, DBL_MAX);
        		for(int iz=1; iz <= grid_z.f_N; iz++) h_mchi2_xz->SetBinContent(ix, iz, DBL_MAX);
   		}

   		for(int iz=1; iz <= grid_z.f_N; iz++){
        		for(int iy=1; iy <= grid_y.f_N; iy++)  h_mchi2_yz->SetBinContent(iy, iz, DBL_MAX);
        		h_chi_delta->SetBinContent(iz, DBL_MAX);
   		}


		for(int ix=0; ix < grid_x.f_N; ix++){
       		    for(int iy=0; iy < grid_y.f_N; iy++){
           		for(int iz=0 ; iz< grid_z.f_N; iz++){
                	    int ip = ix*grid_y.f_N*grid_z.f_N + iy*grid_z.f_N + iz; // index of grid point
		            std::vector<double> point = m_vec_grid[ip];


                            //marginalized minimum
                            //conditional operator, saver the smaller chi.
                	    if(vec_chi[ip]< h_mchi2_xy->GetBinContent(ix+1, iy+1)){
                        	 h_mchi2_xy->SetBinContent(ix+1, iy+1, vec_chi[ip]);
                         	//std::cout << "chi2 value: " << chi[ip] << std::endl;

		 	    }
                	    if(vec_chi[ip]< h_mchi2_xz->GetBinContent(ix+1, iz+1)) h_mchi2_xz->SetBinContent(ix+1, iz+1, vec_chi[ip]);
                	    if(vec_chi[ip]< h_mchi2_yz->GetBinContent(iy+1, iz+1)) h_mchi2_yz->SetBinContent(iy+1, iz+1, vec_chi[ip]);


               		    //global minimum
                	    if(point[2] == temp_best_point[2]) h_gchi2_xy->Fill(point[0], point[1], vec_chi[ip]);
                	    if(point[1] == temp_best_point[1]) h_gchi2_xz->Fill(point[0], point[2], vec_chi[ip]);
               		    if(point[0] == temp_best_point[0]) h_gchi2_yz->Fill(point[1], point[2], vec_chi[ip]);

                 	    //marginalize two parameters
                	    if(vec_chi[ip] < h_chi_delta->GetBinContent(iz+1)) h_chi_delta->SetBinContent(iz+1, vec_chi[ip]);
           		}
        	    }
   		}
		
		h_mchi2_xy->Write(); h_mchi2_xz->Write(); h_mchi2_yz->Write();
		h_gchi2_xy->Write(); h_gchi2_xz->Write(); h_gchi2_yz->Write();
		h_chi_delta->Write();

	   } //end of the m_bool_poly_grid loop
	   fout->Close();
	}//end of 3 dimension case
    } //end of map check

    return 0;
}


//2D interpolation
//x, y vector should only be in increasing order
//NOTE::elements of the vector should first change with vector 'x', then change with vector 'y'!! Order is important here!
TH2D* SBNsinglephoton::Do2DInterpolation(int inter_number, std::vector<double>& x, std::vector<double>& y, std::vector<double>& value){
	double x_min = x[0];
	double x_max = x[x.size()-1];
	double y_min = y[0];
	double y_max = y[y.size() -1];
	//interpolation step size
	double x_step = fabs(x_max - x_min)/double(inter_number -1);
	double y_step = fabs(y_max - y_min)/double(inter_number -1);

	TH2D* h_inter = new TH2D("h_inter", "h_inter", inter_number, x_min, x_max, inter_number, y_min, y_max);
	
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
		  h_inter->Fill(xi, yj, zij);
	      }
	}

	//freee up the pointers
	gsl_spline2d_free(spline);
        gsl_interp_accel_free(xacc);
        gsl_interp_accel_free(yacc);

	return h_inter;	
}


int SBNsinglephoton::PrintOutFitInfo(std::map<int, std::vector<double>>& inmap, std::string intag, bool print_all){
	std::map<int, std::vector<double>> map_chi = inmap;
        std::map<int, std::vector<double>>::iterator itmap = map_chi.begin();

        int best_point = itmap->first;
        std::vector<double> vec_chi = itmap->second;
        double chi_min = *std::min_element(vec_chi.begin(), vec_chi.end());

	//best-fit point index at two grid
	int m_poly_grid_index = best_point/m_num_total_gridpoints;
	int m_grid_index = best_point%m_num_total_gridpoints;

	//print out info
	std::vector<double> m_point = m_vec_grid[m_grid_index];
	std::cout << intag<<": Best chi value is " << chi_min << " at flat point";
	for(int i=0; i<m_point.size(); i++)
	   std::cout << ", " << m_grid.f_dimensions[i].f_name << ": "<< m_point[i] ;
	if(m_bool_poly_grid){
	   std::vector<double> m_poly_point = m_vec_poly_grid[m_poly_grid_index];
	   
	   std::cout << ", and energy/momentum dependet shift for NCpi0 non-coherent component with";
	   for(int i=0 ; i< m_poly_point.size();i++)
		std::cout << ", "<< i+1 <<"nd order factor "<< m_poly_point[i];
	}
	std::cout << ". "<< std::endl;

	if(print_all){
	   std::cout << intag<< "========================Below are detailed dchi and coordinate values==========================" << std::endl;
	   for(int i=0;i<vec_chi.size(); i++){
		int temp_poly_grid_index = i/m_num_total_gridpoints;
		int temp_grid_index = i%m_num_total_gridpoints;
		std::cout << "Chi:" << vec_chi[i] << " Point: ";
		std::vector<double> temp_point = m_vec_grid[temp_grid_index];
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


int SBNsinglephoton::SetFlatFullFracCovarianceMatrix(double flat){

	m_full_fractional_covariance_matrix = new TMatrixT<double>(num_bins_total, num_bins_total);
	for(int i=0;i<num_bins_total; i++) (*m_full_fractional_covariance_matrix)(i,i) = flat*flat;

	return 0;	
}

int SBNsinglephoton::SetStatOnly(){
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
	for(int i=0 ; i<m_grid.f_num_dimensions; i++){
	    if(m_grid.f_dimensions[i].f_has_constrain){
		temp_constrain_param.push_back(m_grid.f_dimensions[i].f_constrain_value);
	    }
	}

	if(m_bool_poly_grid){

           for(int i=0; i< m_poly_grid.f_num_dimensions; i++){
              if(m_poly_grid.f_dimensions[i].f_has_constrain){
		     temp_constrain_param.push_back(m_poly_grid.f_dimensions[i].f_constrain_value);
              }
	   }
	}
	
	return this->ModifyCV(infactor, temp_constrain_param);
}

int SBNsinglephoton::ModifyCV(double infactor, std::vector<double> param){


	if(!m_bool_cv_spectrum_loaded){
                this->LoadCV();
        }


	std::cout<< "SBNsinglephoton::ModifyCV\t|| Start modifying CV spectrum" <<std::endl;
	int index = 0;
	//apply flat normalization
	for(int i=0 ; i<m_grid.f_num_dimensions; i++){
            if(m_grid.f_dimensions[i].f_has_constrain  && (index < param.size()) ){
                m_cv_spectrum->Scale(m_grid.f_dimensions[i].f_name, param[index]);
		index++;
            }
        }

	//energy/momentum dependent scaling
	if(m_bool_poly_grid){
	   std::ostringstream prescale_tag;
	   std::vector<double> temp_scale_parameter;
           for(int i=0; i< m_poly_grid.f_num_dimensions; i++){
              if(m_poly_grid.f_dimensions[i].f_has_constrain && (index < param.size())  ){   
                     prescale_tag << "_" << std::fixed<< std::setprecision(3) << param[index];
		     temp_scale_parameter.push_back(param[index]);
		     index++;
              }
           }

	   if(temp_scale_parameter.size() != 0){
	       std::string temp_filename = tag+"_PreScaled"+prescale_tag.str()+".SBNspec.root";
	       //check if a file exits
	       if(gSystem->AccessPathName(temp_filename.c_str())){
		   std::cout << "SBNsinglephoton::ModifyCV\t|| Prescaled file doesn't exist, start generating it..." << std::endl;
		   this->OpenFiles();
		   this->PreScaleSpectrum(xmlname, temp_scale_parameter);
		   this->CloseFiles();
	       }

	       std::cout <<"SBNsinglephoton::ModifyCV\t|| Adding pre-scaled spectrum to the CV." << std::endl;
	       SBNspec temp_prescale(temp_filename.c_str(), xmlname, false);
	       m_cv_spectrum->Add(&temp_prescale); 
	   }
	}

	if(m_bool_data_spectrum_loaded) m_cv_spectrum->Scale("NCDeltaRadOverlayLEE", (infactor-1)*0.5);
        m_bool_modify_cv = true;
	m_cv_delta_scaling = infactor;

	return 0;
}

SBNspec SBNsinglephoton::GeneratePointSpectra(int np){
	int m_poly_grid_index = np/m_num_total_gridpoints;
        int m_grid_index = np%m_num_total_gridpoints;

	SBNspec spec_rgrid_part((tag+"_CV.SBNspec.root").c_str(), xmlname, false);

	spec_rgrid_part.Scale("BNBext", 0.0);
	spec_rgrid_part.Scale("Dirt", 0.0);

	std::vector<double> temp_p_grid = m_vec_grid[m_grid_index];
	for(int i=0;i<m_grid.f_num_dimensions; i++){
	    if(m_grid.f_dimensions[i].f_name == "NCDeltaRadOverlayLEE" && temp_p_grid[i]<0 ){
		spec_rgrid_part.Scale("NCDeltaRadOverlayLEE", 0.0);
		spec_rgrid_part.Scale("NCDeltaRadOverlaySM", 1+temp_p_grid[i]*2.0);
	    }
	    else
		spec_rgrid_part.Scale(m_grid.f_dimensions[i].f_name, temp_p_grid[i]);
	}

	if(m_bool_poly_grid){
	    double non_coh_factor;   
	    for(int i=0;i<m_grid.f_num_dimensions; i++){
		if(m_grid.f_dimensions[i].f_name == "NCPi0NotCoh"){
		    non_coh_factor=temp_p_grid[i];
		    spec_rgrid_part.Scale("NCPi0NotCoh", 0.0);
		}
	    }

	    std::vector<double> temp_p_polygrid=m_vec_poly_grid[m_poly_grid_index];
	    std::ostringstream prescale_tag;
    //generate tag for prescaled spectra
	    prescale_tag<< "Flat_"<< std::fixed<<std::setprecision(3) <<non_coh_factor << "_PreScaled";
            for(int i=0; i< temp_p_polygrid.size(); i++){
                prescale_tag << "_" << std::fixed<< std::setprecision(3) << temp_p_polygrid[i];
            }

	    this->OpenFiles();
	    m_bool_cv_spectrum_generated=true;
	    this->PreScaleSpectrum(xmlname, non_coh_factor, temp_p_polygrid);
	    this->CloseFiles();

	    SBNspec spec_polygrid_part((tag+"_"+prescale_tag.str()+".SBNspec.root").c_str(), xmlname, false); 
	
	    spec_rgrid_part.Add(&spec_polygrid_part);
	}


	return spec_rgrid_part;
}
