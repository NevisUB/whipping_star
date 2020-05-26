#include "SBNfeld.h"
using namespace sbn;


NeutrinoModel SBNfeld::convert3p1(std::vector<double> ingrid){
    if(ingrid.size()!=3){
        std::cout<<"ERROR: Currently must have a 3-grid, for Dm ue4 and um4. grid.size(): "<<ingrid.size()<<std::endl;
        exit(EXIT_FAILURE);
    }

    NeutrinoModel signalModel(pow(10,ingrid[0]),pow(10,ingrid[1]),pow(10,ingrid[2]));

    return signalModel;
}

int SBNfeld::GenerateOscillatedSpectra(){

    //This will loop over the predefined grid for oscillation studies, and calculate it

    int n_mass = m_grid.f_dimensions[0].GetNPoints();
    bool cv_runone = true;
    std::cout<<"Beginining to Generate all oscillated spectra for : "<<n_mass<<" Mass splittings."<<std::endl;

    for(size_t t =0; t< n_mass; t++){

        std::cout<<"SBNfeld::GenerateOcillatedSpectra()\t\t||\t\t On Mass Splitting number "<<t<<" which is "<<m_grid.f_dimensions[0].GetPoint(t)<<std::endl;
        //need to convert from this gridpoints to a neutrinoModel (Blarg, don't like this step but needed at the moment)
        NeutrinoModel this_model(pow(10,m_grid.f_dimensions[0].GetPoint(t)),1,1); // = this->convert3p1(m_vec_grid.at(t)); 
        this_model.Printall();

        // on construction it makes 3 SBNspecs, 1 sin amp, 1 sin2 amp, 1 CV oscilatted
        SBNgenerate gen(this->xmlname,this_model);

        // Write them to file
        gen.WritePrecomputedOscSpecs(this->tag);

        if(cv_runone){
            gen.spec_central_value.WriteOut(this->tag+"_CV");
            cv_runone= false;
        }
    }
}


int SBNfeld::GenerateScaledSpectra(){

    //Rather than oscillations, we just want to pick 1 subchannel, and scale it up and down. The channel to scale is defined in xml, scale_signal="true"
    // Basically, we take the core spectra, scale it by each gridpoint, and push that back

    std::cout<<"Beginining to Generate all Scaled spectra."<<std::endl;
    m_cv_spec_grid.clear();
    m_cv_spec_grid.resize(m_grid.f_num_total_points);
    m_core_spectrum->CalcFullVector();


    for(size_t t =0; t < m_grid.f_num_total_points; t++){

        std::cout<<"SBNfeld::GenerateScaledSpectra()\t\t||\t\t On scaling point "<<t<<" which is "<<m_grid.f_dimensions[0].GetPoint(t)<<std::endl;

        m_cv_spec_grid[t] = new SBNspec(m_core_spectrum->full_vector,  m_core_spectrum->xmlname, t, false);
        m_cv_spec_grid[t]->Scale(m_subchannel_to_scale, m_grid.f_dimensions[0].GetPoint(t));

        m_cv_spec_grid[t]->CalcFullVector();
        m_cv_spec_grid[t]->CollapseVector();
    }
}


int SBNfeld::GenerateBackgroundSpectrum(){

    std::cout<<"SBNfeld::GenerateBackgroundSpectrum()\t\t||\t\t Generating a background 3+0 spectra"<<std::endl;
    NeutrinoModel this_model = this->convert3p1(m_vec_grid[0]); 
    NeutrinoModel background_only_model(this_model.mNu[0],0.0,0.0); // quirk, this works better

    m_core_spectrum->LoadModel(background_only_model);

    //Is this a troublesome line?!? Shouldn't be right?!
    //m_core_spectrum->SetAppMode();

    std::vector<double> ans = m_core_spectrum->Oscillate(this->tag, false);
    SBNspec background(ans, m_core_spectrum->xmlname,false);
    background.Scale("fullosc",0.0);
    background.WriteOut(this->tag+"_BKG_ONLY");

    std::cout<<"SBNfeld::GenerateBackgroundSpectrum()\t\t||\t\t Done."<<std::endl;
    return 0;

}

int SBNfeld::SetCoreSpectrum(std::string file){

    m_core_spectrum= new SBNosc(file,this->xmlname);
    m_bool_core_spectrum_set = true;
    return 0;
}

int SBNfeld::AddFlatDetSystematic(double percent){

    std::cout<<"Adding a "<<percent<<" percent to the diagonal "<<percent*percent<<std::endl;
    for(int j=0; j< m_full_fractional_covariance_matrix->GetNrows(); j++){
        std::cout<<j<<" before "<<(*m_full_fractional_covariance_matrix)(j,j);
        (*m_full_fractional_covariance_matrix)(j,j)+=percent*percent;   
        std::cout<<" after "<<(*m_full_fractional_covariance_matrix)(j,j)<<std::endl;

    }

    return 0 ;
}

int SBNfeld::SetFractionalCovarianceMatrix(std::string filename,std::string matrixname){
    TFile * fsys = new TFile(filename.c_str(),"read");
    m_full_fractional_covariance_matrix = (TMatrixD*)fsys->Get(matrixname.c_str());
    fsys->Close();

    return 0;
}

int SBNfeld::SetFractionalCovarianceMatrix(TMatrixT<double> * in){

    m_full_fractional_covariance_matrix = in;
    return 0;
}

int SBNfeld::SetEmptyFractionalCovarianceMatrix(){

    m_full_fractional_covariance_matrix = new TMatrixT<double>(num_bins_total,num_bins_total);
    m_full_fractional_covariance_matrix->Zero();
    return 0;
}

int SBNfeld::SetRandomSeed(double sin){
    m_random_seed = sin;
    std::cout<<"SBNfeld::SetRandomSeed()\t\t||\t\tSetting random seed to: "<<sin<<std::endl;
    return 0;
}

int SBNfeld::LoadPreOscillatedSpectrum(int which_pt){

    std::cout<<"Loading Just Grid Pt: "<<which_pt<<" ";
    for(int k=0; k<m_vec_grid[which_pt].size();k++){
        std::cout<<" "<<m_vec_grid[which_pt][k];
    }
    std::cout<<std::endl;

    //need to convert from this gridpoints to a neutrinoModel (Blarg, don't like this step but needed at the moment)
    NeutrinoModel this_model = this->convert3p1(m_vec_grid[which_pt]); 

    // this_model.Printall();
    //And load thus model into our spectra. At this point its comuted all the necessary mass-splittins and which frequencies they are
    m_core_spectrum->LoadModel(this_model);
    //m_core_spectrum->SetAppMode();

    //And apply this oscillaion! Adding to it the bkgSpec that it was initilised with.
    //NOTE we want to return the FULL spectrum, not compressed so we can calculate the covariance matrix, hense the false in this Oscilate
    std::vector<double> ans = m_core_spectrum->Oscillate(this->tag, false);
    std::cout<<"Spectrum: ";
    for(int p=0; p<ans.size();p++){
        std::cout<<" "<<ans[p];
    }
    SBNspec * thispoint = new SBNspec(ans, m_core_spectrum->xmlname,which_pt, false);
    
    thispoint->ScaleAll(global_scale);
    thispoint->CollapseVector();

    //make a print out of this exact spectrum as compared to the "core" spectrum
    std::string tlog  = std::to_string(m_vec_grid[which_pt][0])+"_"+std::to_string(m_vec_grid[which_pt][1])+"_"+std::to_string(m_vec_grid[which_pt][2]);
    m_core_spectrum->CompareSBNspecs(thispoint,tlog); 

    return 0;
}

int SBNfeld::LoadPreOscillatedSpectra(){
    //This is where we load the precalculated spectra and assemble into a actuall oscillate spectrum
    //we have a core spectrum loaded on the start of the SBNfeld class. This is the fullosc+intrinsics...etc..
    // i.e m_cv_spec_grid;

    m_cv_spec_grid.clear();
    std::cout<<"Beginining to Grab all oscillated spectra"<<std::endl;
    m_cv_spec_grid.resize(m_num_total_gridpoints);

    for(size_t t =0; t < m_num_total_gridpoints; t++){
        std::cout<<"SBNfeld::LoadPreOcillatedSpectra()\t\t||\t\t On spectrum "<<t<<"/"<<m_num_total_gridpoints;
        for(int k=0; k<m_vec_grid[t].size();k++){
            std::cout<<" "<<m_vec_grid[t][k];
        }
        std::cout<<std::endl;

        //need to convert from this gridpoints to a neutrinoModel (Blarg, don't like this step but needed at the moment)
        NeutrinoModel this_model = this->convert3p1(m_vec_grid[t]); 

        // this_model.Printall();
        //And load thus model into our spectra. At this point its comuted all the necessary mass-splittins and which frequencies they are
        m_core_spectrum->LoadModel(this_model);
        //m_core_spectrum->SetAppMode();

        //And apply this oscillaion! Adding to it the bkgSpec that it was initilised with.
        //NOTE we want to return the FULL spectrum, not compressed so we can calculate the covariance matrix, hense the false in this Oscilate
        std::vector<double> ans = m_core_spectrum->Oscillate(this->tag, false);
        std::cout<<"Spectrum: ";
        for(int p=0; p<ans.size();p++){
            std::cout<<" "<<ans[p];
        }
        std::cout<<std::endl;
        if(m_vec_grid[t][0] == 0.0603 && m_vec_grid[t][2] == -0.869) m_cv_spec_grid[t] = new SBNspec(ans, m_core_spectrum->xmlname,-11, false);
        //if(m_vec_grid[t][0] == 0.04 && m_vec_grid[t][1] == -1.562) m_cv_spec_grid[t] = new SBNspec(ans, m_core_spectrum->xmlname,-11, false);
	else m_cv_spec_grid[t] = new SBNspec(ans, m_core_spectrum->xmlname, t, false);
        m_cv_spec_grid[t]->ScaleAll(global_scale);
        m_cv_spec_grid[t]->CollapseVector();
	//added by guanqun to write out oscillated spectra
	if(m_vec_grid[t][0] == 0.0603 && m_vec_grid[t][2] == -0.869) m_cv_spec_grid[t]->WriteOut("numu_data_spectrum");	
	//if(m_vec_grid[t][0] == 0.04 && m_vec_grid[t][1] == -1.562) m_cv_spec_grid[t]->WriteOut("test_oscillated");	
	// end of the code added by guanqun.

        if(m_bool_print_comparasons && t ==490){// t==1668 
            //make a print out of this exact spectrum as compared to the "core" spectrum
            std::string tlog  = std::to_string(m_vec_grid[t][0])+"_"+std::to_string(m_vec_grid[t][1])+"_"+std::to_string(m_vec_grid[t][2]);
            m_core_spectrum->CompareSBNspecs(m_cv_spec_grid[t],tlog); 
        }
    }

    return 0;
}




int SBNfeld::LoadBackgroundSpectrum(){
    m_background_spectrum= new SBNosc(this->tag+"_BKG_ONLY.SBNspec.root",this->xmlname);
    m_background_spectrum->ScaleAll(global_scale);
    m_bool_background_spectrum_set = true;

    m_background_spectrum->CollapseVector();
    m_tvec_background_spectrum = new TVectorT<double>(m_background_spectrum->full_vector.size(), &(m_background_spectrum->full_vector)[0]);

    if(m_background_spectrum->full_vector.size() !=  m_full_fractional_covariance_matrix->GetNcols()){

        std::cout<<"SBNfeld::LoadBackgroundSpectrum || ERROR!! background spectrum is of length : "<<m_background_spectrum->full_vector.size()<<" and frac matrix is of size "<<m_full_fractional_covariance_matrix->GetNcols()<<std::endl;
        exit(EXIT_FAILURE);
    }


    m_background_chi = new SBNchi(*m_background_spectrum, *m_full_fractional_covariance_matrix, this->xmlname, false);
    return 0;
}

int SBNfeld::LoadBackgroundSpectrum(std::string filein){
    m_background_spectrum= new SBNosc(filein,this->xmlname);
    m_bool_background_spectrum_set = true;

    m_background_spectrum->CollapseVector();
    m_tvec_background_spectrum = new TVectorT<double>(m_background_spectrum->full_vector.size(), &(m_background_spectrum->full_vector)[0]);
    m_background_chi = new SBNchi(*m_background_spectrum, *m_full_fractional_covariance_matrix, this->xmlname, false);
    return 0;
}


int SBNfeld::CalcSBNchis(){
    //This is where we will calculate a vector of SBNchi's
    //i.e m_sbnchi_grid;

    TMatrixT<double> stat_only_matrix(num_bins_total,num_bins_total);
    stat_only_matrix.Zero();

    for(size_t t =0; t < m_num_total_gridpoints; t++){
        //Setup a SBNchi for this true point
        std::cout<<"Setting up grid point SBNchi "<<t<<std::endl;

        if(m_bool_stat_only){
            m_sbnchi_grid.push_back(new SBNchi(*m_cv_spec_grid.at(t), stat_only_matrix, this->xmlname, false, m_random_seed)) ;

        }else{
            m_sbnchi_grid.push_back(new SBNchi(*m_cv_spec_grid.at(t), *m_full_fractional_covariance_matrix, this->xmlname, false, m_random_seed)) ;
        }


    }

    return 0;
}


int SBNfeld::FullFeldmanCousins(){
    int num_universes = m_num_universes;

    //Ok take the background only spectrum and form a background only covariance matrix. CalcCovarianceMatrix includes stats
    TMatrixT<double> background_full_covariance_matrix = m_sbnchi_grid[0]->CalcCovarianceMatrix(m_full_fractional_covariance_matrix, *m_tvec_background_spectrum);
    TMatrixT<double> background_collapsed_covariance_matrix(m_background_spectrum->num_bins_total_compressed, m_background_spectrum->num_bins_total_compressed);
    m_sbnchi_grid[0]->CollapseModes(background_full_covariance_matrix, background_collapsed_covariance_matrix);    
    TMatrixT<double> inverse_background_collapsed_covariance_matrix = m_sbnchi_grid[0]->InvertMatrix(background_collapsed_covariance_matrix);   

    TFile *fout =  new TFile(("SBNfeld_output_"+tag+".root").c_str(),"recreate");
    fout->cd();


    for(size_t t =0; t < m_num_total_gridpoints; t++){

        SBNspec * true_spec = m_cv_spec_grid.at(t);
        SBNchi  * true_chi = m_sbnchi_grid.at(t);          

        std::cout<<"Starting on point "<<t<<"/"<<m_num_total_gridpoints<<std::endl;
        double sum = std::accumulate(true_spec->full_vector.begin(),true_spec->full_vector.end(),0.0);
        std::cout<<"Sum of events "<<sum<<std::endl;
        std::vector<double> vec_delta_chi(num_universes,0);
        std::vector<double> vec_chi_min(num_universes,0);
        double delta_chi_critical = DBL_MIN;
        //Some output informatation
        double tree_delta_chi = 0;
        double tree_chi_min = 0;
        int tree_bf_grid=0;

        fout->cd();
        TTree t_outtree(("ttree_"+std::to_string(t)).c_str(),("ttree_"+std::to_string(t)).c_str());
        t_outtree.Branch("delta_chi2",&tree_delta_chi);
        t_outtree.Branch("chi2_min",&tree_chi_min);
        t_outtree.Branch("bf_gridpoint",&tree_bf_grid);       


        for(size_t i=0; i< num_universes; i++){

            //step 0. Make a fake-data-experimet for this point, drawn from covariance
            std::vector<float>  fake_data= true_chi->SampleCovariance(true_spec);
         
            //Whats the Best Fit Point 
            std::vector<double> ans = this->PerformIterativeFit(fake_data,t,inverse_background_collapsed_covariance_matrix);

            //step 4 calculate the delta_chi for this universe
            vec_delta_chi[i] = ans[1]-ans[2];
            vec_chi_min[i] = ans[2];

            tree_delta_chi = vec_delta_chi[i];
            tree_chi_min = vec_chi_min[i];
            tree_bf_grid = (int)ans[0]; 

            t_outtree.Fill();
        }

        double tmin = DBL_MAX;
        double tmax = DBL_MIN;
        for(double&v:vec_delta_chi){
            tmin=std::min(tmin,v);
            tmax=std::max(tmax,v);
        }

        std::cout<<"For this point, minimum delta chi is "<<tmin<<" max is "<<tmax<<std::endl;
        TH1D h_delta_chi(("dcuni_"+std::to_string(t)).c_str(),("dcuni_"+std::to_string(t)).c_str(),50,0.0,tmax*1.0);  // This will store all the delta_chi's from each universe for this g_true point
        for(double&v:vec_delta_chi) h_delta_chi.Fill(v);


        double cmin = DBL_MAX;
        double cmax = DBL_MIN;
        for(double&v:vec_chi_min){
            cmin=std::min(cmin,v);
            cmax=std::max(cmax,v);
        }

        std::cout<<"For this point, minimum chi is "<<cmin<<" max is "<<cmax<<std::endl;

        TH1D h_chi_min(("cmuni_"+std::to_string(t)).c_str(),("cmuni_"+std::to_string(t)).c_str(),50,0.0,cmax*1.0);  // This will store all the chi_mins from each universe for this g_true point
        for(double&v:vec_chi_min){
            h_chi_min.Fill(v);
        }

        /*
        //Now lets do a simple fit to a chi^2 
        std::string f_name = "fchi_"+std::to_string(t);
        //TF1 *fchi = new TF1(f_name.c_str(),[&](double*x, double *p){ return p[0]*gsl_ran_chisq_pdf(x[0],p[1]); },0,tmax,2); // gsl does not like this lambda
        TF1 *fchi = new TF1(f_name.c_str(),[&](double*x, double *p){ return p[0]*ROOT::Math::chisquared_pdf(x[0],p[1]); },0,tmax,2);
        fchi->SetParameter(0,1.0); 
        fchi->SetParameter(1,2.0); 
        fchi->SetParName(0,"norm");
        fchi->SetParName(1,"NDOF");
        fchi->SetParLimits(1,0,10.0);
        t_outtree.Fit(f_name.c_str(),"delta_chi2",("1.0/"+std::to_string((double)vec_delta_chi.size())).c_str(),"M");
        double fitted_ndof = fchi->GetParameter(1);
        std::cout<<"FIT NDOF: "<<fitted_ndof<<" "<<fchi->GetParameter(0)<<std::endl;
        */


        fout->cd();
        h_delta_chi.Write();
        h_chi_min.Write();
        t_outtree.Write();


        std::cout<<"This Gridpoint has a max and min Chi^2_min of ("<<cmin<<","<<cmax<<")"<<std::endl;
        std::cout<<"This Gridpoint has a max and min DeltaChi of ("<<tmin<<","<<tmax<<")"<<std::endl;

        //Step 3: We now have a distrubution of Delta_chi's based on many pseudo-fake data experiments. This should approximate a 
        //chi^2 distrubution with 2 DOF but there can be quite large variations. 

        double sig1 = 0.5-(0.6827)/2.0;
        double sig2 = 0.5-(0.9545)/2.0;
        std::vector<double> prob_values = {0.01,sig2,0.05,0.1,sig1,0.5,0.6,1-sig1,0.9,0.95,1-sig2,0.99};
        std::vector<double> delta_chi_quantiles(prob_values.size());	
        std::vector<double> chi_min_quantiles(prob_values.size());	

        h_delta_chi.ComputeIntegral(); 
        std::cout<<"DeltaChi Histogram has "<<h_delta_chi.GetEntries()<<" entries, with an integral of "<<h_delta_chi.Integral()<<std::endl;
        h_delta_chi.GetQuantiles(prob_values.size(), &delta_chi_quantiles[0], &prob_values[0]);


        h_chi_min.ComputeIntegral(); 
        std::cout<<"ChiMin Histogram has "<<h_chi_min.GetEntries()<<" entries, with an integral of "<<h_chi_min.Integral()<<std::endl;
        h_chi_min.GetQuantiles(prob_values.size(), &chi_min_quantiles[0], &prob_values[0]);


        std::cout<<"Grid Point: "<<t;
        for(auto &p: m_vec_grid.at(t)){
            std::cout<<" "<<p;
        }
        std::cout<<std::endl;

        for(int i=0; i< prob_values.size(); i++){
            std::cout<<"--delta_quantile "<<delta_chi_quantiles[i]<<" chimin_quantile "<<chi_min_quantiles[i]<<" @prob "<<prob_values[i]<<std::endl;
        }

        std::cout<<"Finished This Point"<<std::endl;
        //delta_chi_critical = delta_chi_quantiles[0];
        //break;
    }

    fout->Close();
    return 0;
};



std::vector<double> SBNfeld::PerformIterativeFit(std::vector<float> &datavec, size_t grid_pt, TMatrixT<double>& inverse_background_collapsed_covariance_matrix){

    double last_chi_min = DBL_MAX;
    int best_grid_point = -99;

    TMatrixT<double> inverse_current_collapsed_covariance_matrix = inverse_background_collapsed_covariance_matrix;  

    SBNspec * grid_spec = m_cv_spec_grid.at(grid_pt); 
    SBNchi  * grid_chi = m_sbnchi_grid.at(grid_pt); 

    for(size_t n_iter = 0; n_iter < m_max_number_iterations; n_iter++){

        //Step 1. What covariance matrix do we use?
        //For first iteration, use the precalculated background only inverse covariance matrix.

        //For all subsequent iterations what is the full covariance matrix? Use the last best grid point.
        if(n_iter!=0){

            //Calculate current full covariance matrix, collase it, then Invert. 
            TMatrixT<double> current_full_covariance_matrix = grid_chi->CalcCovarianceMatrix(m_full_fractional_covariance_matrix, m_cv_spec_grid[best_grid_point]->full_vector);
            TMatrixT<double> current_collapsed_covariance_matrix(num_bins_total_compressed, num_bins_total_compressed);
            grid_chi->CollapseModes(current_full_covariance_matrix, current_collapsed_covariance_matrix);    
            inverse_current_collapsed_covariance_matrix = grid_chi->InvertMatrix(current_collapsed_covariance_matrix);   
        }

        //Step 2.0 Find the global_minimum_for this universe. Integrate in SBNfit minimizer here, a grid scan for now.
        double chi_min = DBL_MAX;
        for(size_t r =0; r < m_num_total_gridpoints; r++){

            //double chi_tmp = this->CalcChi(datavec, m_cv_spec_grid[r]->collapsed_vector,  inverse_current_collapsed_covariance_matrix, true);
            double chi_tmp = this->CalcChi(datavec, m_cv_spec_grid[r]->collapsed_vector,  inverse_current_collapsed_covariance_matrix, false);

            if(chi_tmp < chi_min){
                best_grid_point = r;
                chi_min = chi_tmp;
            }
        }

        if(n_iter!=0){

            //std::cout<<"On iter: "<<n_iter<<" of uni "<<i<<"/"<<num_universes<<" w/ chi^2: "<<chi_min<<" lastchi^2: "<<last_chi_min<<" diff() "<<fabs(chi_min-last_chi_min)<<" tol: "<<chi_min_convergance_tolerance<<" best_grid_point: "<<best_grid_point<<std::endl;

            //Step 3.0 Check to see if mgrid_chi for this particular fake_data  has converged sufficiently
            if(fabs(chi_min-last_chi_min)< m_chi_min_convergance_tolerance){
                last_chi_min = chi_min;
                break;
            }
        }else{
            //std::cout<<"On iter: "<<n_iter<<" chi^2: "<<chi_min<<std::endl; 
        }

        last_chi_min = chi_min;
    }

    //Now use the curent_iteration_covariance matrix to also calc this_chi here for the delta.
    //double this_chi   = this->CalcChi(datavec, grid_spec->collapsed_vector,inverse_current_collapsed_covariance_matrix, true);
    double this_chi   = this->CalcChi(datavec, grid_spec->collapsed_vector,inverse_current_collapsed_covariance_matrix, false);


    //returns the BF grid, the chi^2 and the minimum_chi at the BF. 
    return {(double)best_grid_point,this_chi,last_chi_min};
}


int SBNfeld::PointFeldmanCousins(size_t grid_pt){

    int num_universes = m_num_universes;

    //Ok take the background only spectrum and form a background only covariance matrix. CalcCovarianceMatrix includes stats
    TMatrixT<double> background_full_covariance_matrix = m_sbnchi_grid[0]->CalcCovarianceMatrix(m_full_fractional_covariance_matrix, *m_tvec_background_spectrum);
    TMatrixT<double> background_collapsed_covariance_matrix(m_background_spectrum->num_bins_total_compressed, m_background_spectrum->num_bins_total_compressed);
    m_sbnchi_grid[0]->CollapseModes(background_full_covariance_matrix, background_collapsed_covariance_matrix);    
    TMatrixT<double> inverse_background_collapsed_covariance_matrix = m_sbnchi_grid[0]->InvertMatrix(background_collapsed_covariance_matrix);   

    if(grid_pt > m_num_total_gridpoints){
        std::cout<<"ERROR! SBNfeld::PointFeldmanCousins() || input grid_point ("<<grid_pt<<") is above num_total_gridpoints ("<<m_num_total_gridpoints<<"). Exiting."<<std::endl;
        exit(EXIT_FAILURE);
    }

    std::cout<<"Starting on grid point "<<grid_pt<<"/"<<m_num_total_gridpoints<<std::endl;

    SBNspec * true_spec = m_cv_spec_grid[grid_pt]; 
    SBNchi  * true_chi = m_sbnchi_grid[grid_pt]; 

    for(size_t i=0; i< num_universes; i++){

        //step 0. Make a fake-data-experimet for this point, drawn from covariance
        std::vector<float> fake_data= true_chi->SampleCovariance(true_spec);
        float last_chi_min = FLT_MAX;
        int best_grid_point = -99;

        TMatrixT<double> inverse_current_collapsed_covariance_matrix = inverse_background_collapsed_covariance_matrix;  

        size_t n_iter = 0;
        for(n_iter = 0; n_iter < m_max_number_iterations; n_iter++){

            //Step 1. What covariance matrix do we use?
            //For first iteration, use the precalculated background only inverse covariance matrix.
            //For all subsequent iterations what is the full covariance matrix? Use the last best grid point.
            if(n_iter!=0){
                //Calculate current full covariance matrix, collase it, then Invert. 
                TMatrixT<double> current_full_covariance_matrix = true_chi->CalcCovarianceMatrix(m_full_fractional_covariance_matrix,m_cv_spec_grid[best_grid_point]->full_vector);
                TMatrixT<double> current_collapsed_covariance_matrix(num_bins_total_compressed, num_bins_total_compressed);
                true_chi->CollapseModes(current_full_covariance_matrix, current_collapsed_covariance_matrix);    
                inverse_current_collapsed_covariance_matrix = true_chi->InvertMatrix(current_collapsed_covariance_matrix);   
            }

            //Step 2.0 Find the global_minimum_for this universe. Integrate in SBNfit minimizer here, a grid scan for now.
            float chi_min = FLT_MAX;

            for(size_t r =0; r < m_num_total_gridpoints; r++){
                float chi_tmp = this->CalcChi(fake_data, m_cv_spec_grid[r]->collapsed_vector,  inverse_current_collapsed_covariance_matrix);
                if(chi_tmp < chi_min){
                    best_grid_point = r;
                    chi_min = chi_tmp;
                }
            }

            if(n_iter!=0){
                //std::cout<<"On iter: "<<n_iter<<" of uni "<<i<<"/"<<num_universes<<" w/ chi^2: "<<chi_min<<" lastchi^2: "<<last_chi_min<<" diff() "<<fabs(chi_min-last_chi_min)<<" tol: "<<m_chi_min_convergance_tolerance<<" best_grid_point: "<<best_grid_point<<std::endl;

                //Step 3.0 Check to see if min_chi for this particular fake_data  has converged sufficiently
                if(fabs(chi_min-last_chi_min)< m_chi_min_convergance_tolerance){
                    last_chi_min = chi_min;
                    break;
                }
            }else{
                //std::cout<<"On iter: "<<n_iter<<" chi^2: "<<chi_min<<std::endl;
            }

            last_chi_min = chi_min;
        }

        //Now use the curent_iteration_covariance matrix to also calc this_chi here for the delta.
        float this_chi = this->CalcChi(fake_data, true_spec->collapsed_vector,inverse_current_collapsed_covariance_matrix);

        //step 4 calculate the delta_chi for this universe

        std::cout<<grid_pt<<" "<<i<<" "<<last_chi_min<<" "<<this_chi-last_chi_min<<" "<<best_grid_point<<" "<<n_iter<<std::endl;
    }

    std::cout<<"Finished Grid Point: "<<grid_pt;
    for(auto &p: m_vec_grid[grid_pt]){
        std::cout<<" "<<p;
    }
    std::cout<<std::endl;


    return 0;
};






float SBNfeld::CalcChi(std::vector<float>& data, std::vector<double>& prediction, TMatrixT<double> & inverse_covariance_matrix){
	return this->CalcChi(data, prediction, inverse_covariance_matrix, false);
};



float SBNfeld::CalcChi(std::vector<float>& data, std::vector<double>& prediction, TMatrixT<double> & inverse_covariance_matrix, bool det_term ){
    float tchi = 0;

    for(int i =0; i<num_bins_total_compressed; i++){
        for(int j =0; j<num_bins_total_compressed; j++){
            tchi += (data[i]-prediction[i])*inverse_covariance_matrix(i,j)*(data[j]-prediction[j] );
        }
    }
	
    if(det_term == true) tchi += log(inverse_covariance_matrix.Determinant());
   
    return tchi;
};


int SBNfeld::GlobalScan(){
    return this->GlobalScan(-1);
}

int SBNfeld::GlobalScan(int which_pt){

    //Ok take the background only spectrum and form a background only covariance matrix. CalcCovarianceMatrix includes stats

    TMatrixT<double> background_full_covariance_matrix = m_sbnchi_grid[0]->CalcCovarianceMatrix(m_full_fractional_covariance_matrix, *m_tvec_background_spectrum);
    TMatrixT<double> background_collapsed_covariance_matrix(m_background_spectrum->num_bins_total_compressed, m_background_spectrum->num_bins_total_compressed);
    m_sbnchi_grid[0]->CollapseModes(background_full_covariance_matrix, background_collapsed_covariance_matrix);    
    TMatrixT<double> inverse_background_collapsed_covariance_matrix = m_sbnchi_grid[0]->InvertMatrix(background_collapsed_covariance_matrix);   


    SBNspec * m_observed_spectrum;
    if(which_pt <0){
        m_observed_spectrum = m_background_spectrum;
    }else{
        m_observed_spectrum = m_cv_spec_grid.at(which_pt);   
    }

    for(size_t t =0; t < m_num_total_gridpoints; t++){

        std::cout<<"Starting on point "<<t<<"/"<<m_num_total_gridpoints<<std::endl;

        SBNspec * true_spec = m_cv_spec_grid.at(t); 
        SBNchi  * true_chi = m_sbnchi_grid.at(t); 
        std::vector<double> true_spec_vec =true_spec->full_vector;

        double chiSq;
        if(which_pt<0){
            chiSq = m_background_chi->CalcChi(true_spec);  //exclusion
        }else{
            chiSq = m_sbnchi_grid[which_pt]->CalcChi(true_spec);  //allowed
        }

        std::cout<<"ANS: "<<t<<" "<<chiSq;
        for(int k=0; k<m_vec_grid[t].size();k++){
            std::cout<<" "<<m_vec_grid[t][k];
        }
        std::cout<<std::endl;
    }




    return 0;
};




int SBNfeld::GlobalScan(SBNspec * observed_spectrum){

    //Ok take the background only spectrum and form a background only covariance matrix. CalcCovarianceMatrix includes stats

    TMatrixT<double> background_full_covariance_matrix = m_sbnchi_grid[0]->CalcCovarianceMatrix(m_full_fractional_covariance_matrix, *m_tvec_background_spectrum);
    TMatrixT<double> background_collapsed_covariance_matrix(m_background_spectrum->num_bins_total_compressed, m_background_spectrum->num_bins_total_compressed);
    m_sbnchi_grid[0]->CollapseModes(background_full_covariance_matrix, background_collapsed_covariance_matrix);    
    TMatrixT<double> inverse_background_collapsed_covariance_matrix = m_sbnchi_grid[0]->InvertMatrix(background_collapsed_covariance_matrix);   

    //first calculate a chi^2 minimum and BF

    observed_spectrum->CollapseVector();
    std::vector<double> ans = this->PerformIterativeFit(observed_spectrum->f_collapsed_vector,0,inverse_background_collapsed_covariance_matrix);
    double chi_min = ans[2];
    size_t bf = ans[0];
    
    std::cout<<"Best fit is point "<<bf<<" chi^2 min of "<<chi_min<<"  is point ";
    for(int k=0; k<m_vec_grid[bf].size();k++){
            std::cout<<" "<<m_vec_grid[bf][k];
    }
    std::cout<<std::endl;


	SBNspec * bf_spec = m_cv_spec_grid[bf];
	
    bf_spec->CollapseVector();
    TVectorT<double>* m_tvec_bf_spectrum = new TVectorT<double>(bf_spec->full_vector.size(), &(bf_spec->full_vector)[0]);


	TMatrixT<double> bf_full_covariance_matrix = m_sbnchi_grid[0]->CalcCovarianceMatrix(m_full_fractional_covariance_matrix, *m_tvec_bf_spectrum);
    	TMatrixT<double> bf_collapsed_covariance_matrix(bf_spec->num_bins_total_compressed, bf_spec->num_bins_total_compressed);
    	m_sbnchi_grid[0]->CollapseModes(bf_full_covariance_matrix, bf_collapsed_covariance_matrix);    
    	TMatrixT<double> inverse_bf_collapsed_covariance_matrix = m_sbnchi_grid[0]->InvertMatrix(bf_collapsed_covariance_matrix);   



    for(size_t t =0; t < m_num_total_gridpoints; t++){
        std::cout<<"Starting on point "<<t<<"/"<<m_num_total_gridpoints<<std::endl;
        //SBNchi *test_chi = m_sbnchi_grid.at(t); 
	SBNspec *test_spec = m_cv_spec_grid.at(t);
	test_spec->CollapseVector();


        //double chiSq = test_chi->CalcChi(observed_spectrum); 
        //double chiSq   = this->CalcChi(observed_spectrum->f_collapsed_vector, test_spec->collapsed_vector, inverse_bf_collapsed_covariance_matrix, true);
        double chiSq   = this->CalcChi(observed_spectrum->f_collapsed_vector, test_spec->collapsed_vector, inverse_bf_collapsed_covariance_matrix, false);


	double deltaChi = chiSq-chi_min;

        std::cout<<"ANS: "<<t<<" "<<chiSq<<" "<<deltaChi;
        for(int k=0; k<m_vec_grid[t].size();k++){
            std::cout<<" "<<m_vec_grid[t][k];
        }
        std::cout<<std::endl;
    }

    std::cout << "Check 1" << std::endl;
    m_cv_spec_grid[bf]->CompareSBNspecs(observed_spectrum, "Ans");
    std::cout << "Check 2" << std::endl;
    return 0;
};



int SBNfeld::SetStatOnly(){
    m_bool_stat_only = true;
    return 0;
}

int SBNfeld::RasterScan(){
    return 0;
};
int SBNfeld::SetNumUniverses(int n){
    m_num_universes = n;
    return 0;
}  



