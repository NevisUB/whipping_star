#include "SBNosc.h"
using namespace sbn;

SBNosc::SBNosc(std::string name, std::string whichxml) : SBNspec(name, whichxml) {
	working_model.zero();
	mass_step_size = 0.04;
	which_mode = BOTH_ONLY;
}
SBNosc::SBNosc(SBNspec & specin): SBNspec(specin){
	working_model.zero();
	which_mode = BOTH_ONLY;
}


SBNosc::SBNosc(std::string name, std::string whichxml, NeutrinoModel in) : SBNosc(name, whichxml) {
	LoadModel(in);

}

int SBNosc::LoadModel(NeutrinoModel in){
	working_model = in;
	calcMassSplittings();
	return 0;
}

/*************************************************
 * for a given working_model.
 * Calculate how many mass splittings, and which "type" it is,
 *
 * **********************************************/

int SBNosc::calcMassSplittings(){
	mass_splittings.clear();

	double fix41=round(log10((working_model.dm41Sq))/mass_step_size)*mass_step_size;
	double fix51=round(log10((working_model.dm51Sq))/mass_step_size)*mass_step_size;
	double fix61=round(log10((working_model.dm61Sq))/mass_step_size)*mass_step_size;

	double round54 = round(log10(fabs(working_model.dm54Sq))/mass_step_size)*mass_step_size;
	double round64 = round(log10(fabs(working_model.dm64Sq))/mass_step_size)*mass_step_size;
	double round65 = round(log10(fabs(working_model.dm65Sq))/mass_step_size)*mass_step_size;

	if(working_model.numsterile == 1)
	{
		mass_splittings.push_back( std::make_pair(fix41,41 ));
	}
       	else if(working_model.numsterile ==2)
	{
		mass_splittings.push_back( std::make_pair (fix41,41));
		mass_splittings.push_back( std::make_pair (fix51,51));
		mass_splittings.push_back( std::make_pair (round54,54));
	}
	else if(working_model.numsterile ==3)
	{
		mass_splittings.push_back( std::make_pair (fix41,41));
		mass_splittings.push_back( std::make_pair (fix51,51));
		mass_splittings.push_back( std::make_pair (fix61,61));

		mass_splittings.push_back( std::make_pair (round54,54));
		mass_splittings.push_back( std::make_pair (round64,64));
		mass_splittings.push_back( std::make_pair (round65,65));
	}




	return 0;
}




int SBNosc::OscillateThis(std::string tag){
	this->CalcFullVector();
	this->CollapseVector();

	calcMassSplittings();

	for(auto ms: mass_splittings){

			//this is wrong
			std::string name_sinsq = tag +"_SINSQ_dm_"+working_model.mass_tag+".SBNspec.root";
			std::string name_sin = tag +"_SIN_dm_"+working_model.mass_tag+".SBNspec.root";


			SBNspec single_frequency(name_sin , xmlname , false);
			SBNspec single_frequency_square(name_sinsq , xmlname ,false);


			if(has_been_scaled){
				single_frequency.Scale(scale_hist_name, scale_hist_val);
				single_frequency_square.Scale(scale_hist_name, scale_hist_val);
			}

			double prob_mumu, prob_ee, prob_mue, prob_mue_sq, prob_muebar, prob_muebar_sq;

			int which_dm = ms.second;

			switch (which_mode)
				{
					case APP_ONLY: //Strictly nu_e app only
						prob_mumu =0;
						prob_ee   =0;
						prob_mue = working_model.oscAmp(2,1,which_dm,1);
						prob_mue_sq = working_model.oscAmp(2,1,which_dm,2);
						prob_muebar = working_model.oscAmp(-2,-1,which_dm,1);
						prob_muebar_sq = working_model.oscAmp(-2,-1,which_dm,2);
						break;
					case DIS_ONLY: //Strictly nu_mu dis only
						prob_mumu = working_model.oscAmp(2,2,which_dm,2);
						prob_ee = 0;
						prob_mue = 0;
						prob_mue_sq =0;
						prob_muebar =0;
						prob_muebar_sq =0;
						break;
					case BOTH_ONLY: // This allows for both nu_e dis/app and nu_mu dis
						prob_mumu = working_model.oscAmp(2,2,which_dm,2);
						prob_ee = working_model.oscAmp(1,1,which_dm,2);
						prob_mue = working_model.oscAmp(2,1,which_dm,1);
						prob_mue_sq = working_model.oscAmp(2,1,which_dm,2);
						prob_muebar = working_model.oscAmp(-2,-1,which_dm,1);
						prob_muebar_sq = working_model.oscAmp(-2,-1,which_dm,2);
						break;
					case WIERD_ONLY: // A strange version where nu_e can appear but not disapear
						prob_mumu =working_model.oscAmp(2,2,which_dm,2);
						prob_ee = 0.0;
						prob_mue = working_model.oscAmp(2,1,which_dm,1);
						prob_mue_sq = working_model.oscAmp(2,1,which_dm,2);
						prob_muebar = working_model.oscAmp(-2,-1,which_dm,1);
						prob_muebar_sq = working_model.oscAmp(-2,-1,which_dm,2);
						break;
					case DISE_ONLY: // A strange version where nu_e can appear but not disapear
						prob_mumu = 0.0;
						prob_ee = working_model.oscAmp(1,1,which_dm,2);
						prob_mue = 0;
						prob_mue_sq = 0;
						prob_muebar = 0;
						prob_muebar_sq = 0;
						break;
				}


			std::cout<<" MassSplit "<<working_model.mass_tag<<" mm: "<<prob_mumu<<" ee: "<<prob_ee<<" mue: "<<prob_mue<<" mueSQ: "<<prob_mue_sq<<" mubar: "<<prob_muebar<<" muebarSQ: "<<prob_muebar_sq<<std::endl;


			for(int i=0; i<num_channels; i++){
				for(int j=0; j<num_subchannels.at(i); j++){
					int osc_pattern = subchannel_osc_patterns.at(i).at(j);
					double osc_amp = 0;
					double osc_amp_sq = 0;
                    std::cout<<i<<" "<<j<<" "<<osc_pattern<<std::endl;
					switch(osc_pattern){
						case 11:
							osc_amp_sq = prob_ee;
							break;
						case -11:
							osc_amp_sq = prob_ee;
							break;
						case 22:
							osc_amp_sq = prob_mumu;
							break;
						case -22:
							osc_amp_sq = prob_mumu;
							break;
						case 21:
							osc_amp = prob_mue;
							osc_amp_sq = prob_mue_sq;
							break;
						case -21:
							osc_amp = prob_muebar;
							osc_amp_sq = prob_muebar_sq;
							break;
						case 0:
						default:
							break;
					}

					single_frequency.Scale(channel_names.at(i)+"_"+subchannel_names.at(i).at(j), osc_amp );
					single_frequency_square.Scale(channel_names.at(i)+"_"+subchannel_names.at(i).at(j), osc_amp_sq );
                    std::cout<<"   "<<channel_names.at(i)<<" "<<subchannel_names.at(i).at(j)<<" "<<osc_amp<<" "<<osc_amp_sq<<std::endl;
				}
			}



			/*single_frequency.Scale("elike_fulloscnue", prob_mue);
			single_frequency.Scale("elike_fulloscbarnue", prob_muebar);
			single_frequency.Scale("elike_intrinsic", prob_ee);
			single_frequency.Scale("elike_mismuon", prob_mumu);
			single_frequency.Scale("elike_misphoton",0.0);
			single_frequency.Scale("elike_dirt",0.0);
			single_frequency.Scale("elike_cosmic",0.0);
			single_frequency.Scale("mlike_intrinsic", prob_mumu);
			single_frequency.Scale("mlike_misncpion",0.0);

			single_frequency_square.Scale("elike_fulloscnue", prob_mue_sq);
			single_frequency_square.Scale("elike_fulloscbarnue", prob_muebar_sq);
			single_frequency_square.Scale("elike_intrinsic", 0.0);
			single_frequency_square.Scale("elike_mismuon", 0.0);
			single_frequency_square.Scale("elike_misphoton",0.0);
			single_frequency_square.Scale("elike_dirt",0.0);
			single_frequency_square.Scale("elike_cosmic",0.0);
			single_frequency_square.Scale("mlike_intrinsic", 0.0);
			single_frequency_square.Scale("mlike_misncpion",0.0);
			*/

			this->Add(&single_frequency);
			this->Add(&single_frequency_square);


	}//Done looping over

	this->CalcFullVector();
	this->CollapseVector();
	this->CalcErrorVector();

	return 0;
};

std::vector<std::vector<double>> SBNosc::Oscillate(std::string tag, double scale){

	std::vector<std::vector<double>> tmp = this->Oscillate(tag);
	for(auto &e: tmp){
	     for(auto &v:e){
		v=v*scale;
	     }
	}

	return tmp;
}

/*
std::vector<double> SBNosc::OscillateWithAmp(double amp, double amp_sq){

		this->CalcFullVector();
		this->CollapseVector();

	std::vector<double > temp = collapsed_vector;


	for(auto ms: mass_splittings){
			char namei[200];

			sPrintf(namei, "%sprecomp/SBN_SIN_%2.2f",data_path.c_str(), ms.first );
			SBNspec single_frequency(namei , xmlname, false);

			sPrintf(namei, "%sprecomp/SBN_SINSQ_%2.2f",data_path.c_str(), ms.first );
			SBNspec single_frequency_square(namei , xmlname, false);

			if(has_been_scaled){
				single_frequency.Scale(scale_hist_name, scale_hist_val);
				single_frequency_square.Scale(scale_hist_name, scale_hist_val);
			}

			double prob_mumu, prob_ee, prob_mue, prob_mue_sq, prob_muebar, prob_muebar_sq;

			int which_dm = ms.second;

			for(int i=0; i<num_channels; i++){
				for(int j=0; j<num_subchannels.at(i); j++){
					int osc_pattern = subchannel_osc_patterns.at(i).at(j);
					double osc_amp = 0;
					double osc_amp_sq = 0;
					switch(osc_pattern){
						case 0:
							break;
						default:
							osc_amp = amp;
							osc_amp_sq = amp_sq;
							break;
					}

					single_frequency.Scale(channel_names.at(i)+"_"+subchannel_names.at(i).at(j), osc_amp);
					single_frequency_square.Scale(channel_names.at(i)+"_"+subchannel_names.at(i).at(j), osc_amp_sq );

				}
			}


			single_frequency.CalcFullVector();
			single_frequency.CollapseVector();


			single_frequency_square.CalcFullVector();
			single_frequency_square.CollapseVector();

			for(int i=0;i<temp.size(); i++){
				temp[i] += single_frequency.collapsed_vector[i];
				temp[i] += single_frequency_square.collapsed_vector[i];
			}

	}//Done looping over

	return temp;
};
*/

std::vector< std::vector<double>> SBNosc::Oscillate(std::string tag){

    return this->Oscillate(tag,true);
}
std::vector<std::vector<double>> SBNosc::Oscillate(std::string tag, bool return_compressed){

     this->CalcFullVector();
     this->CollapseVector();
     this->CalcErrorVector();

     calcMassSplittings();

     std::vector<double> temp;
     std::vector<double> temp_err;
     std::vector<double> temp_err_square;
    
    if(return_compressed){
        temp = collapsed_vector;
	temp_err = collapsed_err_vector;
    }else {
        temp = full_vector;
	temp_err = full_err_vector;
    }

  
    // this vector saves the sum of squares of errors
    temp_err_square.clear();
    for(auto &v: temp_err){
	temp_err_square.push_back(pow(v, 2.0));
    }


	for(auto ms: mass_splittings){

			std::string name_sinsq = tag +"_SINSQ_dm_"+working_model.mass_tag+".SBNspec.root";
			std::string name_sin = tag +"_SIN_dm_"+working_model.mass_tag+".SBNspec.root";

			SBNspec single_frequency(name_sin , xmlname , false);
			SBNspec single_frequency_square(name_sinsq , xmlname ,false);

			if(has_been_scaled){
				single_frequency.Scale(scale_hist_name, scale_hist_val);
				single_frequency_square.Scale(scale_hist_name, scale_hist_val);
			}

    			single_frequency.CalcFullVector();
			single_frequency_square.CalcFullVector();
            		//single_frequency.PrintFullVector();
            		//single_frequency_square.PrintFullVector();


			double prob_mumu, prob_ee, prob_mue, prob_mue_sq, prob_muebar, prob_muebar_sq;

			int which_dm = ms.second;

			switch (which_mode)
				{
					case APP_ONLY: //Strictly nu_e app only
						prob_mumu =0;
						prob_ee   =0;
						prob_mue = working_model.oscAmp(2,1,which_dm,1);
						prob_mue_sq = working_model.oscAmp(2,1,which_dm,2);
						prob_muebar = working_model.oscAmp(-2,-1,which_dm,1);
						prob_muebar_sq = working_model.oscAmp(-2,-1,which_dm,2);
						break;
					case DIS_ONLY: //Strictly nu_mu dis only
						prob_mumu = working_model.oscAmp(2,2,which_dm,2);
						prob_ee = 0;
						prob_mue = 0;
						prob_mue_sq =0;
						prob_muebar =0;
						prob_muebar_sq =0;
						break;
					case BOTH_ONLY: // This allows for both nu_e dis/app and nu_mu dis
						prob_mumu = working_model.oscAmp(2,2,which_dm,2);
						prob_ee = working_model.oscAmp(1,1,which_dm,2);
						prob_mue = working_model.oscAmp(2,1,which_dm,1);
						prob_mue_sq = working_model.oscAmp(2,1,which_dm,2);
						prob_muebar = working_model.oscAmp(-2,-1,which_dm,1);
						prob_muebar_sq = working_model.oscAmp(-2,-1,which_dm,2);
						break;
					case WIERD_ONLY: // A strange version where nu_e can appear but not disapear
						prob_mumu =working_model.oscAmp(2,2,which_dm,2);
						prob_ee = 0.0;
						prob_mue = working_model.oscAmp(2,1,which_dm,1);
						prob_mue_sq = working_model.oscAmp(2,1,which_dm,2);
						prob_muebar = working_model.oscAmp(-2,-1,which_dm,1);
						prob_muebar_sq = working_model.oscAmp(-2,-1,which_dm,2);
						break;
					case DISE_ONLY: // A strange version where nu_e can appear but not disapear
						prob_mumu = 0.0;
						prob_ee = working_model.oscAmp(1,1,which_dm,2);
						prob_mue = 0;
						prob_mue_sq = 0;
						prob_muebar = 0;
						prob_muebar_sq = 0;
						break;

				}


			for(int i=0; i<num_channels; i++){
				for(int j=0; j<num_subchannels.at(i); j++){
					int osc_pattern = subchannel_osc_patterns.at(i).at(j);
					double osc_amp = 0;
					double osc_amp_sq = 0;
					switch(osc_pattern){
						case 11:
							osc_amp_sq = prob_ee;
							break;
						case -11:
							osc_amp_sq = prob_ee;
							break;
						case 22:
							osc_amp_sq = prob_mumu;
							break;
						case -22:
							osc_amp_sq = prob_mumu;
							break;
						case 21:
							osc_amp = prob_mue;
							osc_amp_sq = prob_mue_sq;
							break;
						case -21:
							osc_amp = prob_muebar;
							osc_amp_sq = prob_muebar_sq;
							break;
						case 0:
						default:
							break;
					}

					single_frequency.Scale(channel_names.at(i)+"_"+subchannel_names.at(i).at(j), osc_amp);
					single_frequency_square.Scale(channel_names.at(i)+"_"+subchannel_names.at(i).at(j), osc_amp_sq );
				}
			}


	    	//std::cout<<"mm: "<<prob_mumu<<" ee: "<<prob_ee<<" mue: "<<prob_mue<<" mueSQ: "<<prob_mue_sq<<" mubar: "<<prob_muebar<<" muebarSQ: "<<prob_muebar_sq<<std::endl;
			/*
			single_frequency.Scale("elike_fulloscnue", prob_mue);
			single_frequency.Scale("elike_fulloscbarnue", prob_muebar);
			single_frequency.Scale("elike_intrinsic", prob_ee);
			single_frequency.Scale("elike_mismuon", prob_mumu);
			single_frequency.Scale("elike_misphoton",0.0);
			single_frequency.Scale("elike_dirt",0.0);
			single_frequency.Scale("elike_cosmic",0.0);
			single_frequency.Scale("mlike_intrinsic", prob_mumu);
			single_frequency.Scale("mlike_misncpion",0.0);

			single_frequency_square.Scale("elike_fulloscnue", prob_mue_sq);
			single_frequency_square.Scale("elike_fulloscbarnue", prob_muebar_sq);
			single_frequency_square.Scale("elike_intrinsic", 0.0);
			single_frequency_square.Scale("elike_mismuon", 0.0);
			single_frequency_square.Scale("elike_misphoton",0.0);
			single_frequency_square.Scale("elike_dirt",0.0);
			single_frequency_square.Scale("elike_cosmic",0.0);
			single_frequency_square.Scale("mlike_intrinsic", 0.0);
			single_frequency_square.Scale("mlike_misncpion",0.0);
			*/


			single_frequency.CalcFullVector();
			single_frequency.CollapseVector();
			single_frequency.CalcErrorVector();

			single_frequency_square.CalcFullVector();
			single_frequency_square.CollapseVector();
			single_frequency_square.CalcErrorVector();

		        if(return_compressed){ 
				for(int i=0;i<temp.size(); i++){
					temp[i] += single_frequency.collapsed_vector[i];
					temp_err_square[i] += pow(single_frequency.collapsed_err_vector[i], 2.0);

					temp[i] += single_frequency_square.collapsed_vector[i];
					temp_err_square[i] += pow(single_frequency_square.collapsed_err_vector[i],2.0);
				}

		        }else{
				for(int i=0;i<temp.size(); i++){
					temp[i] += single_frequency.full_vector[i];
					temp_err_square[i] += pow(single_frequency.full_err_vector[i], 2.0);

					temp[i] += single_frequency_square.full_vector[i];
					temp_err_square[i] += pow(single_frequency_square.full_err_vector[i], 2.0);
				}
			}

	}//Done looping over


 	for(int i=0; i< temp.size(); i++)
		temp_err[i] = sqrt(temp_err_square[i]);

        std::vector< std::vector<double>> bundle{temp, temp_err};
	return bundle;
};




int SBNosc::SetMode(int in){
	which_mode = in;

return in;
}

void SBNosc::SetAppMode(){
	std::cout << "Running APP MODE!!!" << std::endl;
	SetMode(APP_ONLY);
}

void SBNosc::SetDisMode(){
	SetMode(DIS_ONLY);
}

void SBNosc::SetBothMode(){
	std::cout << "Running BOTH MODE!!!" << std::endl;
	SetMode(BOTH_ONLY);
}

void SBNosc::SetWierdMode(){
	SetMode(WIERD_ONLY);
}

void SBNosc::SetDisEMode(){
	SetMode(DISE_ONLY);
}
