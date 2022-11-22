#include "SBNconfig.h"
using namespace sbn;


// Constructor that takes the contents of a file.
SBNconfig::SBNconfig(const char * filedata, bool isverbose, bool useuniverse = true) {
    use_universe = useuniverse;  
    LoadFromXML(filedata, isverbose, use_universe);
}//end constructor


//standard constructor given an .xml
SBNconfig::SBNconfig(std::string whichxml, bool isverbose, bool useuniverse): xmlname(whichxml) {

    otag = "SBNconfig::SBNconfig\t||\t";
    is_verbose = isverbose;
    use_universe = useuniverse;  

    std::string line, text;
    std::ifstream in(whichxml);
    while(std::getline(in, line))  text += line + "\n";
    const char* xmldata = text.c_str();
    LoadFromXML(xmldata, isverbose, use_universe);

}//end constructor

SBNconfig::SBNconfig(std::string whichxml, bool isverbose): SBNconfig(whichxml, isverbose, true){} 
SBNconfig::SBNconfig(std::string whichxml): SBNconfig(whichxml, true, true) {}

//Constructor to build a SBNspec from scratch, not reccomended often
SBNconfig::SBNconfig(std::vector<std::string> modein, std::vector<std::string> detin, std::vector<std::string> chanin, std::vector<std::vector<std::string>> subchanin, std::vector<std::vector<double>> binin){

    otag = "SBNconfig::SBNconfig\t|| ";
    is_verbose = true;

    num_detectors = detin.size();
    num_channels = chanin.size();
    num_modes = modein.size();

    if(subchanin.size() != chanin.size()){
        std::cout<<otag<<"ERROR SUBCHAN.size() != chanin.size()"<<std::endl;
        exit(EXIT_FAILURE);
    }

    for(auto sb: subchanin){
        num_subchannels.push_back( sb.size());
    }

    for(auto bn: binin){
        num_bins.push_back(bn.size()-1);
    }

    this->CalcTotalBins();

    xmlname = "NULL";


    //the xml names are the way we track which channels and subchannels we want to use later
    mode_names = modein;
    detector_names = detin;
    channel_names = chanin;
    subchannel_names = subchanin;

    for(auto c: chanin){
        channel_bool.push_back(true);
    }
    for(auto m: modein){
        mode_bool.push_back(true);
    }
    for(auto d: detin){
        detector_bool.push_back(true);
    }
    for(auto c: chanin){
        std::vector<bool> tml;
        for(int i=0; i< num_subchannels.size(); i++){
            tml.push_back(true);
        }
        subchannel_bool.push_back(tml);
    }

    //self explanatory
    bin_edges = binin;

    for(auto binedge: bin_edges){
        std::vector<double> binwidth;
        for(int b = 0; b<binedge.size()-1; b++){
            binwidth.push_back(fabs(binedge.at(b)-binedge.at(b+1)));
        }
        bin_widths.push_back(binwidth);
    }

    // The order is IMPORTANT its the same as defined in xml
    for(auto m: mode_names){
        for(auto d: detector_names){
            for(int c = 0; c< num_channels; c++){
                for(auto sb: subchannel_names.at(c)){
                    std::string tmp = m +"_"+ d+"_"+channel_names.at(c)+"_"+sb;
                    fullnames.push_back(tmp);
                }
            }
        }
    }


    for(int i=0; i< num_bins_total; i++){
        used_bins.push_back(i);
    }


};


int SBNconfig::CalcTotalBins(){
    // These variables are important
    // They show how big each mode block and decector block are, for any given number of channels/subchannels
    // both before and after compression!

    //needs to be calculated AFTER usage bool removal above

    num_bins_detector_block = 0;
    num_bins_detector_block_compressed = 0;

    for(auto i: channel_used){
        num_bins_detector_block += num_subchannels[i]*num_bins[i];
        num_bins_detector_block_compressed += num_bins[i];
    }

    num_bins_mode_block = num_bins_detector_block*num_detectors;
    num_bins_mode_block_compressed = num_bins_detector_block_compressed*num_detectors;

    num_bins_total = num_modes*num_bins_mode_block;
    num_bins_total_compressed = num_modes*num_bins_mode_block_compressed;

    return 0;
}


int SBNconfig::LoadFromXML(const char * filedata, bool isverbose, bool useuniverse = true){
    otag = "SBNconfig::LoadFromXML\t||\t";
    is_verbose = isverbose;
    log<LOG_DEBUG>(L"-------------------------------------------------------------------------");

    has_oscillation_patterns = false;
    use_universe = useuniverse;  

    //max subchannels 100?
    subchannel_plotnames.resize(100);
    subchannel_datas.resize(100);
    subchannel_names.resize(100);
    subchannel_bool.resize(100);
    subchannel_used.resize(100);
    subchannel_osc_patterns.resize(100);
    char *end;

    //Setup TiXml documents
    TiXmlDocument doc;
    bool loadOkay = doc.Parse(filedata, 0, TIXML_ENCODING_UTF8);

    try{
        if(loadOkay) log<LOG_INFO>(L"%1% || Correctly loaded and parsed XML, continuing") % __func__;
        else throw 404;    
    }
    catch (int ernum) {
        log<LOG_ERROR>(L"%1% || ERROR: Failed to load XML configuration file. @ line %2% in %3% ") % __func__ % __LINE__  % __FILE__;
        log<LOG_ERROR>(L"This generally means broken brackets or attribute syntax in xml itself.");
        log<LOG_ERROR>(L"Terminating.");
        exit(EXIT_FAILURE);
    }


    TiXmlHandle hDoc(&doc);
    // we have Modes, Detectors, Channels
    TiXmlElement *pMode, *pDet, *pChan, *pCov, *pMC, *pData,*pPOT, *pWeiMaps, *pList, *pSpec, *pShapeOnlyMap;


    //Grab the first element. Note very little error checking here! make sure they exist.
    pMode = doc.FirstChildElement("mode");
    pDet =  doc.FirstChildElement("detector");
    pChan = doc.FirstChildElement("channel");
    pCov  = doc.FirstChildElement("covariance");
    pMC   = doc.FirstChildElement("MultisimFile");
    pData   = doc.FirstChildElement("data");
    pPOT    = doc.FirstChildElement("plotpot");
    pWeiMaps = doc.FirstChildElement("WeightMaps");
    pList = doc.FirstChildElement("variation_list");
    pSpec = doc.FirstChildElement("varied_spectrum");
    pShapeOnlyMap = doc.FirstChildElement("ShapeOnlyUncertainty");

    if(!pMode){
        log<LOG_ERROR>(L"%1% || ERROR: Need at least 1 mode defined in xml.@ line %2% in %3% ") % __func__ % __LINE__  % __FILE__;
        log<LOG_ERROR>(L"Terminating.");
        exit(EXIT_FAILURE);
    }else{
        while(pMode){
            // What modes are we running in (e.g nu, nu bar, horn current=XXvolts....) Can have as many as we want
            const char* mode_name= pMode->Attribute("name");
            if(mode_name==NULL){
                log<LOG_ERROR>(L"%1% || Modes need a name! Please define a name attribute for all modes. @ line %2% in %3% ") % __func__ % __LINE__  % __FILE__;
                log<LOG_ERROR>(L"Terminating.");
                exit(EXIT_FAILURE);
            }else{
                mode_names.push_back(mode_name);
            }

            const char* mode_plotname= pMode->Attribute("plotname");
            if(mode_plotname==NULL){
                mode_plotnames.push_back(mode_names.back());
            }else{
                mode_plotnames.push_back(mode_plotname);
            }

            const char* use_mode = pMode->Attribute("use");
            if(use_mode==NULL){
                mode_bool.push_back(1);
            }else{
                mode_bool.push_back(strtod(use_mode,&end) );
            }

            pMode = pMode->NextSiblingElement("mode");
            log<LOG_DEBUG>(L"%1% || Loading Mode %2% with use_bool %3%  ") % __func__ % mode_names.back().c_str() % mode_bool.back();

        }
    }

    // How many detectors do we want!
    if(!pDet){
        log<LOG_ERROR>(L"%1% || ERROR: Need at least 1 detector defined in xml.@ line %2% in %3% ") % __func__ % __LINE__  % __FILE__;
        log<LOG_ERROR>(L"Terminating.");
        exit(EXIT_FAILURE);

    }else{

        while(pDet){

            const char* detector_name= pDet->Attribute("name");
            if(detector_name==NULL){
                log<LOG_ERROR>(L"%1% || ERROR: Need all detectors to have a name attribute @ line %2% in %3% ") % __func__ % __LINE__  % __FILE__;
                log<LOG_ERROR>(L"Terminating.");
                exit(EXIT_FAILURE);
            }else{
                detector_names.push_back(detector_name);
            }

            const char* detector_plotname = pDet->Attribute("plotname");
            if(detector_plotname==NULL){
                detector_plotnames.push_back(detector_names.back());
            }else{
                detector_plotnames.push_back(detector_plotname);
            }


            const char* use_detector = pDet->Attribute("use");
            if(use_detector==NULL){
                detector_bool.push_back(1);
            }else{
                detector_bool.push_back(strtod(use_detector,&end) );
            }

            pDet = pDet->NextSiblingElement("detector");
            log<LOG_DEBUG>(L"%1% || Loading Det %2% with use_bool %3%  ") % __func__ % detector_names.back().c_str() % detector_bool.back();

        }
    }

    //How many channels do we want! At the moment each detector must have all channels
    int nchan = 0;
    if(!pChan){
        log<LOG_ERROR>(L"%1% || ERROR: Need at least 1 channel defined in xml.@ line %2% in %3% ") % __func__ % __LINE__  % __FILE__;
        log<LOG_ERROR>(L"Terminating.");
        exit(EXIT_FAILURE);
    }else{


        while(pChan){
            // Read in how many bins this channel uses

            const char* channel_name= pChan->Attribute("name");
            if(channel_name==NULL){
                log<LOG_ERROR>(L"%1% || ERROR: Need all channels to have names in xml.@ line %2% in %3% ") % __func__ % __LINE__  % __FILE__;
                log<LOG_ERROR>(L"Terminating.");
                exit(EXIT_FAILURE);
            }else{
                channel_names.push_back(channel_name);
            }

            const char* channel_plotname= pChan->Attribute("plotname");
            if(channel_plotname==NULL){
                channel_plotnames.push_back(channel_names.back());
            }else{
                channel_plotnames.push_back(channel_plotname);
            }


            const char* channel_unit= pChan->Attribute("unit");
            if(channel_unit==NULL){
                channel_units.push_back("");
            }else{
                channel_units.push_back(channel_unit);
            }

            const char* channel_bool_tmp= pChan->Attribute("use");
            if(channel_bool_tmp==NULL){
                channel_bool.push_back(1);
            }else{
                channel_bool.push_back(strtod(channel_bool_tmp,&end));
            }

            log<LOG_DEBUG>(L"%1% || Loading Channel %2% with use_bool %3%  ") % __func__ % channel_names.back().c_str() % channel_bool.back();


            // What are the bin edges and bin widths (bin widths just calculated from edges now)
            TiXmlElement *pBin = pChan->FirstChildElement("bins");
            std::stringstream iss(pBin->Attribute("edges"));

            double number;
            std::vector<double> binedge;
            std::vector<double> binwidth;
            std::string binstring = "";
            while ( iss >> number ){
                binedge.push_back( number );
                binstring+=" "+std::to_string(number);
            }

            log<LOG_DEBUG>(L"%1% || Loading Bins with edges %2%  ") % __func__ % binstring.c_str();

            for(int b = 0; b<binedge.size()-1; b++){
                binwidth.push_back(fabs(binedge.at(b)-binedge.at(b+1)));
            }

            num_bins.push_back(binedge.size()-1);

            bin_edges.push_back(binedge);
            bin_widths.push_back(binwidth);


            // Now loop over all this channels subchanels. Not the names must be UNIQUE!!
            TiXmlElement *pSubChan;

            pSubChan = pChan->FirstChildElement("subchannel");
            int nsubchan=0;
            while(pSubChan){

                const char* subchannel_name= pSubChan->Attribute("name");
                if(subchannel_name==NULL){
                    log<LOG_ERROR>(L"%1% || ERROR: Subchannels need a name in xml.@ line %2% in %3% ") % __func__ % __LINE__  % __FILE__;
                    log<LOG_ERROR>(L"Terminating.");
                    exit(EXIT_FAILURE);

                }else{
                    subchannel_names[nchan].push_back(subchannel_name);
                    log<LOG_DEBUG>(L"%1% || Subchannel Starting:  %2%") % __func__ % subchannel_names.at(nchan).back().c_str() ;

                }


                const char* subchannel_plotname= pSubChan->Attribute("plotname");
                if(subchannel_plotname==NULL){
                    subchannel_plotnames[nchan].push_back(subchannel_names[nchan].back());
                }else{
                    subchannel_plotnames[nchan].push_back(subchannel_plotname);
                }

                const char* subchannel_data= pSubChan->Attribute("data");
                if(subchannel_data==NULL){
                    subchannel_datas[nchan].push_back(0);
                }else{
                    subchannel_datas[nchan].push_back(1);
                }


                const char* subchannel_bool_tmp= pSubChan->Attribute("use");
                if(subchannel_bool_tmp==NULL){
                    subchannel_bool[nchan].push_back(1);
                }else{
                    subchannel_bool[nchan].push_back(strtod(subchannel_bool_tmp,&end));
                }

                //0 means dont oscillate, 11 means electron disapearance, -11 means antielectron dis..etc..
                if(pSubChan->Attribute("osc"))
                {
                    has_oscillation_patterns = true;
                    subchannel_osc_patterns.at(nchan).push_back(strtod(pSubChan->Attribute("osc"), &end));
                }else{
                    has_oscillation_patterns = false;
                    subchannel_osc_patterns.at(nchan).push_back(0);
                }

                log<LOG_DEBUG>(L"%1% || Subchannel %2% with bool %3% and osc pattern %4% and isdata %5%") % __func__ % subchannel_names.at(nchan).back().c_str() % subchannel_bool.at(nchan).back() % subchannel_osc_patterns.at(nchan).back() % subchannel_datas.at(nchan).back();

                nsubchan++;
                pSubChan = pSubChan->NextSiblingElement("subchannel");
            }
            num_subchannels.push_back(nsubchan);

            nchan++;
            pChan = pChan->NextSiblingElement("channel");
        }
    }

    while(pPOT){
        const char* inplotpot = pPOT->Attribute("value");
        if(inplotpot==NULL){
            plot_pot = 1.0 ;
        }else{
            plot_pot = strtod(inplotpot,&end);
        }
        pPOT = pPOT->NextSiblingElement("plotpot");
    }

    // if wea re creating a covariance matrix using a ntuple and weights, here is the info
    if(pMC){
        while(pMC)
        {

            const char* tree = pMC->Attribute("treename");
            if(tree==NULL){
                log<LOG_ERROR>(L"%1% || ERROR: You must have an associated root TTree name for all MonteCarloFile tags.. eg. treename='events' @ line %2% in %3% ") % __func__ % __LINE__  % __FILE__;
                log<LOG_ERROR>(L"Terminating.");
                exit(EXIT_FAILURE);
            }else{
                montecarlo_name.push_back(tree);
            }


            const char* file = pMC->Attribute("filename");
            if(file==NULL){
                log<LOG_ERROR>(L"%1% || ERROR: You must have an associated root TFile name for all MonteCarloFile tags.. eg. filename='my.root' @ line %2% in %3% ") % __func__ % __LINE__  % __FILE__;
                log<LOG_ERROR>(L"Terminating.");
                exit(EXIT_FAILURE);
            }else{
                montecarlo_file.push_back(file);
            }


            const char* maxevents = pMC->Attribute("maxevents");
            if(maxevents==NULL){
                montecarlo_maxevents.push_back(1e16);
            }else{
                montecarlo_maxevents.push_back(strtod(maxevents,&end) );
            }


            const char* scale = pMC->Attribute("scale");
            if(scale==NULL){
                montecarlo_scale.push_back(1.0);
            }else{
                montecarlo_scale.push_back(strtod(scale,&end) );
            }

            const char* inpot = pMC->Attribute("pot");
            if(inpot==NULL){
                montecarlo_pot.push_back(-1.0);
            }else{
                montecarlo_pot.push_back(strtod(inpot,&end) );
            }

            const char* isfake = pMC->Attribute("fake");
            if(isfake==NULL){
                montecarlo_fake.push_back(false);
            }else{
                montecarlo_fake.push_back(true);
            }

    
            log<LOG_DEBUG>(L"%1% || MultisimFile %2%, treename: %3%  ") % __func__ % montecarlo_file.back().c_str() % montecarlo_name.back().c_str();


            /*
            //Currently take all parameter variations in at once. Depreciated code. 
            TiXmlElement *pParams = pMC->FirstChildElement("parameters");
            std::stringstream sss(pParams->Attribute("names"));

            std::vector<std::string> vstring;
            std::string nam;
            while ( sss >> nam) vstring.push_back( nam );

            parameter_names.push_back(vstring);
            */

            //Here we can grab some friend tree information
            TiXmlElement *pFriend;
            pFriend = pMC->FirstChildElement("friend");
            while(pFriend){


                std::string ffname;
                const char* friend_filename = pFriend->Attribute("filename");
                if(friend_filename==NULL){
                    ffname = montecarlo_file.back();
                }else{
                    ffname = friend_filename;
                }

                if(montecarlo_file_friend_treename_map.count(montecarlo_file.back())>0){

                    (montecarlo_file_friend_treename_map[montecarlo_file.back()]).push_back( pFriend->Attribute("treename") );
                    (montecarlo_file_friend_map[montecarlo_file.back()]).push_back(ffname);

                }else{
                    std::vector<std::string> temp_treename = {pFriend->Attribute("treename")};
                    std::vector<std::string> temp_filename = {ffname};

                    montecarlo_file_friend_treename_map[montecarlo_file.back()] = temp_treename;
                    montecarlo_file_friend_map[montecarlo_file.back()] = temp_filename;
                }
                pFriend = pFriend->NextSiblingElement("friend");
            }

            TiXmlElement *pBranch;
            pBranch = pMC->FirstChildElement("branch");

            std::vector<BranchVariable*> TEMP_branch_variables;
            while(pBranch){

                const char* bnam = pBranch->Attribute("name");
                const char* btype = pBranch->Attribute("type");
                const char* bhist = pBranch->Attribute("associated_subchannel");
                const char* bsyst = pBranch->Attribute("associated_systematic");
                const char* bcentral = pBranch->Attribute("central_value");
                const char* bwname = pBranch->Attribute("eventweight_branch_name");
                const char* badditional_weight = pBranch->Attribute("additional_weight");

                if(bwname== NULL){
                    log<LOG_WARNING>(L"%1% || WARNING: No eventweight branch name passed, defaulting to 'weights' @ line %2% in %3% ") % __func__ % __LINE__  % __FILE__;
                    montecarlo_eventweight_branch_names.push_back("weights");
                }else{
                    log<LOG_DEBUG>(L"%1% || Setting eventweight branch name %2%") %__func__ % bnam;
                    montecarlo_eventweight_branch_names.push_back(std::string(bwname));
                }

                if(bnam == NULL){
                    log<LOG_ERROR>(L"%1% || ERROR!: Each branch must include the name of the branch to use. @ line %2% in %3% ") % __func__ % __LINE__  % __FILE__;
                    log<LOG_ERROR>(L"%1% || ERROR!: e.g name = 'ereco' @ line %2% in %3% ") % __func__ % __LINE__  % __FILE__;
                    log<LOG_ERROR>(L"Terminating.");
                    exit(EXIT_FAILURE);
                }

                if(btype == NULL){
                    log<LOG_WARNING>(L"%1% || WARNING: No branch type has been specified, assuming double. @ line %2% in %3% ") % __func__ % __LINE__  % __FILE__;
                    btype= "double";
                }

                if(bhist == NULL){
                    log<LOG_ERROR>(L"%1% || Each branch must have an associated_subchannel to fill! On branch %4% : @ line %2% in %3% ") % __func__ % __LINE__  % __FILE__ % bnam;
                    log<LOG_ERROR>(L"%1% || e.g associated_subchannel='mode_det_chan_subchannel ") % __func__ % __LINE__  % __FILE__;
                    log<LOG_ERROR>(L"Terminating.");
                    exit(EXIT_FAILURE);
                }


                if(bsyst == NULL){
                    if(use_universe == false){
                        log<LOG_WARNING>(L"%1% || WARNING: No root file with unique systematic variation is provided ") % __func__;
                        log<LOG_ERROR>(L"%1% || ERROR! please provide what systematic variation this file correpsonds to!") % __func__;
                        log<LOG_ERROR>(L"Terminating.");
                        exit(EXIT_FAILURE);
                    }
                    systematic_name.push_back("");
                }else{
                    systematic_name.push_back(bsyst);	

                }


                if(badditional_weight == NULL){
                    montecarlo_additional_weight_bool.push_back(0);
                    montecarlo_additional_weight_names.push_back("");
                }else{
                    montecarlo_additional_weight_names.push_back(badditional_weight);
                    montecarlo_additional_weight_bool.push_back(1);
                    log<LOG_DEBUG>(L"%1% || Setting an additional weight for branch %2% using the branch %3% as a reweighting.") % __func__ % bnam %badditional_weight;

                }



                if((std::string)btype == "double"){
                    if(use_universe){
                        TEMP_branch_variables.push_back( new BranchVariable_d(bnam, btype, bhist ) );
                    } else  if((std::string)bcentral == "true"){
                        TEMP_branch_variables.push_back( new BranchVariable_d(bnam, btype, bhist,bsyst, true) );
                       log<LOG_DEBUG>(L"%1% || Setting as  CV for det sys.") % __func__ ;
                    } else {
                        TEMP_branch_variables.push_back( new BranchVariable_d(bnam, btype, bhist,bsyst, false) );
                       log<LOG_DEBUG>(L"%1% || Setting as individual (not CV) for det sys.") % __func__ ;
                    }
            }else{
                log<LOG_ERROR>(L"%1% || ERROR: currently only double, allowed for input branch variables (sorry!) i @ line %2% in %3% ") % __func__ % __LINE__  % __FILE__ % bnam;
                log<LOG_ERROR>(L"Terminating.");
                exit(EXIT_FAILURE);
            }

            std::string oscillate = "false";
            if(pBranch->Attribute("oscillate")!=NULL){
                oscillate =pBranch->Attribute("oscillate");
            }	

            if(oscillate == "false"){
                log<LOG_DEBUG>(L"%1% || Oscillations are OFF ") % __func__ ;
                TEMP_branch_variables.back()->SetOscillate(false);
            }else if(oscillate == "true"){
                log<LOG_DEBUG>(L"%1% || Oscillations are Set to  ON ") % __func__;
                TEMP_branch_variables.back()->SetOscillate(true);
                TEMP_branch_variables.back()->true_param_name = pBranch->Attribute("true_param_name");
                if(pBranch->Attribute("true_L_name") != NULL){
                    //for oscillation that needs both E and L
                    TEMP_branch_variables.back()->true_L_name = pBranch->Attribute("true_L_name");
                    log<LOG_DEBUG>(L"%1% || Oscillations using true param name:   %2% and baseline %3% ") % __func__ % pBranch->Attribute("true_param_name") % pBranch->Attribute("true_L_name") ;
                }else{
                    //for oscillations that only needs E, such as an energy-dependent scaling for single photon NCpi0!
                    log<LOG_DEBUG>(L"%1% || Oscillations using  Energy only dependent oscillation ( or shift/normalization)  %2% ") % __func__ % pBranch->Attribute("true_param_name") ;
                }
            }else{
                log<LOG_DEBUG>(L"%1% || Do Not Oscillate  ") % __func__  ;
                TEMP_branch_variables.back()->SetOscillate(false);
            }

            log<LOG_DEBUG>(L"%1% || Associated subchannel: %2% ") % __func__ % bhist;

            pBranch = pBranch->NextSiblingElement("branch");
            }
            branch_variables.push_back(TEMP_branch_variables);
            //next file
            pMC=pMC->NextSiblingElement("MultisimFile");
        }
    }

    if(!pList){
       log<LOG_DEBUG>(L"%1% || No Whitelist or Blacklist set, including ALL variations by default.") % __func__  ;
    }else{
        while(pList){

            TiXmlElement *pWhiteList = pList->FirstChildElement("whitelist");
            while(pWhiteList){
                std::string wt = std::string(pWhiteList->GetText());
                variation_whitelist[wt] = true; 
                log<LOG_DEBUG>(L"%1% || Whitelisting variations: %2%") % __func__ % wt.c_str() ;
                pWhiteList = pWhiteList->NextSiblingElement("whitelist");
            }

            TiXmlElement *pBlackList = pList->FirstChildElement("blacklist");
            while(pBlackList){
                std::string bt = std::string(pBlackList->GetText());
                variation_blacklist[bt] = true; 
                log<LOG_DEBUG>(L"%1% || Blacklisting variations: %2%") % __func__ % bt.c_str() ;
                pBlackList = pBlackList->NextSiblingElement("blacklist");
            }
            pList = pList->NextSiblingElement("variation_list");
        }
    }


    //weightMaps
    if(!pWeiMaps){
        log<LOG_DEBUG>(L"%1% || WeightMaps not set, all weights for all variations are 1 (individual branch weights still apply)") % __func__  ;
    }else{
        while(pWeiMaps){


            TiXmlElement *pVariation;
            pVariation = pWeiMaps->FirstChildElement("variation");
            while(pVariation){


                const char* w_pattern = pVariation->Attribute("pattern");
                const char* w_formula = pVariation->Attribute("weight_formula");
                const char* w_use = pVariation->Attribute("use");
                const char* w_mode = pVariation->Attribute("mode");

                if(w_pattern== NULL){
                    std::cout<<otag<<" ERROR! No pattern passed for this variation in WeightMaps'"<<std::endl;
                    exit(EXIT_FAILURE);
                }else{
                    if(is_verbose)std::cout<<otag<<" Loading WeightMaps Variation Pattern : "<<w_pattern<<std::endl;
                    weightmaps_patterns.push_back(std::string(w_pattern));
                }


                if(w_formula== NULL){
                    std::cout<<otag<<"Warning! No formula passed for this variation in WeightMaps. Setting to 1."<<std::endl;
                    weightmaps_formulas.push_back("1");
                }else{
                    if(is_verbose)std::cout<<otag<<" with associated WeightMaps Variation formula : "<<w_formula<<std::endl;
                    weightmaps_formulas.push_back(std::string(w_formula));
                }

                if(w_use== NULL){
                    weightmaps_uses.push_back("true");
                }else{
                    //if(is_verbose)std::cout<<otag<<" Loading WeightMaps Variation BlackList/WhiteList : "<<w_use<<std::endl;
                    weightmaps_uses.push_back(std::string(w_use));
                }

                if(w_mode== NULL){
                    if(is_verbose){
                        std::cout<<otag<<" No mode passed for this variation in WeightMaps'"<<std::endl;
                        std::cout<<otag<<" Assuming its the default multisim;"<<std::endl;
                    }
                    weightmaps_mode.push_back("multisim");
                }else{
                    if(is_verbose)std::cout<<otag<<" Loading WeightMaps Mode  : "<<w_mode<<std::endl;
                    std::string mode = std::string(w_mode);
                    if(mode=="multisim" || mode=="minmax"){
                        weightmaps_mode.push_back(mode);
                    }else{
                        std::cout<<otag<<" ERROR! The mode passed in is "<<mode<<" but only allowed is multisim or minmax.'"<<std::endl;
                        exit(EXIT_FAILURE);
                    }
                }

                pVariation = pVariation->NextSiblingElement("variation");
            }

            pWeiMaps=pWeiMaps->NextSiblingElement("WeightMaps");
        }
    }


    if(!pShapeOnlyMap){
        if(is_verbose){
            std::cout<<otag<<" Not setting up for shape-only covariance matrix generation. MAKE SURE this is what you want if you're generating covariance matrix!!!!" << std::endl;
            std::cout<<otag<<" SAFELY IGNORE above message if you're not generating covariance matrix" << std::endl; 
        }
    }else{
        while(pShapeOnlyMap){

            std::string pshapeonly_systematic_name = std::string(pShapeOnlyMap->Attribute("name"));
            const char* pshapeonly_systematic_use = pShapeOnlyMap->Attribute("use");
            bool pshapeonly_systematic_use_bool = true;

            if(pshapeonly_systematic_use == NULL || std::string(pshapeonly_systematic_use) == "true"){
                std::cout << otag << " Setting up shape-only covariance matrix for systematic: " << pshapeonly_systematic_name << std::endl;
            }else if(std::string(pshapeonly_systematic_use) == "false"){
                std::cout << otag << " Setting up shape-only covariance matrix for systematic: " << pshapeonly_systematic_name << "? False" << std::endl;
                pshapeonly_systematic_use_bool = false;
            }else{
                std::cout << otag << " INVALID argument received for Attribute use of ShapeOnlyUncertainty element for systematic: " << pshapeonly_systematic_name << std::endl;
                std::cout << otag << " Default it to true" << std::endl;
            }

            TiXmlElement *pSubchannel;
            pSubchannel = pShapeOnlyMap->FirstChildElement("subchannel");	

            while(pshapeonly_systematic_use_bool && pSubchannel){

                const char* pshapeonly_subchannel_name = pSubchannel->Attribute("name");
                const char* pshapeonly_subchannel_use = pSubchannel->Attribute("use");

                if(pshapeonly_subchannel_use && std::string(pshapeonly_subchannel_use) == "false" ){
                    std::cout << otag << " Subchannel " << std::string(pshapeonly_subchannel_name) << " is not included " << std::endl;
                }else{
                    std::cout << otag << " Subchannel " << std::string(pshapeonly_subchannel_name) << " is included " << std::endl;
                    shapeonly_listmap[pshapeonly_systematic_name].emplace_back(pshapeonly_subchannel_name);
                }

                pSubchannel = pSubchannel->NextSiblingElement("subchannel");
            }

            pShapeOnlyMap = pShapeOnlyMap->NextSiblingElement("ShapeOnlyUncertainty");
        }

        std::cout<<otag<< " Finish setting up systematics and subchannels for shape-only covariance matrix generation " << std::endl;
    }


    while(pSpec){
        const char* swrite_out = pSpec->Attribute("writeout");
        const char* swrite_out_tag = pSpec->Attribute("writeout_tag");
        const char* sform_matrix = pSpec->Attribute("form_matrix");	

        if( std::string(swrite_out) == "true") write_out_variation = true;
        if(write_out_variation){
            if(swrite_out_tag == NULL) write_out_tag="UNSET";
            else write_out_tag = std::string(swrite_out_tag);
        }

        if( std::string(sform_matrix) == "false") form_covariance = false;

        pSpec = pSpec->NextSiblingElement("varied_spectrum");
    }
    if(is_verbose)  std::cout << otag << " Mode for varied spectra (if needed):   write out: " << (write_out_variation ? "true" : "false") << ", tag: " << write_out_tag << " | form covariance matrix: " << (form_covariance ? "true" : "false") << std::endl;


    //Where is the "data" folder that keeps pre-computed spectra and rootfiles
    //Will eventuall just configure this in CMake
    while(pData){
        data_path = pData->Attribute("path");
        pData = pData->NextSiblingElement("data");
        if(is_verbose)std::cout<<otag<<"data path loaded as: "<<data_path<<std::endl;
    }


    // Where is the covariance matrix you want to use, and whats its name in the root file.
    while(pCov){
        correlation_matrix_rootfile = data_path + pCov->Attribute("file");
        correlation_matrix_name = pCov->Attribute("name");
        pCov = pCov->NextSiblingElement("covariance");
    }


    if(is_verbose) std::cout<<otag<<"|| Calculating used things"<<std::endl;

    // so num_channels here is number of TOTAL channels in xml.
    num_channels = channel_names.size();
    num_modes = mode_names.size();
    num_detectors  = detector_names.size();
    //Calculate bin_widths from bin_edges

    num_subchannels_xml = num_subchannels;
    num_channels_xml = num_channels;
    num_modes_xml = num_modes;
    num_detectors_xml = num_detectors;





    // here we run through every combination, and make note when (and where binwise) all the subchannels that are turned on are.
    std::string tempn;
    int indexcount = 0;


    for(int im = 0; im < num_modes; im++){

        for(int id =0; id < num_detectors; id++){

            for(int ic = 0; ic < num_channels; ic++){

                for(int sc = 0; sc < num_subchannels.at(ic); sc++){

                    tempn = mode_names[im] +"_" +detector_names[id]+"_"+channel_names[ic]+"_"+subchannel_names[ic][sc];
                    if(is_verbose)std::cout<<otag<<""<<tempn<<" "<<im<<" "<<id<<" "<<ic<<" "<<sc<<std::endl;

                    map_subchannel_plotnames[tempn] = subchannel_plotnames[ic][sc];

                    // This is where you choose NOT to use some fields
                    if(mode_bool[im] && detector_bool[id] && channel_bool[ic] && subchannel_bool[ic][sc]){

                        fullnames.push_back(tempn);
                        vec_is_data.push_back(subchannel_datas[ic][sc]);
                        for(int k = indexcount; k < indexcount+num_bins.at(ic); k++){
                            used_bins.push_back(k);
                            //std::cout<<"USED: "<<k<<std::endl;

                        }
                    }
                    std::vector<int> tvec = {indexcount, indexcount+num_bins.at(ic)-1};

                    map_tag_to_covariance_index[tempn] = 	tvec;
                    indexcount = indexcount + num_bins.at(ic);


                }
            }
        }
    }


    if(is_verbose) std::cout<<otag<<"There are "<< fullnames.size() << " used subchannel histograms" << std::endl;


    //For here on down everything is derivable, above is just until I actually Get config working.
    log<LOG_DEBUG>(L"%1% || Starting on MC file parameters.") % __func__;

    num_modes=0;
    for(int i=0;i<mode_bool.size(); i++){	if(mode_bool.at(i)){num_modes++; mode_used.push_back(i);}	}

    num_detectors = 0;
    for(int i=0; i<detector_bool.size(); i++){ if(detector_bool.at(i)) {num_detectors++; detector_used.push_back(i);}	}

    for(int i=0; i< num_channels; i++){
        num_subchannels.at(i) = 0;
        for(int j=0; j<subchannel_bool[i].size(); j++){ 
            if(subchannel_bool[i][j]){
                num_subchannels[i]++;
                subchannel_used[i].push_back(j);	
            }
        }
    }
    //This needs to be above num_channel recalculation;

    num_channels = 0;
    for(int i=0; i< channel_bool.size(); i++){
        if( channel_bool.at(i)){
            num_channels++;
            channel_used.push_back(i);
        }
    }



    log<LOG_DEBUG>(L"%1% || calculating total bins.") % __func__;
    this->CalcTotalBins();


    if(is_verbose){
        std::cout<<otag<<"Checking number of XX"<<std::endl;
        std::cout<<otag<<"--> num_modes: "<<num_modes<<" out of "<<num_modes_xml<<std::endl;
        std::cout<<otag<<"--> num_detectors: "<<num_detectors<<" out of "<<num_detectors_xml<<std::endl;
        std::cout<<otag<<"--> num_channels: "<<num_channels<<" out of "<<num_channels_xml<<std::endl;
        for(auto i: channel_used){
            std::cout<<otag<<"----> num_subchannels: "<<num_subchannels.at(i)<<" out of "<<num_subchannels_xml.at(i)<<std::endl;
            std::cout<<otag<<"----> num_bins: "<<num_bins.at(i)<<std::endl;
        }

        std::cout<<otag<<"--> num_bins_detector_block: "<<num_bins_detector_block<<std::endl;
        std::cout<<otag<<"--> num_bins_detector_block_compressed: "<<num_bins_detector_block_compressed<<std::endl;
        std::cout<<otag<<"--> num_bins_mode_block: "<<num_bins_mode_block<<std::endl;
        std::cout<<otag<<"--> num_bins_mode_block_compressed: "<<num_bins_mode_block_compressed<<std::endl;


        std::cout<<otag<<"--> num_bins_total: "<<num_bins_total<<std::endl;
        std::cout<<otag<<"--> num_bins_total_compressed: "<<num_bins_total_compressed<<std::endl;


    }

    //Now delete all info corresponding to NON-USED channels & subchannels.


    //num_subchannels, num_bins, bin_edges, bin_widths, channel_names, subchannel_names;
    auto temp_num_subchannels = num_subchannels;
    auto temp_num_bins= num_bins;
    auto temp_subchannel_names = subchannel_names;
    auto temp_channel_names = channel_names;
    auto temp_channel_units = channel_units;
    auto temp_bin_edges = bin_edges;
    auto temp_bin_widths = bin_widths;
    auto temp_detector_names = detector_names;
    auto temp_mode_names = mode_names;
    auto temp_mode_bool = mode_bool;
    auto temp_channel_bool= channel_bool;
    auto temp_detector_bool = detector_bool;
    auto temp_subchannel_bool = subchannel_bool;
    auto temp_subchannel_osc_patterns = subchannel_osc_patterns;
    auto temp_subchannel_plotnames = subchannel_plotnames;

    num_subchannels.clear();
    num_bins.clear();
    channel_names.clear();
    channel_units.clear();
    bin_edges.clear();
    bin_widths.clear();
    detector_names.clear();
    mode_names.clear();

    mode_bool.clear();
    channel_bool.clear();
    detector_bool.clear();

    subchannel_names.clear();
    subchannel_names.resize(num_channels);
    subchannel_bool.clear();
    subchannel_bool.resize(num_channels);
    subchannel_plotnames.clear();
    subchannel_plotnames.resize(num_channels);
    subchannel_osc_patterns.clear();
    subchannel_osc_patterns.resize(num_channels);

    int ic=0;
    for(int c: channel_used){
        //if(is_verbose){std::cout<<otag<<"Adding channel: "<<c<<std::endl;}
        num_subchannels.push_back( temp_num_subchannels.at(c));
        num_bins.push_back( temp_num_bins.at(c));
        channel_names.push_back( temp_channel_names.at(c));
        channel_units.push_back(temp_channel_units.at(c));
        bin_edges.push_back( temp_bin_edges.at(c));
        bin_widths.push_back( temp_bin_widths.at(c));

        channel_bool.push_back(temp_channel_bool.at(c));

        for(int sc: subchannel_used[c]){
            subchannel_bool[ic].push_back(temp_subchannel_bool.at(c).at(sc));
            subchannel_names[ic].push_back(temp_subchannel_names.at(c).at(sc));
            subchannel_osc_patterns[ic].push_back(temp_subchannel_osc_patterns.at(c).at(sc));
            subchannel_plotnames[ic].push_back(temp_subchannel_plotnames[c][sc]);
        }
        ic++;
    }
    for(int d: detector_used){
        detector_names.push_back(temp_detector_names.at(d));
        detector_bool.push_back(temp_detector_bool.at(d));
        //if(is_verbose) std::cout<<otag<<"Using Detector: "<<detector_names.back()<<std::endl;
    }

    for(int m: mode_used){
        mode_names.push_back(temp_mode_names.at(m));
        mode_bool.push_back(temp_mode_bool.at(m));
    }



    if(is_verbose) {

        //	std::cout<<"--> num_channels: "<<num_channels<<" channel_bool.size(): "<<channel_bool.size()<<" channel_names.size(): "<<channel_names.size()<<std::endl;
        //	std::cout<<"--> num_modes: "<<num_modes<<" mode_bool.size(): "<<mode_bool.size()<<" mode_names.size(): "<<mode_names.size()<<std::endl;
        //	std::cout<<"--> num_detectors: "<<num_detectors<<" detector_bool.size(): "<<detector_bool.size()<<" detector_names.size(): "<<detector_names.size()<<std::endl;
        for(int i=0; i< num_channels; i++){
            //		std::cout<<"--> num_subchannels: "<<num_subchannels.at(i)<<" subchannel_bool.size(): "<<subchannel_bool.at(i).size()<<" subchannel_names.at(i).size(): "<<subchannel_names.at(i).size()<<std::endl;
        }
    }



    a_num_bins = num_bins.data();
    a_num_subchannels = num_subchannels.data();



    if(is_verbose){std::cout<<otag<<"Done!"<<std::endl;}
    if(is_verbose){std::cout<<otag<<"---------------------------------------------------------------"<<std::endl;}

    return 0;
}
