#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <vector>
#include <unistd.h>
#include <getopt.h>
#include <cstring>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TString.h"
#include "TNtuple.h"
#include "TChain.h"
#include "TMath.h"
#include "TSystem.h"
#include "TMatrixT.h"
#include "TRandom.h"
#include "TError.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TGraph.h"

#include "params.h"
#include "ngrid.h"
#include "SBNconfig.h"
#include "SBNchi.h"
#include "SBNspec.h"
#include "SBNosc.h"
#include "SBNfit.h"
#include "SBNfit3pN.h"
#include "SBNcovariance.h"
#include "loghelper.h"
log_level_t GLOBAL_LEVEL = LOG_DEBUG;


#define no_argument 0
#define required_argument 1
#define optional_argument 2

using namespace sbn;

/*************************************************************
 *************************************************************
 *		BEGIN example.cxx
 ************************************************************
 ************************************************************/
int main(int argc, char* argv[])
{

    std::string xml = "example.xml";
    bool print_mode = false;

    /*************************************************************
     *************************************************************
     *		Command Line Argument Reading
     ************************************************************
     ************************************************************/
    const struct option longopts[] =
    {
        {"xml", 		required_argument, 	0, 'x'},
        {"print", 		no_argument, 		0, 'p'},
        {0,			no_argument, 		0,  0},
    };

    int iarg = 0;
    opterr=1;
    int index;

    while(iarg != -1)
    {
        iarg = getopt_long(argc,argv, "x:p", longopts, &index);

        switch(iarg)
        {
            case 'x':
                xml = optarg;
                break;
            case 'p':
                print_mode=true;
                break;

            case '?':
            case 'h':
                std::cout<<"Allowed arguments:"<<std::endl;
                std::cout<<"\t-x\t--xml\t\tInput .xml file for SBNconfig"<<std::endl;
                std::cout<<"\t-p\t--print\t\tRuns in print mode, making a lot of plots of canvases and Variations."<<std::endl;
                return 0;
        }
    }

    /*************************************************************
     *************************************************************
     *			Main Program Flow
     ************************************************************
     ************************************************************/
    time_t start_time = time(0);

    //a tag to identify outputs
    std::string tag = "EXAMPLE1";


    NGrid mygrid;
    mygrid.AddDimension("dm", -4, 4, 0.5);
    mygrid.AddDimension("sin",-4, 0, 0.25);
    mygrid.AddDimension("dcp", 2.1, 5.4, 0.631);


    std::vector<std::vector<double>> grid = mygrid.GetGrid();
    std::cout<<"Total pts: "<<mygrid.f_num_total_points<<std::endl;
    std::cout<<mygrid.f_dimensions[0].f_name<<" "<<mygrid.f_dimensions[0].f_N<<std::endl;
    std::cout<<mygrid.f_dimensions[1].f_name<<" "<<mygrid.f_dimensions[1].f_N<<std::endl;
    std::cout<<mygrid.f_dimensions[2].f_name<<" "<<mygrid.f_dimensions[2].f_N<<std::endl;

    std::cout<<"We have "<<grid.size()<< " points with : "<<grid[0].size()<<" dimensions"<<std::endl;

    for(int i=0; i<grid.size(); i++){
        std::cout<<i<<" (";
        for(int d=0; d< grid[i].size(); d++){
            std::cout<<grid[i][d]<<" ";
        }
        std::cout<<")"<<std::endl;
    }


    return 0;


    std::cout<<"Begining Covariance Calculation for tag: "<<tag<<std::endl;

    //Create a SBNcovariance object initilizing with the inputted xml
    //This will load all the files and weights as laid out
    SBNcovariance example_covar(xml);

    //Form the covariance matrix from loaded weights and MC events
    example_covar.FormCovarianceMatrix(tag);

    if(print_mode){
        //and make some plots of the resulting things
        //Will be outputted in the form: SBNfit_covariance_plots_TAG.root
        example_covar.PrintMatricies(tag);

        //Will be outputted in the form: SBNfit_variation_plots_TAG.root
        example_covar.PrintVariations(tag);
    }

    std::cout << "Total wall time: " << difftime(time(0), start_time)/60.0 << " Minutes.\n";
    return 0;



    int nb = 4;
    int nc = 1;
    TMatrixD inv_C(nb,nb);
    inv_C.Zero();
    inv_C(0,0)=6400+400;
    inv_C(0,1)=9600;
    inv_C(0,2)=912000;
    inv_C(0,3)=456000;
    inv_C(1,0)=9600;
    inv_C(1,1)=14400+600;
    inv_C(1,2)=1368000;
    inv_C(1,3)=684000;
    inv_C(2,0)=912000;
    inv_C(2,1)=1368000;
    inv_C(2,2)=144000000+60000;
    inv_C(2,3)=72000000;
    inv_C(3,0)=456000;
    inv_C(3,1)=684000;
    inv_C(3,2)=72000000;
    inv_C(3,3)=26000000+30000;
   
    inv_C.Print();
    inv_C.Invert();
    
    float mc[4] ={400,600,60000,30000};
    float data[4] ={400,600,72000,36000};

    SBNconstraint cc(inv_C,mc,data,nb,nc);
    cc.print();


}
