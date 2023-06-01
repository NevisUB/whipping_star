#ifndef NGRID_H_
#define NGRID_H_

#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <map>
#include <time.h>
#include <stdexcept>
#include <cstdlib>

#include "TMatrixT.h"

//--- local header ----
#include "tinyxml.h"

struct SBNfitResult{

    size_t bf_pt;
    double chi_min;
    std::vector<double> chi_vec;
    TMatrixT<double> bf_covariance_matrix;
    TMatrixT<double> bf_inverse_covariance_matrix;

    SBNfitResult(std::vector<double> in_chi_vec, size_t in_bf_pt, double in_chi_min, TMatrixT<double> in_bf_covariance_matrix, TMatrixT<double> in_bf_inverse_covariance_matrix) : chi_vec(in_chi_vec), bf_pt(in_bf_pt), chi_min(in_chi_min),  bf_covariance_matrix(in_bf_covariance_matrix) ,  bf_inverse_covariance_matrix(in_bf_inverse_covariance_matrix) {};

};



struct NGridDimension{
    std::string f_name;
    double f_min;
    double f_max;
    double f_step;
    double f_fixed_value;
    bool f_is_fixed;
    int f_N;
	
    bool f_has_constrain;  
    double f_constrain_value;
    double f_constrain_range;

    std::vector<double> f_points;
    std::vector<double> f_edges;  //edges that can be used in creating histograms

    NGridDimension(std::string name, double min, double max, double step) : f_name(name), f_min(min), f_max(max), f_step(step) {
        f_N = ceil(fabs(f_min-f_max)/step);
        f_points.resize(f_N);
        this->CalcGrid();
        f_is_fixed = false;
	f_has_constrain = false;
    };


    NGridDimension(std::string name, double min, double max, double step, double constrained_value, double constrain_range):f_name(name), f_min(min), f_max(max), f_step(step), f_constrain_value(constrained_value), f_constrain_range(constrain_range){
	f_N = ceil(fabs(f_min-f_max)/step);
	f_points.resize(f_N);
	this->CalcGrid();
	f_is_fixed=false;
	f_has_constrain=true;
    }


    NGridDimension(std::string name, double min, double max, double step, double constrained_value):f_name(name), f_min(min), f_max(max), f_step(step), f_constrain_value(constrained_value){
	f_N = ceil(fabs(f_min-f_max)/step);
	f_points.resize(f_N);
	this->CalcGrid();
	f_is_fixed=false;
	f_has_constrain=true;
	f_constrain_range=-999;
    }


    NGridDimension(std::string name, double val) : f_name(name), f_fixed_value(val), f_is_fixed(true){
        f_N = 1;
        f_step = 0.0;
        f_max = f_fixed_value;
        f_min = f_fixed_value;
        f_points.resize(f_N);
        this->CalcGrid();
	f_has_constrain = false;
    }


    int GetNPoints() const {return f_N;};

    void CalcGrid(){
        for(int i=0; i<f_N; i++){
            f_points[i]= f_min + i*f_step;
        }
	f_edges = f_points;
	f_edges.push_back(f_points.back() + f_step);
        return;
    }

    inline const std::string& GetName() const {return f_name; }
    inline double GetPoint(int n){ return f_points[n];    };
    inline const std::vector<double>& GetPoints() const { return f_points; }
    inline const std::vector<double>& GetEdges() const {return f_edges; }
    inline double GetMax() const {return f_max;}
    inline double GetMin() const {return f_min;}
    inline double GetStep() const { return f_step;}
};

struct NGrid{
    int f_num_dimensions = 0;
    int f_num_total_points = 1;

    std::vector<NGridDimension> f_dimensions;

    //constructors 
    NGrid(){}

    // given provided configuration xml, grab the grids in it
    NGrid(std::string xml){

	std::cout << "Load Grid using configuration xml: " << xml << std::endl;

	TiXmlDocument doc(xml.c_str());
	bool loadOkay = doc.LoadFile();
        if(!loadOkay){
	    throw std::runtime_error("\tERROR: fail to load configuration XML, generally means broken .xml brackets or attribute syntax.");
        }
        TiXmlHandle hDoc(&doc);
	TiXmlElement *pGrid;
	pGrid = doc.FirstChildElement("GridDimension");
	
	while(pGrid){
	    const char* dim_name = pGrid->Attribute("name");
	    const char* dim_fixval = pGrid->Attribute("fix_val");
	    const char* dim_lowedge = pGrid->Attribute("min");
	    const char* dim_upedge = pGrid->Attribute("max");
	    const char* dim_stepsize = pGrid->Attribute("step");
	    char* pEnd; 

            if(dim_name==NULL){
                std::cerr << "\tERROR! Grid dimension need a name! Please kindly define a name for it"<<std::endl;
		print_config_rule();
		throw std::runtime_error("Invalid Syntax to setup Grid");
            }	
	    std::cout << "--> On Grid dimension " << dim_name;

	    if(dim_fixval){
		std::cout << " with fixed value : " << dim_fixval << std::endl;
	 	this->AddFixedDimension(std::string(dim_name), strtod(dim_fixval, &pEnd));	
	    }else{

		if(dim_lowedge == NULL || dim_upedge == NULL || dim_stepsize == NULL){
		    print_config_rule();
		    throw std::runtime_error("Invalid Syntax to setup Grid");
		}

		std::cout << " with min: " << dim_lowedge << " max: " << dim_upedge << " and step size: " << dim_stepsize << std::endl;
		this->AddDimension(std::string(dim_name), strtod(dim_lowedge, &pEnd),  strtod(dim_upedge, &pEnd),  strtod(dim_stepsize, &pEnd));
	    }

            pGrid = pGrid->NextSiblingElement("GridDimension");
        }

	std::cout << "Finish reading grids from xml... Total of " << f_num_dimensions << " dimensions."  << std::endl;

    }

    void AddDimension(std::string name, double min, double max, double step ){
        f_dimensions.emplace_back( NGridDimension(name,min,max,step));
        f_num_dimensions++;
        f_num_total_points *= f_dimensions.back().GetNPoints();
        return;
    }

    void AddDimension(std::string name, std::string grid_scan){
           std::vector<double> vect;
           std::stringstream ss(grid_scan);

            double number;
            while ( ss >> number ) vect.push_back( number );

 
        double min = vect[0];
        double max = vect[1];
        double step = fabs(max-min)/vect[2];
        if(min>=max){
            std::cout<<"ERROR! min grid value ("<<min<<")  is larger than max ("<<max<<") in grid string "<<grid_scan<<std::endl;
            exit(EXIT_FAILURE);
        }
        std::cout<<"NGrid definied with a min value of "<<min<<" a max value of "<<max<<" with "<<vect[2]<<" steps of size "<<step<<std::endl;

        f_dimensions.emplace_back( NGridDimension(name,min,max,step));
        f_num_dimensions++;
        f_num_total_points *= f_dimensions.back().GetNPoints();
        return;
    }

    void AddFixedDimension(std::string name, double val){
        f_dimensions.emplace_back(NGridDimension(name,val));
        f_num_dimensions++;
        return;
    }


    void AddConstrainedDimension(std::string name, double min, double max, double step, double constrained_value, double constrain_range){
	f_dimensions.emplace_back(NGridDimension(name, min, max, step, constrained_value, constrain_range));
	f_num_dimensions++;
	f_num_total_points *= f_dimensions.back().GetNPoints();
	return ;
    }

    void AddConstrainedDimension(std::string name, double min, double max, double step, double constrained_value){
	f_dimensions.emplace_back(NGridDimension(name, min, max, step, constrained_value));
	f_num_dimensions++;
	f_num_total_points *= f_dimensions.back().GetNPoints();
	return ;
    }




    // get a vector of Ndim grid points
    // NGrid added later will be varied first.
    std::vector<std::vector<double>> GetGrid(){
        std::vector<std::vector<double>> grid;            

        //count from 0 to f_num_total_points
        for(int i=0; i<f_num_total_points; i++){

            std::vector<double> point(f_num_dimensions,-99);

            //were going to take each number and write each digit in base K where K is that vectors length
            int divisor=1;
            for(int j =f_num_dimensions-1 ;j>=0; j--){

                int this_index =  (i/divisor)%f_dimensions[j].GetNPoints();

                point[j] = f_dimensions[j].GetPoint(this_index);

                //in order so that each digit is written in the correct base, modify divisor here
                divisor=divisor*f_dimensions[j].GetNPoints();
            }

            grid.push_back(point);
        }



        return grid;
    }



    int Print(){

        std::vector<std::vector<double>> grid = this->GetGrid();
        std::cout<<"Total pts: "<<this->f_num_total_points<<std::endl;
        for(int i=0; i< f_num_dimensions; i++){
            if(this->f_dimensions[i].f_is_fixed){
                std::cout<<this->f_dimensions[i].f_name<<" (fixed) "<<this->f_dimensions[i].f_N<<std::endl;
            }else{
                std::cout<<this->f_dimensions[i].f_name<<" "<<this->f_dimensions[i].f_N<<std::endl;
            }
        }
        std::cout<<"We have "<<grid.size()<< " points with : "<<grid[0].size()<<" dimensions"<<std::endl;
       
        /*
        for(int i=0; i< grid.size(); i++){
            std::cout<<i;
            for(int j=0; j< grid[i].size(); j++){
                std::cout<<" "<<grid[i][j];
            }
            std::cout<<std::endl;

        }
        */

        return 0;

    }

    void print_config_rule(){
	std::cerr << "\n"<< std::string(30, '_') <<  std::endl;
	std::cerr << "\tRULEs for xml-based grid setup" << std::endl;
	std::cerr << "\tPlease specify the name, lower, upper range for the grid, as well as the stepsize. Please follow style below: " << std::endl; 
	std::cerr << "\t\tOption 1: name=\"xx\" min=\"xx\" max=\"xx\" step=\"xx\" to setup evenly distributed grid" << std::endl;
	std::cerr << "\t\tOption 2: name=\"xx\" fix_val=\"xx\" to setup fixed value for this dimension" << std::endl;
	std::cerr << "\n"<< std::string(30, '_') <<  std::endl;
 	return;	
    }
};

#endif
