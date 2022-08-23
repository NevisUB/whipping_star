#ifndef SBNprob_H_
#define SBNprob_H_

#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <array>
#include <map>
//#include <globes/globes.h>   /* GLoBES library */

#include <prob.h>


namespace sbn{

	struct SBNprob{
		public:	
			SBNprob(int);
			SBNprob(int, std::vector<double>,std::vector<double>, std::vector<double>);

			int init();

			double t12,t13,t23,t14,t24,t34;
			double dm21, dm31,dm41;
			double ms1, ms2,ms3,ms4;
			double d13,d24,d34;

			double potential_cc,potential_nc;
			bool use_matter_effect;
			bool use_nc_matter_effect;
			bool use_antineutrino_mode;

			int dimension;
			double convert_ev_to_gev;	
			double convert_km_to_inverse_gev;
			double rho;
			double degree;
			double num_neutrinos;
			double conversion_parameter;

			double GaussianPDF(double x, double mean, double sigma);

			ComplexMatrix hamiltonian;
			ComplexMatrix hamil_kin;
			ComplexMatrix potential;
			ComplexMatrix utvu;
			ComplexMatrix U;
			ComplexMatrix u_conj;	
			
			//New methodology
			ComplexMatrix u_h0_ut;

			int SetMatterEffect(bool);
			int SetAntiNeutrinoMode(bool);
			int SetNCMatterEffect(bool);
			int SetParameters( std::vector<double>,std::vector<double>, std::vector<double>);

			double ProbabilityVacuumExact(int a, int b ,double E, double L);
			double ProbabilityVacuumExact(int a, int b, int nuornubar, double E, double L );
			
			double ProbabilityMatterExact(int a, int b ,int nuornubar, double E, double L);
			double ProbabilityMatterExact(int a, int b ,double E, double L);

			//double ProbabilityGlobes(int a, int b, int panti, double E, double L );

			double ProbabilityMatterExactSmear(int, int ,double, double, double p, double n);

			int PlotProbabilityMatter(int a, int b, double Emin, double Emax, double L, double percen, double n, std::ofstream *filestream);

			std::vector<double> GetTernaryPoints(int start_flavour, double Energy, double Length);


	};



}

#endif
