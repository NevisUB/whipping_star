#include "SBNprob.h"
using namespace sbn;

SBNprob::SBNprob(int dim,std::vector<double> angles, std::vector<double> phases, std::vector<double> mass ) :hamiltonian(dim), hamil_kin(dim), potential(dim), utvu(dim), U(dim), u_conj(dim),u_h0_ut(dim){
	dimension=dim;
	num_neutrinos = dim;
	degree = 3.14159/180.0;
	conversion_parameter = 5.06842*1e9;//this is frim km to inv ev
	
	this->SetParameters(angles,phases,mass);
	rho = 2.8;

	use_matter_effect = true;
	use_nc_matter_effect = true;
	use_antineutrino_mode = false;
	this->init();
}
int SBNprob::SetParameters(std::vector<double> angles, std::vector<double> phases, std::vector<double> mass){

	t12=angles.at(0)*degree;
	t23=angles.at(1)*degree;
	t13=angles.at(2)*degree;
	
	t14=angles.at(3)*degree;
	t24=angles.at(4)*degree;
	t34=angles.at(5)*degree;

	d13=phases.at(0)*degree;
	d24=phases.at(1)*degree;
	d34=phases.at(2)*degree;

	dm21=mass.at(0);
	dm31=mass.at(1);
	dm41=mass.at(2);

	ms1 = fabs(dm21);
	ms2 = fabs(dm21)+dm21;
	ms3 = fabs(dm21)+dm31;
	ms4 = dm41; 

	this->init();

		return 0;
}

int SBNprob::init(){

	//convert_ev_to_gev = 1e-9;	
	
	double Ye=0.4957;
	//formula from KOPP theisis
	potential_cc= 7.56e-14*rho*Ye;// want potential to be in this unitless way for now;


	potential_nc= -0.5*potential_cc;

	if(!use_nc_matter_effect){
		potential_nc =0.0;
	}

	ComplexMatrix R12(num_neutrinos);
	ComplexMatrix R13(num_neutrinos);
	ComplexMatrix R23(num_neutrinos);
	ComplexMatrix R14(num_neutrinos);
	ComplexMatrix R24(num_neutrinos);
	ComplexMatrix R34(num_neutrinos);

	R12.SetRotation(1,2,t12);
	R13.SetComplexRotation(1,3,t13,d13);
	R23.SetRotation(2,3,t23);

	R14.SetRotation(1,4,t14);
	R24.SetComplexRotation(2,4,t24,d24);
	R34.SetComplexRotation(3,4,t34,d34);


	U.SetIdentity();
	U.mult(&R34);
	U.mult(&R24);
	U.mult(&R23);
	U.mult(&R14);
	U.mult(&R13);
	U.mult(&R12);


	//std::cout<<"PRINT U"<<std::endl;

	for(int i=0; i<4; i++){
	//std::cout<<"("<<U.real(i,0)<<" "<<U.imag(i,0)<<") ("<<U.real(i,1)<<" "<<U.imag(i,1)<<") ("<<U.real(i,2)<<" "<<U.imag(i,2)<<") ("<<U.real(i,3)<<" "<<U.imag(i,3)<<")"<<std::endl; 

	}

	u_conj=U;
	u_conj.HermitianConjugate();

	std::vector<double> masses = {0,0.5*dm21,0.5*dm31,0.5*dm41};
	//std::cout<<"Msq "<<0<<" "<<dm21<<" "<<dm31<<" "<<dm41<<std::endl;
	hamil_kin.SetDiagonal(masses);

	//So we have U, Udagger and H

	u_h0_ut = U;
	u_h0_ut.mult(&hamil_kin);
	u_h0_ut.mult(&u_conj);

return 0;
}

double SBNprob::ProbabilityMatterExact(int a, int b, double E, double L){
	double ans;
	if(use_antineutrino_mode){
		ans = ProbabilityMatterExact(a,b,-1,E,L);
	}else{
		ans = ProbabilityMatterExact(a,b,1,E,L);
	}
	return ans; 
}

double SBNprob::ProbabilityMatterExact(int a, int b, int nuornubar, double Energy, double Length ){
	
	for(int i=0;i<4; i++){
	for(int j=0;j<4; j++){
	//	std::cout<<"UHoTUtr "<<i<<" "<<j<<" real "<<u_h0_ut.real(i,j)<<" imag "<<u_h0_ut.imag(i,j)<<std::endl;
	}
	}

	double E = Energy*1e9; //in eV
	double L = Length*conversion_parameter;// in eV^-1 also 
	
	//std::cout<<"LCONVR "<<L<<std::endl;
	
	ComplexMatrix S(num_neutrinos);

	hamiltonian = u_h0_ut;

	//Are we working with antineutrinos here?
	if(nuornubar<0){
		hamiltonian.conj();
	
		potential.real(0,0)= -1.0*potential_cc;
		potential.real(dimension-1,dimension-1)=-potential_nc;
	}else{
		potential.real(0,0)= potential_cc;
		potential.real(dimension-1,dimension-1)=-potential_nc;
	}

	hamiltonian.mult(1.0/E,1.0/E);


	
	/*for(int i=0;i<4; i++){
	for(int j=0;j<4; j++){
		std::cout<<"Hamil "<<i<<" "<<j<<" real "<<hamiltonian.real(i,j)<<" imag "<<hamiltonian.imag(i,j)<<std::endl;
	}
	}*/

	if(use_matter_effect){
		hamiltonian.add(&potential);
	}
	//std::cout<<"Hamil after adding V, "<<hamiltonian.real(0,0)<<" "<<hamiltonian.imag(0,0)<<std::endl;
	
	//Blarg, hermitian->antihermitian...  Using the fact Exp[-I M] = Cos[M]-I Sin[M], even for matricies
	std::vector<double> eigenval;
	ComplexMatrix eigenvec(num_neutrinos);
	ComplexMatrix eigenvecTr(num_neutrinos);




	//hamiltonian has units of inverse km at this point I believe. 
	//Going to find the eigenvalues and eigenvectors of this hamiltonian. This destorys temp4eigen at the moment.
	ComplexMatrix temp4eigen(num_neutrinos);
	temp4eigen = hamiltonian;
	temp4eigen.GetEigenStuff(&eigenval, &eigenvec);

	//std::cout<<"Eigen: "<<eigenval.at(0)<<" "<<eigenval.at(1)<<" "<<eigenval.at(2)<<" "<<eigenval.at(3)<<std::endl;
	
	eigenvecTr = eigenvec;
	eigenvecTr.HermitianConjugate();


	//Now calculate the S matrix in the mass basis in matter. Its diagonal here by definitoon;
	//its zero frm its constructer already
	for(int i=0; i<S.dimension; i++){
		double cphase = -L*eigenval.at(i);
		S.real(i,i) = cos(cphase); 
		S.imag(i,i) = sin(cphase); 
	}

	//Now lets transform back to the flavour basis using Q:=eigenvec
	//T0 is just a multplicative aid
	ComplexMatrix T0(num_neutrinos);
	T0 = eigenvec;
	T0.mult(&S);
	T0.mult(&eigenvecTr);

	S=T0;

	double re = S.real(abs(a),abs(b));
	double im = S.imag(abs(a),abs(b));

	return re*re+im*im;


};




double SBNprob::ProbabilityVacuumExact(int a, int b, double E, double L ){
	return ProbabilityVacuumExact(a,b,1,E,L);
}

double SBNprob::ProbabilityVacuumExact(int a, int b, int nuornubar, double E, double L ){
	use_matter_effect =false;
	double ans = ProbabilityMatterExact(a,b,nuornubar,E,L);
	use_matter_effect = true;
	return ans;
}



double SBNprob::ProbabilityMatterExactSmear(int a, int b, double E, double L ,double percen, double n){



	double sigma = percen*E/sqrt(E)+0.05*E;

	double low = E-4*sigma;
	if (low<0) low=0.01;

	double high = E+4*sigma;

	double step = fabs(low-high)/n;

	//std::cout<<E<<" low: "<<low<<" high: "<<high<<" step: "<<step<<" sigma: "<<sigma<<std::endl;
	double avg_prob = 0.0;
	for(double et = low; et<=high; et+=step){
		double tmp1 =	ProbabilityMatterExact(a,b,et,L);
		double tmp2 =	GaussianPDF(et, E, sigma);
		double tmp3 =	tmp1*tmp2*step;
		//std::cout<<E<<" "<<et<<" "<<tmp1<<" "<<tmp2<<" "<<tmp3<<" "<<avg_prob<<std::endl;
		avg_prob += tmp3;	
	}

	return avg_prob;

};


double SBNprob::GaussianPDF(double x, double mean, double sigma){
	double pi=3.14159;
	return 1.0/(sqrt(2*pi)*sigma)*exp( -pow(x-mean,2.0)/(2.0*sigma*sigma));


}

int SBNprob::PlotProbabilityMatter(int a, int b, double EminT, double EmaxT, double L, double percen, double n, std::ofstream *filestream){

	std::cout<<"Starting "<<a<<" "<<b<<" plot. L: "<<L<<" VNC: "<<potential_nc<<" VCC: "<<potential_cc<<" useMatter? "<<use_matter_effect<<" Neutrino or Anti?: "<<use_antineutrino_mode<<std::endl;

	double Emin = 0.5*EminT;//(1-0.5*percen)*EminT;
	double Emax = 1.5*EmaxT;//(1.5*percen)*EmaxT;

	int ua=abs(a);
	int ub=abs(b);

	std::vector<double> prob_vec_exact;
	std::vector<double> E_vec;
	std::vector<double> step_vec;
	double step = fabs(Emin-Emax)/n;

	for(double ee=Emin; ee<=Emax; ee+=step){
		double Ee = pow(10,ee);
		prob_vec_exact.push_back( ProbabilityMatterExact(ua,ub,Ee,L));
		E_vec.push_back(Ee);
	}

	step_vec.push_back(  fabs( E_vec.at(0)-E_vec.at(1)));
	for(int i=1; i< E_vec.size(); i++){
		step_vec.push_back( fabs( E_vec.at(i)-E_vec.at(i-1)));
	}

	for(int i=0; i< prob_vec_exact.size(); i++){
		if(E_vec.at(i) < pow(10,EminT) || E_vec.at(i)> pow(10,EmaxT)) continue;

		double avg_prob = 0.0;
		double sigma = 	percen*E_vec.at(i)+0.05;///sqrt(E);

		for(int j=0; j< prob_vec_exact.size(); j++){
			
			double tmp1 =   prob_vec_exact.at(j);
			double tmp2 =	GaussianPDF(E_vec.at(j), E_vec.at(i), sigma);
			avg_prob += tmp1*tmp2*step_vec.at(j);	
		//	std::cout<<E_vec.at(i)<<" "<<prob_vec_exact.at(i)<<" "<<E_vec.at(j)<<" "<<prob_vec_exact.at(j)<<" "<<sigma<<" "<<tmp1<<" "<<tmp2<<" "<<avg_prob<<std::endl;
		}
		*filestream<<a<<b<<" "<<L<<" "<<E_vec.at(i)<<" "<<prob_vec_exact.at(i)<<" "<<avg_prob<<std::endl;
	}


	return 0;
}


/*
double SBNprob::ProbabilityGlobes(int a, int b, int panti, double E, double L ){
	double ans = 0;
	{
	char* dunno;
	glbInit(dunno);


	glb_params globes_params = glbAllocParams();
	glbDefineParams(globes_params,t12,t13, t23, d13,dm21,dm31);
	glbSetDensityParams(globes_params,1.0,GLB_ALL);
	glbSetOscillationParameters(globes_params);

  	glbSetRates();
  
	ans =  glbConstantDensityProbability(a+1,b+1,panti, E,L, rho);

	glbFreeParams(globes_params);

	}
	return ans;

}
*/

int SBNprob::SetMatterEffect(bool in){
	use_matter_effect = in;
	this->init();
	return 0;
}

int SBNprob::SetNCMatterEffect(bool in){
	use_nc_matter_effect = in;
	this->init();
	return 0;
}

int SBNprob::SetAntiNeutrinoMode(bool in){
       use_antineutrino_mode = in;
       //d13 = -d13;
       //d24 = -d24;
       //d34 = -d34;
       this->init();
       return 0;

}

std::vector<double> SBNprob::GetTernaryPoints(int start_flavour, double Energy, double Length){

	double p_nue = this->ProbabilityMatterExact(start_flavour, 1, 1, Energy,Length);
	double p_numu = this->ProbabilityMatterExact(start_flavour, 2, 1, Energy,Length);
	double p_nutau = this->ProbabilityMatterExact(start_flavour, 3, 1, Energy,Length);

	std::vector<double> ans = { 0.5*(2*p_numu+p_nutau)/(p_nue+p_numu+p_nutau), sqrt(3)/2*p_nutau/(p_nue+p_numu+p_nutau)};

	return ans;
	
}





SBNprob::SBNprob(int dim) : hamiltonian(dim), hamil_kin(dim), potential(dim), utvu(dim), U(dim), u_conj(dim), u_h0_ut(dim) {
	dimension=dim;
	num_neutrinos = dim;
	degree = 3.14159/180.0;
	conversion_parameter = 5.06842;

	t12 = 30*degree;
	t23 = 44*degree;
	t13 = 8*degree;
	t14 = 15*degree;
	t24 = 10*degree;
	t34 = 20*degree;

	d34=0;
	d24=0;
	d13=0;

	dm21=7.5*pow(10,-5);//7.5e-5;
	dm31=2.552*pow(10,-3);//e-3;
	dm41=1;

	rho = 2.8;

	use_matter_effect = true;
	use_nc_matter_effect = true;
	this->init();
}


