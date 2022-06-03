#include "SBNfit3pN.h"
using namespace sbn;

/****************************************************************
 *		Generic 3+N 
 * *************************************************************/


SBNfit3pN::SBNfit3pN(SBNosc inBk, SBNosc inSg, int npa) : SBNfit(inBk,inSg,npa), signal_osc_spectrum(inSg) {

}


double SBNfit3pN::MinimizerCalcChi(const double * X){
	num_func_calls++;
	SBNosc tempOsc = signal_osc_spectrum; 

		double imn[3] = {X[0],X[1],X[2]};
		double iue[3] = {X[3],X[4],X[5]};
		double ium[3] = {X[6],X[7],X[8]};
		double iph[3] = {X[9],X[10],X[11]};
		NeutrinoModel signalModel(imn,iue,ium,iph);
					
		tempOsc.LoadModel(signalModel);
		std::vector<std::vector<double>> ans = tempOsc.Oscillate("tag");
	

		last_calculated_chi =this->CalcChi(ans[0]);
	return last_calculated_chi;

}

/****************************************************************
 *			3+1	Only 
 * *************************************************************/

SBNfit3p1::SBNfit3p1(SBNosc inBk, SBNosc inSg, int npa) : SBNfit(inBk,inSg,npa), signal_osc_spectrum(inSg) {

}

double SBNfit3p1::MinimizerCalcChi(const double * X){
	num_func_calls++;
	SBNosc tempOsc = signal_osc_spectrum; 
	
	NeutrinoModel signalModel(X[0],X[1],X[2]);
	signalModel.numsterile=1;		

	tempOsc.LoadModel(signalModel);

	std::vector<std::vector<double>> ans = tempOsc.Oscillate("tag");
	
	last_calculated_chi =this->CalcChi(ans[0]);

	return last_calculated_chi;
}

