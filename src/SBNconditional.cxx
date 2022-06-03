#include "SBNconditional.h"
using namespace sbn;

std::vector<TMatrixT<double>> sbn::splitCovariance(TMatrixT<double> & input, int start_cons_point){ 
//        std::cout<<"SPLOT "<<input.GetNcols()<<std::endl;
        std::cout<<"sig binning "<<0<<" "<<start_cons_point-1<<" "<<0<<" "<<start_cons_point-1<<std::endl;
        TMatrixT<double> sig_sig = input.GetSub(0,start_cons_point-1,0,start_cons_point-1); 
        TMatrixT<double> cons_cons = input.GetSub(start_cons_point, input.GetNrows()-1 ,start_cons_point,input.GetNrows()-1); 
        TMatrixT<double> sig_cons = input.GetSub(start_cons_point,input.GetNrows()-1,0,start_cons_point-1); 
        TMatrixT<double> cons_sig = input.GetSub(0,start_cons_point-1,start_cons_point, input.GetNrows()-1); 


        std::vector<TMatrixT<double>> ans = {sig_sig,sig_cons,cons_sig,cons_cons};
        return ans;
}


TMatrixT<double> sbn::getConstrainedCovariance(std::vector<TMatrixT<double>>& v_mat){
        
        TMatrixT<double> cons_invert = v_mat.back();
        cons_invert.Invert();

        TMatrixT<double> ans = v_mat.front()- v_mat[2]*cons_invert*v_mat[1];
        return ans;
}

std::vector<double> sbn::getConstrainedPrediction(std::vector<TMatrixT<double>>& v_mat, std::vector<double> pred, std::vector<double> data, int start_cons_point){

	std::vector<double> ans(pred.begin(), pred.begin()+start_cons_point);   //nue prediction

	//difference between data&MC in numu bins
	std::vector<double> numu_diff;
	for(int i=start_cons_point; i< data.size(); i++){
		numu_diff.push_back(data[i]-pred[i]);
	}

	TMatrixT<double> cons_invert = v_mat.back();
        cons_invert.Invert();
	TMatrixT<double> corr_matrix = v_mat[2]*cons_invert;

	for(int i=0; i<start_cons_point; i++){	
	    for(int k=0; k<numu_diff.size();k++){
		ans[i] += corr_matrix(i, k)*numu_diff[k];	
	    }
	}
	return ans;
}
