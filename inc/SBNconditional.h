#ifndef SBNCONDO_H_
#define SBNCONDO_H_

#include <cmath>
#include <vector>
#include <iostream>
#include "SBNspec.h"
#include "SBNchi.h"
#include "SBNconfig.h"

#include "TH1.h"
#include "TH2.h"
#include "TMatrixT.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TStyle.h"
#include "TLine.h"
#include "TLatex.h"
#include "TText.h"
#include "TMath.h"
#include "TGraph.h"

#include "TMath.h"
#include <ctime>
#include "params.h"

#include "TDecompChol.h"
#include "TDecompSVD.h"
#include "TMatrixDEigen.h"
#include "TMatrixDSymEigen.h"

#include "ngrid.h"
#include <gsl/gsl_randist.h>

namespace sbn{

    std::vector<TMatrixT<double>> splitCovariance(TMatrixT<double> & input, int start_cons_point);
    TMatrixT<double> getConstrainedCovariance(std::vector<TMatrixT<double>>& v_mat);
    std::vector<double> getConstrainedPrediction(std::vector<TMatrixT<double>>& v_mat, std::vector<double> pred, std::vector<double> data, int start_cons_point);

}
#endif
