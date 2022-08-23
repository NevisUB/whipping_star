#ifndef HIVEP_H_
#define HIVEP_H_

#include <cmath>
#include <vector>
#include <iostream>
#include "SBNconfig.h"
#include <TH1D.h>
#include <string>
#include <TF1.h>
#include <THStack.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TLine.h>
#include <TStyle.h>
#include <TMatrixT.h>
//#include <TROOT.h>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <numeric>

#include <ctime>
#include <TFile.h>
#include "params.h"
#include <TRandom3.h>

#include <vector>
#include <string>
#include <set>
#include <numeric>
#include <iostream>
#include <algorithm>
#include <numeric>

#include "TTreeFormula.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TPad.h"
#include "TH1.h"
#include "TH2.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLatex.h"
#include "THStack.h"
#include "TGraph.h"
#include "TLine.h"
#include "TMVA/Types.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Reader.h"
#include "TFriendElement.h"
#include "TText.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TVectorT.h"
#include "TEntryList.h"
#include "TObject.h"

int hivePlotStack(std::vector<TH1D> vec_th1 , TH1 * tsum, TMatrixD * covar_collapsed, std::string tag, std::string unit,std::vector<TColor*> cols, std::vector<std::string> &names, std::vector<int> fillstyles );

TText * drawPrelim(double x, double y);
TText * drawPrelim(double x, double y,double s);
TText * drawPrelim(double x, double y,double s, std::string in);
TText * drawPrelim(double x, double y, std::string in);


#endif
