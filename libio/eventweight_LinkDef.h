#include "vector"
#include "map"
#include "string"

#ifdef __CINT__
#pragma link C++ class string+;
//#pragma link C++ class string::*;
#pragma link C++ class vector<float>+;
#pragma link C++ class vector<double>+;
//#pragma link C++ class vector<Bool_t>::*;
#pragma link C++ class pair<string, vector<float>>+;
#pragma link C++ class pair<string, vector<double>>+;
#pragma link C++ class map<string, vector<float> >+;
#pragma link C++ class map<string, vector<double> >+;
#endif

