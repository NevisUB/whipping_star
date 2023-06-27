#include "denan.h"
#include "denan_LinkDef.h"
 
double DeNan(double x, double replace){
    if(x!=x || isinf(x))return replace;
    return x;
}



