#include "SBNchi.h"

#ifdef USE_GPU
#include "openacc.h"
#include "openacc_curand.h"
#endif
using namespace sbn;


std::vector<TMatrixT<double>> splitNormShape(TMatrixT<double> & Min,std::vector<double> & fullvec){
    int ncol = Min.GetNrows();
    std::vector<TMatrixT<double>> ans;
    for(int i=0; i<3;i++){
        ans.push_back(TMatrixT<double>(ncol,ncol));
        ans.back().Zero();
    }

    // ans[0] is shape, ans[1] is mixed, ans[2] is norm

    double Nt = std::accumulate(fullvec.begin(), fullvec.end(),0.0);

    for(int i=0; i<ncol; i++){
        for(int j=0; j<ncol; j++){

            ans[0](i,j)  = Min(i,j);
            ans[1](i,j)  = 0.0;
            ans[2](i,j)  = 0.0;
            for(int k=0; k<ncol;k++){
                ans[0](i,j)-fullvec[j]/Nt*Min(i,k) - fullvec[i]/Nt*Min(k,j);
                ans[1](i,j)+fullvec[j]/Nt*Min(i,k) + fullvec[i]/Nt*Min(k,j);

                for(int l=0; l<ncol; l++){
                    ans[0](i,j) += fullvec[i]*fullvec[j]/(Nt*Nt)*Min(k,l);
                    ans[1](i,j) += -2*fullvec[i]*fullvec[j]/(Nt*Nt)*Min(k,l);
                    ans[2](i,j) += fullvec[i]*fullvec[j]/(Nt*Nt)*Min(k,l);
                }

            };
        }
    }

    return ans;
}


/***********************************************
 *		Constructors
 * ********************************************/



SBNchi::SBNchi(std::string xml) : SBNconfig(xml,false){};
SBNchi::SBNchi(SBNspec in, TMatrixT<double> matrix_systematicsin) : SBNchi(in,matrix_systematicsin,false){}
SBNchi::SBNchi(SBNspec in, TMatrixT<double> matrix_systematicsin, bool is_verbose) : SBNchi(in,matrix_systematicsin,in.xmlname, is_verbose){}
SBNchi::SBNchi(SBNspec in, TMatrixT<double> matrix_systematicsin, std::string inxml, bool is_verbose) : SBNchi(in, matrix_systematicsin, inxml,  is_verbose,-1){}
SBNchi::SBNchi(SBNspec in, TMatrixT<double> matrix_systematicsin, std::string inxml, bool is_verbose, double random_seed) : SBNconfig(inxml, is_verbose), core_spectrum(in){

    last_calculated_chi = -9999999;
    is_stat_only= false;


    pseudo_from_collapsed = false;
    matrix_collapsed.ResizeTo(num_bins_total_compressed, num_bins_total_compressed);
    matrix_systematics.ResizeTo(num_bins_total, num_bins_total);
    matrix_fractional_covariance.ResizeTo(num_bins_total, num_bins_total);

    TMatrixD m = matrix_systematicsin;
    for (int i = 0; i < used_bins.size(); i++){
        TMatrixDColumn(m,i) = TMatrixDColumn(m,used_bins.at(i));
        m.ResizeTo(used_bins.size(),matrix_systematicsin.GetNcols());
    }

    m_cmin = -999;
    m_cmax = -999;

    matrix_fractional_covariance = m;
    matrix_systematics.Zero();
    max_sample_chi_val =150.0;
    m_tolerance = 1e-12;
    

    this->InitRandomNumberSeeds(random_seed);
    this->ReloadCoreSpectrum(&core_spectrum);
}

//Alternative constrctors
SBNchi::SBNchi(SBNspec in, std::string newxmlname) : SBNconfig(newxmlname), core_spectrum(in){
    is_stat_only = false;

    matrix_collapsed.ResizeTo(num_bins_total_compressed, num_bins_total_compressed);
    matrix_fractional_covariance.ResizeTo(num_bins_total, num_bins_total);

    if(fullnames.size() !=in.fullnames.size()){
        std::cerr<<"ERROR: SBNchi::SBNchi | Selected covariance matrix and background spectrum are different sizes!"<<std::endl;
        exit(EXIT_FAILURE);
    }else{
        for(int i=0; i< fullnames.size(); i++){
            if(fullnames[i]!=in.fullnames[i]){
                std::cerr<<"ERROR: SBNchi::SBNchi | Spectrum and Covariance matrix have different (or different order) subchannels!"<<std::endl;
                exit(EXIT_FAILURE);
            }
        }
    }
    m_cmin = -999;
    m_cmax = -999;


    m_tolerance = 1e-12;;
    pseudo_from_collapsed = false;
    max_sample_chi_val =150.0;
    matrix_fractional_covariance = FillSystematicsFromXML();
    last_calculated_chi = -9999999;
    core_spectrum.CollapseVector();
    matrix_systematics.ResizeTo(num_bins_total,num_bins_total);
    matrix_systematics.Zero();
    matrix_systematics=matrix_fractional_covariance;

    this->InitRandomNumberSeeds();

    this->ReloadCoreSpectrum(&in);
}

SBNchi::SBNchi(SBNspec in) : SBNchi(in,true){}


SBNchi::SBNchi(SBNspec in, bool is_is_stat_only): SBNconfig(in.xmlname), core_spectrum(in), is_stat_only(is_is_stat_only){
    last_calculated_chi = -9999999;

    matrix_collapsed.ResizeTo(num_bins_total_compressed, num_bins_total_compressed);
    matrix_systematics.ResizeTo(num_bins_total, num_bins_total);
    matrix_fractional_covariance.ResizeTo(num_bins_total, num_bins_total);
    m_cmin = -999;
    m_cmax = -999;

    m_tolerance = 1e-12;
    max_sample_chi_val =150.0;
    this->InitRandomNumberSeeds();

    pseudo_from_collapsed = false;

    if(is_is_stat_only){
        matrix_fractional_covariance.Zero();
        matrix_systematics.Zero();

    }else{
        matrix_fractional_covariance = FillSystematicsFromXML();
        matrix_systematics.Zero();
    }

    this->ReloadCoreSpectrum(&core_spectrum);

}


/***********************************************
 *		Rest for now
 * ********************************************/
void SBNchi::InitRandomNumberSeeds(){
    this->InitRandomNumberSeeds(-1);
}

void SBNchi::InitRandomNumberSeeds(double seed){

    if(seed<0){
        rangen_twister = new std::mt19937(random_device_seed());
        rangen_linear = new std::minstd_rand(random_device_seed());
        rangen_carry = new std::ranlux24_base(random_device_seed());
        rangen = new TRandom3(0);
    }else{
        rangen_twister = new std::mt19937(seed);
        rangen_linear = new std::minstd_rand(seed);
        rangen_carry = new std::ranlux24_base(seed);
        rangen = new TRandom3(seed);
    }

    m_dist_normal=new std::normal_distribution<float>;
    std::normal_distribution<float> dtemp(0.0,1.0);
    m_dist_normal->param(dtemp.param());
}


int SBNchi::ReloadCoreSpectrum(SBNspec *bkgin){
    otag = "SBNchi::ReloadCoreSpectrum\t|| ";
    is_verbose = true;

    bool is_fractional = true;
    cholosky_performed = false;

    if(is_verbose) std::cout<<otag<<"Begininning to reload core spec! First Set new core spec"<<std::endl;
    core_spectrum = *bkgin;
    core_spectrum.CollapseVector();

    if(is_verbose) std::cout<<otag<<" || Clear all previous chi^2 data"<<std::endl;
    vec_last_calculated_chi.clear();
    vec_last_calculated_chi.resize(num_bins_total_compressed, std::vector<double>( num_bins_total_compressed,0) );

    //Reset matrix_systematics to fractional
    if(is_verbose) std::cout<<otag<<" Reseting matrix_systematics to matrix_fractional_covariance"<<std::endl;
    matrix_systematics = matrix_fractional_covariance;

    if(matrix_systematics.GetNcols()!=num_bins_total){
        std::cout<<otag<<"ERROR: trying to pass a matrix to SBNchi that isnt the right size"<<std::endl;
        std::cout<<otag<<"ERROR: num_bins_total: "<<num_bins_total<<" and matrix is: "<<matrix_systematics.GetNcols()<<std::endl;
        exit(EXIT_FAILURE);
    }

    if(is_verbose)std::cout<<otag<<"Go from fracCovariance to fullCovariance. matrix_systematics.GetNcols(): "<<matrix_systematics.GetNcols()<<" matrix_systematics.GetNrows(): "<<matrix_systematics.GetNrows()<<" core->fullvec.size(): "<<core_spectrum.full_vector.size()<<std::endl;
    // systematics per scaled event
    for(int i =0; i<matrix_systematics.GetNcols(); i++)
    {
        //std::cout<<"KRAK: "<<core_spectrum.full_vector.at(i)<<std::endl;
        for(int j =0; j<matrix_systematics.GetNrows(); j++)
        {
            if(is_fractional){
                if(std::isnan(matrix_systematics(i,j)) ){
                    matrix_systematics(i,j) = 0;
                }else{
                    matrix_systematics(i,j) = matrix_systematics(i,j)*core_spectrum.full_vector.at(i)*core_spectrum.full_vector.at(j);
                }
                //if(i==j) matrix_systematics(i,j) += pow(core_spectrum.full_err_vector[i],2); //double counting see below.
            }
        }
    }

    //add MC intrinsic error
    if(is_verbose) std::cout<< otag << "Add MC intrinsic error to the full covariance matrix" << std::endl;
    core_spectrum.CalcErrorVector();
    for(int i=0; i<num_bins_total; i++){
         matrix_systematics(i,i) += pow(core_spectrum.full_err_vector.at(i), 2.0); 
    }

    if(is_verbose)std::cout<<otag<<"Filling stats into cov matrix"<<std::endl;
    // Fill stats from the back ground vector
    TMatrixT <double> Mstat(num_bins_total, num_bins_total);
    FillStatsMatrix(Mstat, core_spectrum.full_vector);
    //FillStatsMatrix(Mstat, core_spectrum.full_err_vector);

    if(Mstat.IsSymmetric()){
        if(is_verbose)std::cout<<otag<<"Stat matrix is symmetric (it is just diagonal core)"<<std::endl;
    }else{
        std::cout<<otag<<"ERROR: SBNchi::FormCovarianceMatrix, stats  is not symmetric!"<<std::endl;
        exit(EXIT_FAILURE);
    }

    /*
       TMatrixD Mcorr = matrix_systematics;
    //TEMP
    TFile *f = new TFile("gloop.root","recreate");
    f->cd();
    for(int i=0; i<Mcorr.GetNrows();i++){
    for(int j=0; j<Mcorr.GetNrows();j++){
    Mcorr(i,j) = matrix_systematics(i,j)/(sqrt(matrix_systematics(i,i))*sqrt(matrix_systematics(j,j)));
    }
    }
    TH2D Hcorr(Mcorr);
    TH2D Hsys(matrix_systematics);

    f->cd();
    Hcorr.Write("corr");
    Hsys.Write("total");
    f->Close();
    */


    //And then define the total covariance matrix in all its glory
    TMatrixT <double> Mtotal(num_bins_total, num_bins_total);
    Mtotal.Zero();

    if(is_stat_only){
        if(is_verbose)std::cout<<otag<<"Using stats only in covariance matrix"<<std::endl;
        Mtotal = Mstat;
    }else{
        if(is_verbose)std::cout<<otag<<" Using stats+sys in covariance matrix"<<std::endl;
        Mtotal = Mstat + matrix_systematics;
    }

    //Also going to do matrix_systematics_collapsed;

    m_matrix_systematics_collapsed.ResizeTo(num_bins_total_compressed, num_bins_total_compressed);
    CollapseModes(matrix_systematics, m_matrix_systematics_collapsed);

    if(is_verbose) std::cout<<otag<<"Mstat: "<<Mstat.GetNrows()<<" x "<<Mstat.GetNcols()<<std::endl;
    if(is_verbose) std::cout<<otag<<"matrix_systematics: "<<matrix_systematics.GetNrows()<<" x "<<matrix_systematics.GetNcols()<<std::endl;
    if(is_verbose) std::cout<<otag<<"Mtotal: "<<Mtotal.GetNrows()<<" x "<<Mtotal.GetNcols()<<std::endl;

    if(Mtotal.IsSymmetric() ){
        if(is_verbose)	std::cout<<otag<<"Total Mstat +matrix_systematics is symmetric"<<std::endl;
    }else{

        double tol = m_tolerance;
        double biggest_deviation = 0;
        int bi=0;
        int bj=0;

        if(is_verbose)std::cout<<otag<<"WARNING: Stats + sys result appears to be not symmetric!"<<std::endl;
        for(int i=0; i<Mtotal.GetNrows(); i++){
            for(int j=0; j<Mtotal.GetNcols(); j++){
                double dev = fabs(Mtotal(i,j)-Mtotal(j,i));
                if(dev>biggest_deviation){
                    biggest_deviation = 2*dev/(fabs(Mtotal(i,j))+fabs(Mtotal(j,i)));
                    bi=i;
                    bj=j;
                }
                if(Mtotal(i,j)!=Mtotal(i,j)){

                    std::cout<<"ERROR: we have NAN's  Better check your inputs."<<std::endl;
                    exit(EXIT_FAILURE);

                }
            }
        }

        if(is_verbose) std::cout<<otag<<"WARNING: Biggest Relative Deviation from symmetry is i:"<<bi<<" j: "<<bj<<" of order "<<biggest_deviation<<" M(j,i)"<<Mtotal(bj,bi)<<" M(i,j)"<<Mtotal(bi,bj)<<std::endl;

        if(biggest_deviation >tol){

            std::cout<<otag<<"ERROR: Thats too unsymettric, killing process. Better check your inputs."<<std::endl;

            exit(EXIT_FAILURE);
        }else{

            if(is_verbose)	std::cout<<otag<<"WARNING: Thats within tolderence. Continuing."<<std::endl;
        }
    }

    TMatrixT<double > Mctotal(num_bins_total_compressed,num_bins_total_compressed);

    CollapseModes(Mtotal, Mctotal);

    matrix_collapsed = Mctotal;

    vec_matrix_collapsed = TMatrixDToVector(Mctotal);
    double invdet=0;

    TMatrixD McI(num_bins_total_compressed,num_bins_total_compressed);
    McI.Zero();

    if(is_verbose) std::cout<<otag<<" About to do a SVD decomposition"<<std::endl;
    TDecompSVD svd(Mctotal);
    if (!svd.Decompose()) {

        std::cout <<otag<<"Decomposition failed, matrix not symettric?, has nans?" << std::endl;
        std::cout<<otag<<"ERROR: The matrix to invert failed a SVD decomp!"<<std::endl;

        for(int i=0; i< num_bins_total_compressed; i++){
            for(int j=0; j< num_bins_total_compressed; j++){
                std::cout<<i<<" "<<j<<" "<<Mctotal(i,j)<<std::endl;
            }
        }

        exit(EXIT_FAILURE);
        return 0;
    } else {
        McI = svd.Invert();
    }
    if( !McI.IsValid()){
        std::cout<<otag<<"ERROR: The inverted matrix isnt valid! Something went wrong.."<<std::endl;
        exit(EXIT_FAILURE);

    }

    if(is_verbose)std::cout<<otag<<"SUCCESS! Inverted."<<std::endl;
    // matrix_systematics.Print();
    // McI.Print();
    //Mtotal.Print();

    vec_matrix_inverted = TMatrixDToVector(McI);

    // test for validity
    bool is_small_negative_eigenvalue = false;
    double tolerence_positivesemi = m_tolerance;


    if(Mctotal.IsSymmetric() ){
        std::cout<<otag<<"Total Mctotal +matrix_systematics is symmetric"<<std::endl;
    }else{
        std::cout<<otag<<"Total Mctotal +matrix_systematics isNOT symmetric"<<std::endl;
        Mctotal.Print();
        for(int i=0; i< Mctotal.GetNrows(); i++){
            for(int j=0; j< Mctotal.GetNcols(); j++){
                if(i==j){
                    continue;
                }
                if(i>j) continue;
                Mctotal(i,j)=Mctotal(j,i);
            }
        }
    }

    if(Mctotal.IsSymmetric() ){
        std::cout<<otag<<"Total Mctotal +matrix_systematics is NOW symmetric"<<std::endl;
    }else{
        std::cout<<otag<<"Total Mctotal +matrix_systematics is Still not. Means something is broken?"<<std::endl;
    }
    //if a matrix is (a) real and (b) symmetric (checked above) then to prove positive semi-definite, we just need to check eigenvalues and >=0;
    TMatrixDEigen eigen (Mctotal);
    TVectorD eigen_values = eigen.GetEigenValuesRe();
    TVectorD eigen_valuesIm = eigen.GetEigenValuesIm();

    for(int i=0; i< eigen_values.GetNoElements(); i++){
        if(eigen_values(i)<0){
            is_small_negative_eigenvalue = true;
            if(fabs(eigen_values(i))> tolerence_positivesemi ){
                std::cout<<otag<<" collapsed covariance matrix contains (at least one)  negative eigenvalue: "<<eigen_values(i)<<std::endl;
                std::cout<<otag<<" full covariance "<<std::endl;
                int cc =core_spectrum.GetHistNumber(i); 
                std::cout<<otag<<" collapsed covariance matrix contains (at least one)  negative eigenvalue: "<<eigen_values(i)<<std::endl;
                std::cout<<otag<<" This occurs at element "<<i<<" in hist "<<core_spectrum.GetHistNumber(i)<<" "<<core_spectrum.fullnames[core_spectrum.GetHistNumber(i)]<<std::endl;
                std::cout<<otag<<" Collllvec: "<<core_spectrum.collapsed_vector[cc]<<std::endl;
                std::cout<<otag<<" Negative: "<<eigen_valuesIm(i)<<std::endl;

                for(int r=0; r <core_spectrum.collapsed_vector.size(); r++) std::cout<<core_spectrum.collapsed_vector[r]<<" ";
                std::cout<<std::endl;
                for(int r=0; r <core_spectrum.collapsed_vector.size(); r++) std::cout<<eigen_values(r)<<" ";
                std::cout<<std::endl;
                for(int r=0; r <core_spectrum.collapsed_vector.size(); r++) std::cout<<eigen_valuesIm(r)<<" ";
                std::cout<<std::endl;

                exit(EXIT_FAILURE);
            }
        }
    }


    if(is_small_negative_eigenvalue){
        if(is_verbose)	std::cout<<otag<<"Generated covariance matrix is (allmost) positive semi-definite. It did contain small negative values of absolute value <= :"<<tolerence_positivesemi<<std::endl;
    }else{
        if(is_verbose)	std::cout<<otag<<"Generated covariance matrix is also positive semi-definite."<<std::endl;
    }

    core_spectrum.CollapseVector();

    return 0;
}



/*********************************************
 *		Different Chi^2 calculations
 * ******************************************/

//Standard chi^2 calculation
double SBNchi::CalcChi(SBNspec *sigSpec){
    double tchi = 0;

    if(sigSpec->collapsed_vector.size()==0){
        //        if(is_verbose)	std::cout<<"WARNING: SBNchi::CalcChi, inputted sigSpec has un-compressed vector, I am doing it now, but this is inefficient!"<<std::endl;
        sigSpec->CollapseVector();
    }

    int k=0;

    for(int i =0; i<num_bins_total_compressed; i++){
        for(int j =0; j<num_bins_total_compressed; j++){
            k++;

            if(i==j && fabs(vec_matrix_inverted.at(i).at(j)) > 1e16 && vec_matrix_inverted[i][j] < 0){
                std::cout<<"ERROR: SBNchi::CalcChi || diagonal of inverse covariance is negative! : "<<vec_matrix_inverted[i][j]<<" @ ("<<i<<","<<j<<")"<<std::endl;
            }
            std::cout<<"HOME: "<<i<<" "<<j<<" "<<core_spectrum.collapsed_vector.at(i)<<" "<<core_spectrum.collapsed_vector.at(j)<<" "<<vec_matrix_inverted.at(i).at(j)<<" "<<-sigSpec->collapsed_vector.at(i)<<" "<<-sigSpec->collapsed_vector.at(j)<<std::endl;
            vec_last_calculated_chi.at(i).at(j) =(core_spectrum.collapsed_vector.at(i)-sigSpec->collapsed_vector.at(i))*vec_matrix_inverted.at(i).at(j)*(core_spectrum.collapsed_vector.at(j)-sigSpec->collapsed_vector.at(j) );
            tchi += vec_last_calculated_chi.at(i).at(j);
        }
    }

    last_calculated_chi = tchi;
    return tchi;
}


float SBNchi::CalcChi(std::vector<float> * sigVec){
    float tchi = 0;

#ifndef _OPENACC
    if(sigVec->size() != num_bins_total_compressed ){
        std::cerr<<"ERROR: SBNchi::CalcChi(std::vector<float>) ~ your inputed vector does not have correct dimensions"<<std::endl;
        std::cerr<<"sigVec.size(): "<<sigVec->size()<<" num_bins_total_compressed: "<<num_bins_total_compressed<<std::endl;
        exit(EXIT_FAILURE);
    }
#endif

    for(int i =0; i<num_bins_total_compressed; i++){
        for(int j =0; j<num_bins_total_compressed; j++){
            tchi += (core_spectrum.collapsed_vector[i]-(*sigVec)[i])*vec_matrix_inverted[i][j]*(core_spectrum.collapsed_vector[j]-(*sigVec)[j] );
        }
    }

    //last_calculated_chi = tchi;
    return tchi;
}


double SBNchi::CalcChi(std::vector<double> * sigVec){

    double tchi = 0;

#ifndef _OPENACC
    if(sigVec->size() != num_bins_total_compressed ){
        std::cerr<<"ERROR: SBNchi::CalcChi(std::vector<double>) ~ your inputed vector does not have correct dimensions"<<std::endl;
        std::cerr<<"sigVec.size(): "<<sigVec->size()<<" num_bins_total_compressed: "<<num_bins_total_compressed<<std::endl;
        exit(EXIT_FAILURE);
    }
#endif

    for(int i =0; i<num_bins_total_compressed; i++){
        for(int j =0; j<num_bins_total_compressed; j++){
            tchi += (core_spectrum.collapsed_vector[i]-(*sigVec)[i])*vec_matrix_inverted[i][j]*(core_spectrum.collapsed_vector[j]-(*sigVec)[j] );
        }
    }

    //last_calculated_chi = tchi;
    return tchi;
}

double SBNchi::CalcChi(double* sigVec){
    double tchi = 0;

    for(int i =0; i<num_bins_total_compressed; i++){
        for(int j =0; j<num_bins_total_compressed; j++){
            tchi += (core_spectrum.collapsed_vector[i]-(sigVec)[i])*vec_matrix_inverted[i][j]*(core_spectrum.collapsed_vector[j]-(sigVec)[j] );
        }
    }

    //last_calculated_chi = tchi;
    return tchi;
}

float SBNchi::CalcChi(float **invert_matrix, float* core, float *sig){
    float tchi = 0;

    for(int i =0; i<num_bins_total_compressed; i++){
        for(int j =0; j<num_bins_total_compressed; j++){

            tchi += (core[i]-sig[i])*invert_matrix[i][j]*(core[j]-sig[j] );
        }
    }

    return tchi;
}

float SBNchi::PoissonLogLiklihood(float * pred, float *data){
    float ans = 0;

    int n = num_bins_total_compressed;
    for(int i =0; i< n; i++){
        ans +=   ( data[i] >0.001 ? 2.0*(pred[i]-data[i]  + data[i]*log(data[i]/pred[i]) ) : 2*pred[i]);
    }
    return ans;
}


double SBNchi::CalcChi(double **invert_matrix, double* core, double *sig){
    double tchi = 0;

    for(int i =0; i<num_bins_total_compressed; i++){
        for(int j =0; j<num_bins_total_compressed; j++){
            tchi += (core[i]-sig[i])*invert_matrix[i][j]*(core[j]-sig[j]);
        }
    }

    return tchi;
}

double SBNchi::CalcChi(TMatrixT<double> M_invert, std::vector<double>& spec, std::vector<double>& data){

	return this->CalcChi(M_invert, spec, data, false);
}

double SBNchi::CalcChi(TMatrixT<double> M_invert, std::vector<double>& spec, std::vector<double>& data, bool print){
      double tchi = 0;
      for(int i=0; i< num_bins_total_compressed; i++){
          for(int j=0; j< num_bins_total_compressed ;j++){
              double this_contribution = M_invert(i,j)*(spec[i]- data[i])*(spec[j]-data[j]);
              tchi += this_contribution;
              if(print) std::cout<<i<<" "<<j<<" "<<this_contribution<<", "<<tchi<<", "<<" Invar: "<<M_invert(i,j)<<", "<<spec[i]<<", "<<data[i]<<", "<<spec[j]<<", "<<data[j]<<std::endl;
          }
      }
      return tchi;
}

//same as above but passing in a vector instead of whole SBNspec
double SBNchi::CalcChi(std::vector<double> sigVec){
    double tchi = 0;

    if(sigVec.size() != num_bins_total_compressed ){
        std::cerr<<"ERROR: SBNchi::CalcChi(std::vector<double>) ~ your inputed vector does not have correct dimensions"<<std::endl;
        std::cerr<<"sigVec.size(): "<<sigVec.size()<<" num_bins_total_compressed: "<<num_bins_total_compressed<<std::endl;
        exit(EXIT_FAILURE);
    }

    for(int i =0; i<num_bins_total_compressed; i++){
        for(int j =0; j<num_bins_total_compressed; j++){
            tchi += (core_spectrum.collapsed_vector[i]-sigVec[i])*vec_matrix_inverted[i][j]*(core_spectrum.collapsed_vector[j]-sigVec[j] );
        }
    }

    last_calculated_chi = tchi;
    return tchi;
}


//a shape-fit chi 
double SBNchi::CalcShapeChi(TMatrixT<double> inverted_matrix, std::vector<double>& data, std::vector<double>& constrained_vec, std::vector<double>& free_vec, bool is_background){
	double tchi = 0.0;
	int matrix_bins = inverted_matrix.GetNcols();    
	// check dimension
	if(data.size() != matrix_bins || constrained_vec.size()!= matrix_bins || free_vec.size()!= matrix_bins){
        	std::cerr<<"ERROR: SBNchi::CalcShapeChi ~ your inputed vector does not have correct dimensions"<<std::endl;
        	std::cerr<<"data.size(): "<<data.size()<<" num_bins of matrix "<< matrix_bins<<std::endl;
        	std::cerr<<"constrained_vec.size(): "<< constrained_vec.size()<<" num_bins of matrix "<< matrix_bins<<std::endl;
        	std::cerr<<"free_vec.size(): "<<free_vec.size()<<" num_bins of matrix: "<< matrix_bins<<std::endl;
        	exit(EXIT_FAILURE);
   	}



	std::vector<double> free_true;
	free_true.clear();
	// if 'free_vec' is the background hypothesis, then you need to normalize it to the total excess
	// if 'free_vec' is signal hypothesis, then no need to normalize it to total excess
	if(is_background){
		double sum_data = std::accumulate(data.begin(), data.end(), 0.0);
		double sum_constrain = std::accumulate(constrained_vec.begin(), constrained_vec.end(), 0.0);
		double sum_free = std::accumulate(free_vec.begin(), free_vec.end(), 0.0);
		for(int i =0; i< free_vec.size(); i++)
			free_true.push_back(free_vec[i]*(sum_data-sum_constrain)/sum_free);
	}	
	else free_true = free_vec;

	for(int i=0 ; i< matrix_bins; i++){
		for(int j=0; j< matrix_bins; j++){
			tchi += inverted_matrix(i,j)*(data[i]-constrained_vec[i]-free_true[i])*(data[j]-constrained_vec[j]-free_true[j]);
		}
	}

	return tchi;
}




//A log-lilihood based one used @ MiniBooNE
double SBNchi::CalcChiLog(SBNspec *sigSpec){
    double tchi = 0;

    if(sigSpec->collapsed_vector.size()==0){
        std::cout<<"WARNING: SBNchi::CalcChi, inputted sigSpec has un-compressed vector, I am doing it now, but this is inefficient!"<<std::endl;
        sigSpec->CollapseVector();
    }

    for(int i =0; i<num_bins_total_compressed; i++){
        for(int j =0; j<num_bins_total_compressed; j++){
            vec_last_calculated_chi.at(i).at(j) =(core_spectrum.collapsed_vector[i]-sigSpec->collapsed_vector[i])*vec_matrix_inverted[i][j]*(core_spectrum.collapsed_vector[j]-sigSpec->collapsed_vector[j] );
            tchi += vec_last_calculated_chi.at(i).at(j);
        }
    }

    double absDetM = log(fabs(matrix_collapsed.Determinant()));

    last_calculated_chi = tchi+absDetM;
    return tchi+absDetM;
}


double SBNchi::CalcChi(SBNspec *sigSpec, SBNspec *obsSpec){
    double tchi=0;
    if(sigSpec->collapsed_vector.size()==0){
        std::cout<<"WARNING: SBNchi::CalcChi, inputted sigSpec has un-compressed vector, I am doing it now, but this is inefficient!"<<std::endl;
        sigSpec->CollapseVector();
    }

    if(obsSpec->collapsed_vector.size()==0){
        std::cout<<"WARNING: SBNchi::CalcChi, inputted obsSpec has un-compressed vector, I am doing it now, but this is inefficient!"<<std::endl;
        obsSpec->CollapseVector();
    }

    for(int i =0; i<num_bins_total_compressed; i++){
        for(int j =0; j<num_bins_total_compressed; j++){
            tchi += (obsSpec->collapsed_vector[i]-sigSpec->collapsed_vector[i])*vec_matrix_inverted[i][j]*(obsSpec->collapsed_vector[j]-sigSpec->collapsed_vector[j] );
        }
    }

    last_calculated_chi = tchi;
    return tchi;
}



/**************************************************************************
 *			Collapsing code
 * ************************************************************************/


//This is the powerhouse, takes each detector matrix filled with num_channels channels of num_subchannels[i] subchannels, and collapses it.
void SBNchi::CollapseSubchannels(TMatrixT <double> & M, TMatrixT <double> & Mc){
    bool debug = false;
    if(debug)	std::cout<<"Starting:M "<<M.GetNcols()<<" "<<M.GetNrows()<<" "<<115<<std::endl;
    if(debug)	std::cout<<"Starting:Mc "<<Mc.GetNcols()<<" "<<Mc.GetNrows()<<" "<<30<<std::endl;

    std::vector<std::vector<TMatrixT<double>>> Summed(num_channels, std::vector<TMatrixT<double>>(num_channels) );	//Initialise a matrix of matricies, to ZERO.
    for(int ic = 0; ic < num_channels; ic++){
        for(int jc =0; jc < num_channels; jc++){
            Summed[ic][jc].ResizeTo(num_bins[jc],num_bins[ic]) ;// This is CORRECT, do not switch (ie Summed[0][1] = size (num_bins[1], num_bins[0])
            Summed[ic][jc] = 0.0;
        }
    }

    int mrow = 0.0;
    int mcol = 0.0;

    for(int ic = 0; ic < num_channels; ic++){ 	 //Loop over all rows
        for(int jc =0; jc < num_channels; jc++){ //Loop over all columns

            if(debug)std::cout<<"Diagonal! : "<<ic<<" "<<jc<<" mcol is: "<<mcol<<" mrow is: "<<mrow<<std::endl;

            for(int m=0; m < num_subchannels[ic]; m++){
                for(int n=0; n< num_subchannels[jc]; n++){ //For each big block, loop over all subchannels summing toGether
                    Summed[ic][jc] +=  M.GetSub(mrow+n*num_bins[jc] ,mrow + n*num_bins[jc]+num_bins[jc]-1, mcol + m*num_bins[ic], mcol+ m*num_bins[ic]+num_bins[ic]-1 );
                }
            }


            mrow += num_subchannels[jc]*num_bins[jc];//As we work our way left in columns, add on that many bins
        }//end of column loop

        mrow = 0; // as we end this row, reSet row count, but jump down 1 column
        mcol += num_subchannels[ic]*num_bins[ic];
    }//end of row loop

    ///********************************* And put them back toGether! ************************//
    Mc.Zero();
    mrow = 0;
    mcol = 0;

    //Repeat again for Contracted matrix
    for(int ic = 0; ic < num_channels; ic++){
        for(int jc =0; jc < num_channels; jc++){

            Mc.SetSub(mrow,mcol,Summed[ic][jc]);
            mrow += num_bins[jc];
        }

        mrow = 0;
        mcol +=num_bins[ic];
    }

    return;
}



//This is the detector layer, Take a given mode and run over each detector V detector sub matrix
void SBNchi::CollapseDetectors(TMatrixT <double> & M, TMatrixT <double> & Mc){

    Mc.Zero();
    int nrow = num_bins_detector_block;// N_e_bins*N_e_spectra+N_m_bins*N_m_spectra;
    int crow = num_bins_detector_block_compressed; //N_e_bins+N_m_bins;

    for(int m =0; m< num_detectors; m++){
        for(int n =0; n< num_detectors; n++){
            TMatrixT<double> imat(nrow,nrow);
            TMatrixT<double> imatc(crow,crow);

            imat = M.GetSub(n*nrow,n*nrow+nrow-1, m*nrow,m*nrow+nrow-1);
            CollapseSubchannels(imat,imatc);
            Mc.SetSub(n*crow,m*crow,imatc);
        }
    }

    return;
}

//This is the Mode layer, Take a given full matrix and runs over each Mode V Mode sub matrix
void SBNchi::CollapseModes(TMatrixT <double> & M, TMatrixT <double> & Mc){

    Mc.Zero();
    int nrow = num_bins_mode_block;// (N_e_bins*N_e_spectra+N_m_bins*N_m_spectra)*N_dets;
    int crow=  num_bins_mode_block_compressed;// (N_e_bins+N_m_bins)*N_dets;

    for(int m =0; m< num_modes ; m++){
        for(int n =0; n< num_modes; n++){

            TMatrixT<double> imat(nrow,nrow);
            TMatrixT<double> imatc(crow,crow);

            imat = M.GetSub(n*nrow,n*nrow+nrow-1, m*nrow,m*nrow+nrow-1);

            CollapseDetectors(imat,imatc);
            Mc.SetSub(n*crow,m*crow,imatc);

        }
    }

    return;
}


TMatrixT<double> SBNchi::CalcNeymanCovarianceMatrix(TMatrixT<double>* frac_covar, std::vector<double>& mc_full, std::vector<double>& mc_full_err,  std::vector<double>& data_full){
	TMatrixT<double> Mfull(frac_covar->GetNcols(), frac_covar->GetNcols());
	TMatrixT<double> Mout(num_bins_total_compressed, num_bins_total_compressed);	

	for(int i =0 ; i<frac_covar->GetNcols(); i++){
		for(int j =0; j< frac_covar->GetNcols(); j++){
			if(std::isnan( (*frac_covar)(i,j) )){
				Mfull(i,j) = 0.0;
			}
			else Mfull(i,j) = (*frac_covar)(i,j) *mc_full[i]*mc_full[j];
			if(i == j){
			      Mfull(i,j) += data_full[i];
			      if(!is_stat_only) Mfull(i,j) += pow(mc_full_err[i], 2.0);
			}
		}
	}

	CollapseModes(Mfull, Mout);
	return Mout;
}

TMatrixT<double> SBNchi::CalcCovarianceMatrix(TMatrixT<double>*M, std::vector<double>& spec, bool add_stats){

    TMatrixT<double> Mout(M->GetNcols(), M->GetNcols() );
    for(int i =0; i<M->GetNcols(); i++)
    {
        for(int j =0; j<M->GetNrows(); j++)
        {
            if(  std::isnan( (*M)(i,j) )){
                Mout(i,j) = 0.0;
            }else{

                Mout(i,j) = (*M)(i,j)*spec[i]*spec[j];
            }
            if(add_stats){  if(i==j) Mout(i,i) += spec[i]; }  //stats part
        }
    }
    return Mout;
}


// add stat matrix, could be used for Neyman or Pearson statistic
// can be used either for collapsed matrix or uncollapsed matrix
TMatrixT<double> SBNchi::AddStatMatrix(TMatrixT<double>*M,  const std::vector<double>& datavec ){

    if(M->GetNcols() != datavec.size()){
	std::cout << "ERROR: AddStatMatrix: size of covariance matrix doesn't match size of the spectrum" << std::endl;
	exit(EXIT_FAILURE);
    }

    TMatrixT<double> Mout(M->GetNcols(), M->GetNcols() );
    for(int i =0; i<M->GetNcols(); i++)
    {
        for(int j =0; j<M->GetNrows(); j++)
        {
	    Mout(i,j) = (*M)(i,j);
            if(i==j) Mout(i,i) += datavec[i];
        }
    }
    return Mout;
}





    
    //generate Pearson covariance matrix(uncollapsed matrix)
TMatrixT<double> SBNchi::CalcCovarianceMatrix(TMatrixT<double>*M, TVectorT<double>& spec, TVectorT<double>& spec_err){

    TMatrixT<double> Mout( M->GetNcols(), M->GetNcols() );
    // systematics per scaled event
    for(int i =0; i<M->GetNcols(); i++)
    {
        //std::cout<<"KRAK: "<<core_spectrum.full_vector.at(i)<<std::endl;
        for(int j =0; j<M->GetNrows(); j++)
        {
            if(  std::isnan( (*M)(i,j))){
                Mout(i,j) = 0.0;
            }else{
                Mout(i,j) = (*M)(i,j)*spec(i)*spec(j);
            }
            if(i==j){
		 Mout(i,i) +=spec(i);
		 //add MC intrinsic error if it's not stat-only case
		 //if(!is_stat_only) Mout(i,i) += pow(spec_err(i), 2.0);
	    }
        }
    }
    return Mout;
}





TMatrixT<double> SBNchi::CalcCovarianceMatrix(TMatrixT<double>*M, std::vector<double>& spec, std::vector<double>& spec_err){

    TMatrixT<double> Mout(M->GetNcols(), M->GetNcols() );

    for(int i =0; i<M->GetNcols(); i++)
    {
        for(int j =0; j<M->GetNrows(); j++)
        {
            if(  std::isnan( (*M)(i,j) )){
                Mout(i,j) = 0.0;
            }else{

                Mout(i,j) = (*M)(i,j)*spec[i]*spec[j];
            }
            if(i==j){
		 Mout(i,i) += spec[i];   //stats part
		 if(!is_stat_only) Mout(i,i) += pow(spec_err[i], 2.0); //MC intrinsic error part 
	    }
        }
    }
    return Mout;
}



//split covariance matrix into shape-only, normalizatoin-only, mixed covariance matrix
TMatrixT<double> SBNchi::SplitCovarianceMatrix(TMatrixT<double>* frac_covar, std::vector<double>& spec, std::vector<double>& spec_err, int which_process){
	
	int matrix_bins = frac_covar->GetNcols();
	if(matrix_bins != spec.size()){
		std::cout << "SBNchi::SplitCovarianceMatrix||\t Dimension of vector "<< spec.size() << " doesn't match dimension of covariance matrix " << matrix_bins << std::endl;
		exit(EXIT_FAILURE);
	}

	TMatrixT<double> full_covar(matrix_bins, matrix_bins);
	TMatrixT<double> Mout(matrix_bins, matrix_bins);

	for(int i=0; i<matrix_bins; i++){
		for(int j=0; j< matrix_bins; j++){
			if( std::isnan( (*frac_covar)(i,j)  )) full_covar(i,j) =0;
			else full_covar(i,j) = (*frac_covar)(i,j)*spec[i]*spec[j]; 
			if((i == j) && (!is_stat_only)) full_covar(i,j) += pow(spec_err[i], 2.0);
		}
	}	

	double N_T = std::accumulate(spec.begin(), spec.end(), 0.0);
	if(N_T == 0){
		Mout.Zero();
		return Mout;
	}
	double f_1 = full_covar.Sum()/pow(N_T, 2.0);
	std::vector<double> P_sum;
	for(int i=0; i< matrix_bins; i++){
                double P_temp = 0;
                for(int j=0; j< matrix_bins; j++) P_temp += full_covar(i,j);
                P_sum.push_back(P_temp);
        }

	switch(which_process)
	{
		case 1:
		//get shape-only covariance matrix
		    std::cout << "SBNchi::SplitCovarianceMatrix||\tGet shape-only covariance matrix" << std::endl;
		    for(int i=0; i<matrix_bins; i++){
			for(int j=0 ;j<matrix_bins; j++)
		 	    Mout(i,j) = full_covar(i,j) - (spec[j]*P_sum[i]+spec[i]*P_sum[j])/N_T + spec[i]*spec[j]*f_1;
		    }
		    break;
		case 2:
		//get mixed covariance matrix
		     std::cout << "SBNchi::SplitCovarianceMatrix||\tGet mixed covariance matrix" << std::endl;
		     for(int i=0; i<matrix_bins; i++){
                        for(int j=0 ;j<matrix_bins; j++)
			    Mout(i,j) = (spec[j]*P_sum[i]+spec[i]*P_sum[j])/N_T - 2.0*spec[i]*spec[j]*f_1;
		     } 
		     break;
		case 3:
		//get normalization only covariance matrix
		    std::cout << "SBNchi::SplitCovarianceMatrix||\tGet normalization-only covariance matrix" << std::endl;
		    for(int i=0; i<matrix_bins; i++){
                        for(int j=0 ;j<matrix_bins; j++){
			    //Mout(i,j) = spec[i]*spec[j]*f_1;
			    Mout(i,j) = f_1;  //fractional normalization only matrix
			}
		    }
		    break;			
		default:
		//return full_covariance matrix 
		     break;
	}
	
	return Mout;
}


//given fractional covariance matrix, MC predicted spctrum, bkgd spectrum
//return  shape only systematic covariance matrix
TMatrixT<double> SBNchi::CalcShapeOnlyCovarianceMatrix(TMatrixT<double> &M, SBNspec *mc, SBNspec* bkg){


	int mc_num_bins = mc->num_bins_total;
	std::vector<double> mc_full = mc->full_vector;
	mc->CalcErrorVector();
	std::vector<double> mc_full_err = mc->full_err_vector;
	std::vector<double> bkgd_full = bkg->full_vector;
		
	TMatrixT<double> full_systematic(mc_num_bins, mc_num_bins);
	TMatrixT<double> full_shape_covar(mc_num_bins,mc_num_bins);

	if(mc_full.size() != M.GetNcols()){
		std::cout << "Dimension of MC  full vector " << mc_num_bins<< " does not match dimension of covariance matrix:" << M.GetNcols() << std::endl;
		exit(EXIT_FAILURE);
	}	

        //fill the usual systematic covariance matrix
	for(int i=0; i< mc_num_bins; i++){
		for(int j=0; j< mc_num_bins; j++){
			if( std::isnan(  M(i,j)  )) full_systematic(i,j) = 0.0;
			else full_systematic(i,j) = M(i,j)*mc_full[i]*mc_full[j];

			if((i==j) && (!is_stat_only)) full_systematic(i,j) += pow(mc_full_err[i], 2.0);
		}
	}

	double sum_bkd = std::accumulate(bkgd_full.begin(), bkgd_full.end(), 0.0);
	if(sum_bkd == 0){
		full_shape_covar.Zero();
		return full_shape_covar;
	}
	
        double N = full_systematic.Sum()/pow(sum_bkd, 2.0);

        //vector of sum over rows for collapsed syst covariance matrix
	std::vector<double> P_sum;
	for(int i=0; i< mc_num_bins; i++){
		double P_temp = 0;
		for(int j=0; j< mc_num_bins; j++) P_temp += full_systematic(i,j);
		P_sum.push_back(P_temp);
	}
	
	//construct shape only systematic covariance matrix
	for(int i=0; i< mc_num_bins; i++){
		for(int j=0 ;j< mc_num_bins; j++){
			//shape only
            		full_shape_covar(i,j) = full_systematic(i,j)+bkgd_full[i]*bkgd_full[j]*N - (bkgd_full[i]*P_sum[j]+bkgd_full[j]*P_sum[i])/sum_bkd;
		}
	}	

	return full_shape_covar;
}


//calculate the shape+mixed part of a covariance matrix
TMatrixT<double> SBNchi::CalcShapeMixedCovarianceMatrix(TMatrixT<double> &M, SBNspec *mc, SBNspec* bkg){


	int mc_num_bins = mc->num_bins_total;
	std::vector<double> mc_full = mc->full_vector;
	mc->CalcErrorVector();
        std::vector<double> mc_full_err = mc->full_err_vector;
	std::vector<double> bkgd_full = bkg->full_vector;
		
	TMatrixT<double> full_systematic(mc_num_bins, mc_num_bins);
	TMatrixT<double> full_shape_covar(mc_num_bins,mc_num_bins);

	if(mc_full.size() != M.GetNcols()){
		std::cout << "Dimension of MC  full vector " << mc_num_bins<< " does not match dimension of covariance matrix:" << M.GetNcols() << std::endl;
		exit(EXIT_FAILURE);
	}	


        //fill the usual systematic covariance matrix
	for(int i=0; i< mc_num_bins; i++){
		for(int j=0; j< mc_num_bins; j++){
			if( std::isnan(  M(i,j)  )) full_systematic(i,j) = 0.0;
			else full_systematic(i,j) = M(i,j)*mc_full[i]*mc_full[j];

			if((i==j) && (!is_stat_only)) full_systematic(i,j)+= pow(mc_full_err[i], 2.0);
		}
	}



	double sum_bkd = std::accumulate(bkgd_full.begin(), bkgd_full.end(), 0.0);
	if(sum_bkd == 0){
		full_shape_covar.Zero();
		return full_shape_covar;
	}
	double N = full_systematic.Sum()/pow(sum_bkd, 2.0);
	//vector of sum over rows for collapsed syst covariance matrix
	std::vector<double> P_sum;
	for(int i=0; i< mc_num_bins; i++){
		double P_temp = 0;
		for(int j=0; j< mc_num_bins; j++) P_temp += full_systematic(i,j);
		P_sum.push_back(P_temp);
	}
	
	//construct shape only systematic covariance matrix
	for(int i=0; i< mc_num_bins; i++){
		for(int j=0 ;j< mc_num_bins; j++){
			full_shape_covar(i,j) = full_systematic(i,j) - bkgd_full[i]*bkgd_full[j]*N;
		}
	}	

	return full_shape_covar;
}


double SBNchi::CalcChi_statonlyCNP(std::vector<double> &pred, std::vector<double>& data){
      std::vector<double> diag(pred.size());
      for(int j =0; j<pred.size(); j++)
       {   
          //diag[j] =  ( data[j] >0.001 ? 3.0/(1.0/data[j] +  2.0/pred[j])  : pred[j]/2.0 );
          diag[j] =  data[j];
      }
      double tchi = 0.0;
      for(int i =0; i<pred.size(); i++){
              tchi += pow(pred[i]-data[i],2)/diag[i];
      }   
      return tchi;
}


//here spec is full vector of MC, spec_collapse is collapsed vector of MC, datavec is collapsed vector of data
TMatrixT<double> SBNchi::CalcCovarianceMatrixCNP(TMatrixT<double> &M, std::vector<double>& spec, std::vector<double>& spec_err, std::vector<double>& spec_collapse, const std::vector<double>& datavec ){

    if(M.GetNcols() != spec.size()){
        std::cout << "ERROR: your input vector does not have the right dimenstion  " << std::endl; 
        std::cout << "Fractional Matrix size :"<< M.GetNcols() << " || Input Full Vector size "<< spec.size() << std::endl;  
        exit(EXIT_FAILURE);
    }

    TMatrixT<double> M_temp(M.GetNcols(), M.GetNcols() );
    TMatrixT<double> Mout(spec_collapse.size(), spec_collapse.size()); //collapsed covariance matrix

    //systematic apart 
    for(int i =0; i<M.GetNcols(); i++)
    {
        for(int j =0; j<M.GetNrows(); j++)
        {
            if(  std::isnan( M(i,j) )) {
                M_temp(i,j) = 0.0;
            }else{

                M_temp(i,j) = M(i,j)*spec[i]*spec[j];
	            
                if( (i==j) && (!is_stat_only)) M_temp(i,j)+= pow(spec_err[i], 2.0);
            }

        }
    }

    CollapseModes(M_temp, Mout);
    //add stats part	
    for(int i=0; i< spec_collapse.size(); i++){
        Mout(i,i) +=   ( datavec[i] >0.001 ? 3.0/(1.0/datavec[i] +  2.0/spec_collapse[i])  : spec_collapse[i]/2.0 ); 
        //Mout(i,i) +=   spec_collapse[i];//( datavec[i] >0.001 ? 3.0/(1.0/datavec[i] +  2.0/spec_collapse[i])  : spec_collapse[i]/2.0 ); 
    }
    return Mout;
}


//return collapsed syst+stat matrix (in CNP method)//float.., currently used. 
TMatrixT<double> SBNchi::CalcCovarianceMatrixCNP(TMatrixT<double>* M, std::vector<double>& spec, std::vector<double>& spec_err, std::vector<double>& spec_collapse, const std::vector<float>& datavec ){

    if(M->GetNcols() != spec.size()){
	 std::cout << "ERROR: your input vector does not have the right dimenstion  " << std::endl; 
	 std::cout << "Fractional Matrix size :"<< M->GetNcols() << " || Input Full Vector size "<< spec.size() << std::endl;  
	 exit(EXIT_FAILURE);
    }

    TMatrixT<double> M_temp(M->GetNcols(), M->GetNcols() );
    TMatrixT<double> Mout(spec_collapse.size(), spec_collapse.size()); //collapsed covariance matrix
  
    //systematic apart 
    for(int i =0; i<M->GetNcols(); i++)
    {
        for(int j =0; j<M->GetNrows(); j++)
        {
            if(  std::isnan( (*M)(i,j) )){
                M_temp(i,j) = 0.0;
            }else{

                M_temp(i,j) = (*M)(i,j)*spec[i]*spec[j];
	            if((i==j) && (!is_stat_only)) M_temp(i,j)+= pow(spec_err[i], 2.0);
            }
        }
    }
  
    CollapseModes(M_temp, Mout);
    //add stats part	
    for(int i=0; i< spec_collapse.size(); i++){
	Mout(i,i) +=   ( datavec[i] >0.001 ? 3.0/(1.0/datavec[i] +  2.0/spec_collapse[i])  : spec_collapse[i]/2.0 ); 
	//Mout(i,i) += datavec[i];
   }
    return Mout;
}

// add stat part to the collapsed systematic covariance matrix: CNP chi
TMatrixT<double> SBNchi::AddStatMatrixCNP(TMatrixT<double>*M, std::vector<double>& spec, const std::vector<double>& datavec ){
    TMatrixT<double> Mout(M->GetNcols(), M->GetNcols() );
    for(int i =0; i<M->GetNcols(); i++)
    {
        for(int j =0; j<M->GetNrows(); j++)
        {
	    Mout(i,j) = (*M)(i,j);
            if(i==j) Mout(i,i) +=   ( datavec[i] >0.001 ? 3.0/(1.0/datavec[i] +  2.0/spec[i])  : spec[i]/2.0 );
        }
    }
    return Mout;
}

//this function is wrong
TMatrixT<double> SBNchi::CalcCovarianceMatrixCNP(TMatrixT<double>*M, std::vector<double>& spec, const std::vector<float>& datavec ){

    TMatrixT<double> Mout(M->GetNcols(), M->GetNcols() );

    for(int i =0; i<M->GetNcols(); i++)
    {
        for(int j =0; j<M->GetNrows(); j++)
        {
            if(std::isnan( (*M)(i,j) )){
                Mout(i,j) = 0.0;
            }else{
                Mout(i,j) = (*M)(i,j)*spec[i]*spec[j];
            }
            if(i==j) Mout(i,i) +=   (datavec[i] >0.001 ? 3.0/(1.0/datavec[i] +  2.0/spec[i])  : spec[i]/2.0 );
        }
    }
    return Mout;
}


float SBNchi::CalcChi_Pearson(float * pred, float* data){

    is_verbose = false;
    TMatrixT<double> inverse_collapsed = m_matrix_systematics_collapsed;
    //Add on Pearson Stats terms proportional to data
    for(int j =0; j<num_bins_total_compressed; j++)
    {
        inverse_collapsed(j,j) +=  pred[j];
    }
    inverse_collapsed = this->InvertMatrix(inverse_collapsed);   

    float tchi = 0.0;
    for(int i =0; i<num_bins_total_compressed; i++){
        for(int j =0; j<num_bins_total_compressed; j++){
            tchi += (pred[i]-data[i])*inverse_collapsed(i,j)*(pred[j]-data[j]);
        }
    }
    return tchi;
}


float SBNchi::CalcChi_CNP(float * pred, float* data){

    is_verbose = false;
    TMatrixT<double> inverse_collapsed = m_matrix_systematics_collapsed;
    //Add on CNP terms proportional to data
    for(int j =0; j<num_bins_total_compressed; j++)
    {
        inverse_collapsed(j,j) +=  ( data[j] >0.001 ? 3.0/(1.0/data[j] +  2.0/pred[j])  : pred[j]/2.0 );
    }

    inverse_collapsed = this->InvertMatrix(inverse_collapsed);   

    float tchi = 0.0;
    for(int i =0; i<num_bins_total_compressed; i++){
        for(int j =0; j<num_bins_total_compressed; j++){
            tchi += (pred[i]-data[i])*inverse_collapsed(i,j)*(pred[j]-data[j]);
        }
    }
    return tchi;
}




TMatrixT<double> SBNchi::CalcCovarianceMatrix(TMatrixT<double>*M, TVectorT<double>& spec,bool add_stats){

    TMatrixT<double> Mout( M->GetNcols(), M->GetNcols() );
    // systematics per scaled event
    for(int i =0; i<M->GetNcols(); i++)
    {
        //std::cout<<"KRAK: "<<core_spectrum.full_vector.at(i)<<std::endl;
        for(int j =0; j<M->GetNrows(); j++)
        {
            if(  std::isnan( (*M)(i,j))){
                Mout(i,j) = 0.0;
            }else{
                Mout(i,j) = (*M)(i,j)*spec(i)*spec(j);
            }
            if(add_stats){           if(i==j) Mout(i,i) +=spec(i);}
        }
    }
    return Mout;
}


TMatrixT<double> SBNchi::InvertMatrix(TMatrixT<double> &M){

    double invdet=0;

    TMatrixT<double> McI(M.GetNrows(),M.GetNrows());
    McI.Zero();

    otag = "SBNchi||\tInvertMatrix: ";

    //check the matrix is symmetric
    if(M.IsSymmetric() ){
        if(is_verbose)  std::cout<<otag<<"Covariance matrix is symmetric"<<std::endl;
    }else{

        //double tol = 1e-13;
        double tol = 1e-7; //run with a relaxed tolerance
        double biggest_deviation = 0;
        int bi =0;
        int bj=0;

        if(is_verbose) std::cout<<otag<<"WARNING: this covariance matrix appears to be not symmetric!"<<std::endl;
        for(int i=0; i<M.GetNrows(); i++){
            for(int j=0; j<M.GetNcols(); j++){
                double dev = fabs(M(i,j)-M(j,i));
                if(dev>biggest_deviation){
                    biggest_deviation = 2*dev/(fabs(M(i,j))+fabs(M(j,i)));
                    bi=i;
                    bj=j;
                }
                if(M(i,j)!=M(i,j)){

                    std::cout<<"ERROR: we have NAN's  Better check your inputs."<<std::endl;
                    exit(EXIT_FAILURE);

                }
            }
        }

        if(is_verbose) std::cout<<otag<<"WARNING: Biggest Relative Deviation from symmetry is i:"<<bi<<" j: "<<bj<<" of order "<<biggest_deviation<<" M(j,i)"<<M(bj,bi)<<" M(i,j)"<<M(bi,bj)<<std::endl;

        if(biggest_deviation >tol){

            std::cout<<otag<<"ERROR: Thats too unsymettric, killing process. Better check your inputs."<<std::endl;
            std::cout<<otag<<"ERROR: Biggest Relative Deviation from symmetry is i:"<<bi<<" j: "<<bj<<" of order "<<biggest_deviation<<" M(j,i)"<<M(bj,bi)<<" M(i,j)"<<M(bi,bj)<<std::endl;

            exit(EXIT_FAILURE);
        }else{

            if(is_verbose)      std::cout<<otag<<"WARNING: Thats within tolderence. Continuing."<<std::endl;
        }
    }

    if(is_verbose) std::cout<<otag<<" About to do a SVD decomposition"<<std::endl;
    TDecompSVD svd(M);

    if (!svd.Decompose()){
        std::cout<<otag<<" (InvertMatrix) Decomposition failed, matrix not symettric?, has nans?" << std::endl;
        std::cout<<otag<<"ERROR: The matrix to invert failed a SVD decomp!"<<std::endl;

        for(int i=0; i< M.GetNrows(); i++){
            for(int j=0; j< M.GetNrows(); j++){
                std::cout<<i<<" "<<j<<" "<<M(i,j)<<std::endl;
            }
        }

        exit(EXIT_FAILURE);

    } else {
        McI = svd.Invert();
    }
    if( !McI.IsValid()){
        std::cout<<otag<<"ERROR: The inverted matrix isnt valid! Something went wrong.."<<std::endl;
        exit(EXIT_FAILURE);

    }


    //check if the matrix is positive, semi-definite;
    bool is_small_negative_eigenvalue = false;
    double tolerence_positivesemi = 1e-5;


    TMatrixDEigen eigen (M);
    TVectorD eigen_values = eigen.GetEigenValuesRe();


    for(int i=0; i< eigen_values.GetNoElements(); i++){
        if(eigen_values(i)<0){
            is_small_negative_eigenvalue = true;
            if(fabs(eigen_values(i))> tolerence_positivesemi ){
                std::cout<<otag<<" covariance matrix contains (at least one)  negative eigenvalue: "<<eigen_values(i)<<std::endl;
                M.Print();
                exit(EXIT_FAILURE);
            }
        }
    }


    if(is_small_negative_eigenvalue){
        if(is_verbose)  std::cout<<otag<<"Covariance matrix is (allmost) positive semi-definite. It did contain small negative values of absolute value <= :"<<tolerence_positivesemi<<std::endl;
    }else{
        if(is_verbose)  std::cout<<otag<<"Covariance matrix is also positive semi-definite."<<std::endl;
    }


    return McI;

}

/*
TMatrixT<double> SBNchi::InvertMatrix(TMatrixT<double> &M){
    //This is an older simplier 
    double invdet=0;

    TMatrixT<double> McI(M.GetNrows(),M.GetNrows());
    McI.Zero();

    if(is_verbose) std::cout<<otag<<" About to do a SVD decomposition"<<std::endl;
    TDecompSVD svd(M);

    if (!svd.Decompose()){
        std::cout<<otag<<" (InvertMatrix) Decomposition failed, matrix not symettric?, has nans?" << std::endl;
        std::cout<<otag<<"ERROR: The matrix to invert failed a SVD decomp!"<<std::endl;

        for(int i=0; i< M.GetNrows(); i++){
            for(int j=0; j< M.GetNrows(); j++){
                std::cout<<i<<" "<<j<<" "<<M(i,j)<<std::endl;
            }
        }

        exit(EXIT_FAILURE);

    } else {
        McI = svd.Invert();
    }
    if( !McI.IsValid()){
        std::cout<<otag<<"ERROR: The inverted matrix isnt valid! Something went wrong.."<<std::endl;
        exit(EXIT_FAILURE);

    }


    return McI;

}
*/



/**************************************************************************
 *			Misc
 * ************************************************************************/

int SBNchi::FillCovarianceMatrix(TMatrixT<double>*in){
    in->ResizeTo(num_bins_total,num_bins_total);
    for(int i=0; i<num_bins_total;i++){
        for(int j=0; j<num_bins_total;j++){
            (*in)(i,j) = matrix_systematics(i,j);
        }
    }

    return 0;
}


int SBNchi::FillCollapsedCovarianceMatrix(TMatrixT<double>*in){
    in->ResizeTo(num_bins_total_compressed,num_bins_total_compressed) ;
    for(int i=0; i<num_bins_total_compressed;i++){
        for(int j=0; j<num_bins_total_compressed;j++){
            (*in)(i,j) = vec_matrix_collapsed.at(i).at(j) - (i==j? core_spectrum.collapsed_vector.at(i) : 0.0 );
        }
    }

    return 0;
}


int SBNchi::FillCorrelationMatrix(TMatrixT<double>*in){
    in->ResizeTo(num_bins_total,num_bins_total) ;
    for(int i=0; i<num_bins_total;i++){
        for(int j=0; j<num_bins_total;j++){
            double val = matrix_systematics(i,j);
            double p1 = sqrt(matrix_systematics(j,j)); 
            double p2 = sqrt(matrix_systematics(i,i)); 
            (*in)(i,j) = val/(p1*p2);
        }
    }
    return 0;
}

int SBNchi::FillCollapsedCorrelationMatrix(TMatrixT<double>*in){
    in->ResizeTo(num_bins_total_compressed,num_bins_total_compressed) ;
    for(int i=0; i<num_bins_total_compressed;i++){
        for(int j=0; j<num_bins_total_compressed;j++){
            double val = vec_matrix_collapsed.at(i).at(j) - (i==j? core_spectrum.collapsed_vector.at(i) : 0.0 ); 
            double p1 = sqrt(vec_matrix_collapsed.at(j).at(j)- core_spectrum.collapsed_vector.at(j)); 
            double p2 = sqrt(vec_matrix_collapsed.at(i).at(i)- core_spectrum.collapsed_vector.at(i)); 
            (*in)(i,j) = val/(p1*p2);
        }
    }

    return 0;
}


int SBNchi::FillFractionalMatrix(TMatrixT<double>*in){
    in->ResizeTo(num_bins_total,num_bins_total) ;
    for(int i=0; i<num_bins_total;i++){
        for(int j=0; j<num_bins_total;j++){
            (*in)(i,j) = ( matrix_systematics(i,j))/(core_spectrum.full_vector.at(i)*core_spectrum.full_vector.at(j));
        }
    }

    return 0;
}


int SBNchi::FillCollapsedFractionalMatrix(TMatrixT<double>*in){
    in->ResizeTo(num_bins_total_compressed,num_bins_total_compressed) ;
    for(int i=0; i<num_bins_total_compressed;i++){
        for(int j=0; j<num_bins_total_compressed;j++){
            (*in)(i,j) = ( vec_matrix_collapsed.at(i).at(j) - (i==j? core_spectrum.collapsed_vector.at(i) : 0.0 ) )/(core_spectrum.collapsed_vector.at(i)*core_spectrum.collapsed_vector.at(j));
        }
    }

    return 0;
}


TMatrixT<double> * SBNchi::GetCollapsedMatrix(){
    TMatrixT<double> * tmp = new TMatrixT<double>(num_bins_total_compressed,num_bins_total_compressed);
    for(int i=0; i<num_bins_total_compressed;i++){
        for(int j=0; j<num_bins_total_compressed;j++){
            (*tmp)(i,j) = vec_matrix_collapsed.at(i).at(j);
        }
    }

    return tmp;
}




void SBNchi::FakeFillMatrix(TMatrixT <double> &M){
    //Fills a square matrix of dim matrix_size with random numbers for now.
    std::uniform_real_distribution<double> dist(0,1);

    int matrix_size=M.GetNrows();
    if(M.GetNrows()!=M.GetNcols()){std::cout<<"#ERROR: not a square matrix!"<<std::endl;}
    for(int i=0; i<matrix_size; i++){
        for (int j = i;j<matrix_size;j++){
            M(i,j)= dist(*rangen_twister);
            M(j,i)=M(i,j);
        }
    }
    return ;
}


std::vector<std::vector<double >> SBNchi::TMatrixDToVector(TMatrixT <double > Min)
{
    int dimension =  Min.GetNrows();

    std::vector<std::vector<double >>  ans(dimension, std::vector<double>(dimension));

    for(int i = 0; i< dimension; i++){
        for(int k = 0; k< dimension; k++){
            ans[i][k]=Min(i,k);
            if(ans[i][k]==-0){
                ans[i][k]=0;
            }
        }
    }
    return ans;
}

void SBNchi::FillStatsMatrix(TMatrixT <double> &M, std::vector<double> diag){
    int matrix_size = M.GetNrows();

    if(matrix_size != diag.size()){std::cerr<<"#ERROR: FillStatsMatrix, matrix not equal to diagonal"<<std::endl;}
    if(M.GetNrows()!=M.GetNcols()){std::cerr<<"#ERROR: not a square matrix!"<<std::endl;}

    M.Zero();

    for(int i=0; i<matrix_size; i++)
    {
        M(i,i) = diag.at(i);
    }

    return ;
}


TMatrixT<double> SBNchi::FillSystMatrix(const TMatrixT<double>& frac_covar, const std::vector<double>& full, const std::vector<double>& full_err){
	return FillSystMatrix(frac_covar, full, full_err, false);
}


//return a full or collapsed systematic covariance matrix
TMatrixT<double> SBNchi::FillSystMatrix(const TMatrixT<double>& frac_covar, const std::vector<double>& full, const std::vector<double>& full_err, bool do_collapse){

	int matrix_size = frac_covar.GetNcols();
	if(matrix_size != full.size()){
		std::cout << "ERROR: FillSystMatrix, matrix has diffrent size as spectrum: "<< matrix_size<< " vs " << full.size() << std::endl;
		exit(EXIT_FAILURE);
	}

	TMatrixT<double> full_syst(matrix_size, matrix_size);
	for(int i=0; i< matrix_size ; i++){
		for(int j=0; j<matrix_size; j++){
			if(std::isnan(frac_covar(i,j))) full_syst(i,j) = 0.0;
			else full_syst(i,j) = frac_covar(i,j)*full[i]*full[j];

			if((i==j) && (!is_stat_only)) full_syst(i,j) += pow(full_err[i], 2.0);
		}
	}
	if(do_collapse == true){
		TMatrixT<double> collapsed_syst(num_bins_total_compressed, num_bins_total_compressed);
		CollapseModes(full_syst, collapsed_syst);
		return collapsed_syst;
	}else{
		return full_syst;
	}
}


TMatrixT<double> SBNchi::FillSystematicsFromXML(){
    return FillSystematicsFromXML(correlation_matrix_rootfile, correlation_matrix_name);
}


TMatrixT<double > SBNchi::FillSystematicsFromXML(std::string rootname, std::string matname){
    //Pretty much obsolete now, should fill directly really.
    std::cout<<"SBNchi::FillSystematicsFromXML || filling from "<<rootname<<std::endl;

    TMatrixT<double> temp2(num_bins_total,num_bins_total);
    TFile *fm= new TFile(rootname.c_str());

    TMatrixT<float> * temp = (TMatrixT <float>* )fm->Get(matname.c_str());
    //TMatrixT<double> * temp = (TMatrixT <double>* )fm->Get(matname.c_str());

    std::vector<std::vector<double>> mcont;

    for(int p:used_bins){
        std::vector<double> tvec;
        for(int u:used_bins){
            tvec.push_back( (*temp)(p,u) );
        }
        mcont.push_back(tvec);
    }

    for(int i =0; i<num_bins_total; i++)
    {
        for(int j =0; j<num_bins_total; j++)
        {
            temp2(i,j)=mcont[i][j];
        }
    }
    delete temp; temp  = nullptr;

    std::cout<<"SBNchi::FillSystematicsFromXML || loaded with dim : "<<temp2.GetNcols()<<" "<<temp2.GetNrows()<<std::endl;

    fm->Close();
    delete fm; fm = nullptr;

    if(temp2.IsSymmetric()){
        if(is_verbose)std::cout<<"Inputted fracCov covariance matrix is symmetric"<<std::endl;
    }else{
        std::cerr<<"ERROR: SBNchi::FillSystematicsFromXML, matrix_systematics input is not symmetric!"<<std::endl;
        //exit(EXIT_FAILURE);
    }

    return temp2;

}





TH2D* SBNchi::GetChiogram(){
    TH2D *tmp = new TH2D("chi-o-gram","chi-o-gram",num_bins_total_compressed,0, num_bins_total_compressed ,num_bins_total_compressed,0, num_bins_total_compressed);

    for(int i =0; i<num_bins_total_compressed; i++){
        for(int j =0; j<num_bins_total_compressed; j++){

            tmp->SetBinContent(i+1, j+1, vec_last_calculated_chi.at(i).at(j));
        }
    }

    return tmp;
}

int SBNchi::PrintMatricies(std::string tag){
    TFile* fout = new TFile(("SBNfit_collapsed_matrix_plots_"+tag+".root").c_str(),"recreate");
    fout->cd();

    gStyle->SetOptStat(0);

    TMatrixD full, frac, corr;
    this->FillCollapsedCovarianceMatrix(&full);
    this->FillCollapsedFractionalMatrix(&frac);
    this->FillCollapsedCorrelationMatrix(&corr);

    corr.Write("collapsed_correlation");

    gStyle->SetPalette(kLightTemperature);

    frac.Write("collapsed_fractional_covariance");
    TH2D h2_frac(frac);
    //h2_frac.Write();
    h2_frac.SetName("frac");
    TCanvas *c_frac = new TCanvas("collapsed fractional covariance matrix");
    TPad* p_frac = (TPad*)c_frac->cd();
    //p_frac->SetLogz();
    c_frac->SetFixedAspectRatio();
    h2_frac.Draw("colz");
    h2_frac.SetTitle("Collapsed fractional covariance matrix");
    h2_frac.GetXaxis()->SetTitle("Reco Bin i");
    h2_frac.GetYaxis()->SetTitle("Reco Bin j");
    if(m_cmin !=-999 || m_cmax !=-999)   h2_frac.GetZaxis()->SetRangeUser(m_cmin,m_cmax);

    c_frac->SetRightMargin(0.150);

    int use_frac =0;

    for(int im =0; im<num_modes; im++){
        for(int id =0; id<num_detectors; id++){
            for(int ic = 0; ic < num_channels; ic++){
                TLine *lv = new TLine(0, num_bins.at(ic)+use_frac, num_bins_total_compressed, num_bins.at(ic)+use_frac);
                TLine *lh = new TLine(num_bins.at(ic)+use_frac,0, num_bins.at(ic)+use_frac, num_bins_total_compressed);
                lv->SetLineWidth(1.5);
                lh->SetLineWidth(1.5);
                use_frac+=num_bins.at(ic);
                lv->Draw();
                lh->Draw();
            }
        }
    }
    c_frac->Write();
    c_frac->SaveAs(("SBNfit_collapsed_fractional_covariance_"+tag+".SBNplot.pdf").c_str(),"pdf");

    for(int i=0; i<h2_frac.GetNbinsX(); i++){
        std::cout<<"Collapsed Frac "<<i<<" "<<sqrt(h2_frac.GetBinContent(i+1,i+1))*100.0<<std::endl;
    }

    for(int i=0; i<core_spectrum.collapsed_vector.size(); i++){
        std::cout<<sqrt(core_spectrum.collapsed_vector.at(i))/core_spectrum.collapsed_vector.at(i)*100.0<<" ";
    }std::cout<<std::endl;


    gStyle->SetPalette(kLightTemperature);

    full.Write("collapsed_covariance");
    TH2D h2_full(full);
    h2_full.SetName("full");
    TCanvas *c_full = new TCanvas("collapsed covariance matrix");
    c_full->cd();
    c_full->SetFixedAspectRatio();
    h2_full.Draw("colz");
    h2_full.SetTitle("Collapsed covariance matrix");
    h2_full.GetXaxis()->SetTitle("Reco Bin i");
    h2_full.GetYaxis()->SetTitle("Reco Bin j");

    c_full->SetRightMargin(0.150);

    int use_full =0;
    for(int im =0; im<num_modes; im++){
        for(int id =0; id<num_detectors; id++){
            for(int ic = 0; ic < num_channels; ic++){
                TLine *lv = new TLine(0, num_bins.at(ic)+use_full, num_bins_total_compressed, num_bins.at(ic)+use_full);
                TLine *lh = new TLine(num_bins.at(ic)+use_full,0, num_bins.at(ic)+use_full, num_bins_total_compressed);
                lv->SetLineWidth(1.5);
                lh->SetLineWidth(1.5);
                use_full+=num_bins.at(ic);
                lv->Draw();
                lh->Draw();
            }
        }
    }
    c_full->Write();
    c_full->SaveAs(("SBNfit_collapsed_covariance_"+tag+".SBNplot.pdf").c_str(),"pdf");

    gStyle->SetPalette(kLightTemperature);
    TH2D h2_corr(corr);
    h2_corr.SetName("corr");
    //h2_corr.Write();
    TCanvas *c_corr = new TCanvas("collapsed correlation matrix");
    c_corr->cd();
    c_corr->SetFixedAspectRatio();
    h2_corr.Draw("colz");
    h2_corr.SetTitle("Collapsed correlation matrix");
    h2_corr.GetXaxis()->SetTitle("Reco Bin i");
    h2_corr.GetYaxis()->SetTitle("Reco Bin j");
    h2_corr.GetZaxis()->SetRangeUser(-1,1);
    c_corr->SetRightMargin(0.150);

    int use_corr =0;
    for(int im =0; im<num_modes; im++){
        for(int id =0; id<num_detectors; id++){
            for(int ic = 0; ic < num_channels; ic++){
                TLine *lv = new TLine(0, num_bins.at(ic)+use_corr, num_bins_total_compressed, num_bins.at(ic)+use_corr);
                TLine *lh = new TLine(num_bins.at(ic)+use_corr,0, num_bins.at(ic)+use_corr, num_bins_total_compressed);
                lv->SetLineWidth(1.5);
                lh->SetLineWidth(1.5);
                use_corr+=num_bins.at(ic);
                lv->Draw();
                lh->Draw();
            }
        }
    }
    c_corr->Write();
    c_corr->SaveAs(("SBNfit_collapsed_correlation_"+tag+".SBNplot.pdf").c_str(),"pdf");



    TCanvas *chiogram = new TCanvas("Chi-o-gram","Chi-o-gram");
    chiogram->cd();
    chiogram->SetFixedAspectRatio();

    TH2D * h_chiogram = (TH2D*)this->GetChiogram();
    h_chiogram->Draw("colz");
    h_chiogram->SetTitle("Reco Bin i");
    h_chiogram->SetTitle("Reco Bin j");

    int use_chio =0;
    for(int im =0; im<num_modes; im++){
        for(int id =0; id<num_detectors; id++){
            for(int ic = 0; ic < num_channels; ic++){
                TLine *lv = new TLine(0, num_bins.at(ic)+use_chio, num_bins_total_compressed, num_bins.at(ic)+use_chio);
                TLine *lh = new TLine(num_bins.at(ic)+use_chio,0, num_bins.at(ic)+use_chio, num_bins_total_compressed);
                lv->SetLineWidth(1.5);
                lh->SetLineWidth(1.5);
                use_chio+=num_bins.at(ic);
                lv->Draw();
                lh->Draw();
            }
        }
    }
    chiogram->Write();

    gStyle->SetPalette(kBird);

    TMatrixD ufull, ufrac, ucorr;
    this->FillCovarianceMatrix(&ufull);
    this->FillFractionalMatrix(&ufrac);
    this->FillCorrelationMatrix(&ucorr);

    plot_one(ufull,"SBNfit_uncollapsed_covariance_matrix_"+tag,fout,true,false,false);
    plot_one(ufrac,"SBNfit_uncollapsed_fractional_covariance_matrix"+tag,fout,true,false,false);
    plot_one(ucorr,"SBNfit_uncollapsed_correlation_matrix"+tag,fout,true,false,true);

    fout->Close();
    return 0;
}
int SBNchi::plot_one(TMatrixD matrix, std::string tag, TFile *fin, bool plot_pdf, bool indiv, bool is_corr){
    fin->cd();
    if(indiv){
        TDirectory *individualDir = fin->GetDirectory("individualDir"); 
        if (!individualDir) { 
            individualDir = fin->mkdir("individualDir");       
        }
        fin->cd(); 
        individualDir->cd();
    }
    if(is_corr) gStyle->SetPalette(kLightTemperature);


    if(is_corr){
        for(int i=0; i<matrix.GetNrows();i++){
            for(int j=0; j<matrix.GetNrows();j++){
                double val = matrix(i,j);
                if(val!=val || isinf(val) || std::isnan(val)){
                    matrix(i,j)=0.0;
                    matrix(i,i)=1.0;
                    matrix(j,j)=1.0;
                }
            }

        }
    }

    TH2D h2_full(matrix);
    h2_full.SetName((tag+"_th2d").c_str());
    TCanvas *c_full = new TCanvas((tag+"_canvas").c_str());
    TPad *p_full = (TPad*)c_full->cd();
    c_full->SetFixedAspectRatio();
    h2_full.Draw("colz");
    h2_full.SetTitle(tag.c_str());
    h2_full.GetXaxis()->SetTitle("Global Bin Number");
    h2_full.GetYaxis()->SetTitle(" ");
    h2_full.GetYaxis()->SetLabelSize(0);
    //p_full->SetLogz();
    if(is_corr){
        h2_full.GetZaxis()->SetRangeUser(0.4,1);
    }
    else{
        h2_full.GetZaxis()->SetRangeUser(-0.25,0.25);
    }    

    c_full->SetFrameFillColor(kWhite);
    c_full->SetFillColor(kWhite);
    p_full->SetFillColor(kWhite);


    c_full->SetRightMargin(0.150);
    c_full->SetLeftMargin(0.10);//0.250
    c_full->SetTopMargin(0.10);
    int use_full =0;

    double percent_left = 0.05;//0.15
    double nice_shift = num_bins_total*0.02;

    for(int im =0; im<num_modes; im++){
        for(int id =0; id<num_detectors; id++){
            for(int ic = 0; ic < num_channels; ic++){
                for(int isc = 0; isc < num_subchannels.at(ic); isc++){


                    std::string mode_det = mode_names[im] +" " +detector_names[id];
                    std::string chan_sub = channel_names[ic]+" "+subchannel_names[ic][isc];


                    TText * tmd = new TText(-num_bins_total*percent_left*0.15, use_full+nice_shift*0.5, (mode_det+" "+chan_sub).c_str() );

                    //TText * tmd = new TText(use_full*1.05, num_bins_total*1.015, chan_sub.c_str());
                    //TText * tcs = new TText(use_full*1.05, num_bins_total*1.055, mode_det.c_str());
                    tmd->SetTextColor(kBlack);
                    //tcs->SetTextColor(kBlack);
                    tmd->SetTextSize(0.03);
                    tmd->SetTextAlign(31);
                    //tcs->SetTextSize(0.03);

                    //dont plot names ta the moment
                    //tmd->Draw();
                    //tcs->Draw();


                    /*
                       TText * tlow_bin = new TText(-num_bins_total*percent_left, use_full+nice_shift*0.5, sbnfit_to_string_prec(bin_edges[ic].front(),0).c_str());
                       TText * thigh_bin = new TText(-num_bins_total*percent_left, (use_full+num_bins[ic])-nice_shift*1.4, sbnfit_to_string_prec(bin_edges[ic].back(),0).c_str());
                       tlow_bin->SetTextSize(0.02);
                       thigh_bin->SetTextSize(0.02);
                       tlow_bin->Draw();
                       thigh_bin->Draw();

                       TText * tunit = new TText(-num_bins_total*percent_left, use_full+0.5*num_bins[ic], channel_units[ic].c_str());
                       tunit->SetTextSize(0.03);
                       tunit->Draw();
                       */

                    if(isc<num_subchannels[ic]-1){
                        TLine *lscv = new TLine(-num_bins_total*percent_left, num_bins.at(ic)+use_full, num_bins_total, num_bins.at(ic)+use_full);
                        TLine *lsch = new TLine(num_bins.at(ic)+use_full,0, num_bins.at(ic)+use_full, num_bins_total*1.045);
                        lscv->SetLineWidth(3);
                        lsch->SetLineWidth(3);
                        lscv->SetLineColor(kRed);
                        lsch->SetLineColor(kRed);
                        lscv->SetLineStyle(9);
                        lsch->SetLineStyle(9);

                        //Going to drop the little ones for now
                        //lscv->Draw();
                        //lsch->Draw();

                        use_full+=num_bins.at(ic);

                    }
                }
                TLine *lv = new TLine(-num_bins_total*percent_left, num_bins.at(ic)+use_full, num_bins_total, num_bins.at(ic)+use_full);
                TLine *lh = new TLine(num_bins.at(ic)+use_full,0, num_bins.at(ic)+use_full, num_bins_total*1.045);
                lv->SetLineWidth(2);
                lh->SetLineWidth(2);
                lv->SetLineColor(kBlack);
                lh->SetLineColor(kBlack);
                use_full+=num_bins.at(ic);
                lv->Draw();
                lh->Draw();

            }
        }
    }


    c_full->Write();
    if(plot_pdf) c_full->SaveAs((tag+".pdf").c_str(),"pdf");




    return 0;
}


int SBNchi::PerformCholoskyDecomposition(SBNspec *specin){
    specin->CalcFullVector();
    specin->CalcErrorVector();
    is_verbose=false;
    double tol = m_tolerance;

    std::cout<<" Starting Cholosky Decomp, tolderance is "<<tol<<" and sample_collapse "<<pseudo_from_collapsed<<std::endl;

    TMatrixD U  = matrix_fractional_covariance;

    for(int i =0; i<U.GetNcols(); i++)
    {
        for(int j =0; j<U.GetNrows(); j++)
        {
            if(std::isnan(U(i,j))){
                U(i,j) = 0;
            }else{
                U(i,j)= U(i,j)*specin->full_vector.at(i)*specin->full_vector.at(j);
            }
            if(i==j && !(is_stat_only) )
            {   
                U(i,j) += pow(specin->full_err_vector[i],2);
            }
        }
    }

    //New bit, do we collapse and sample from collapsed or full! Debugging April2020 Collab Meeting @ Zarkos info from MiniBooNE ERA
    int n_t = (pseudo_from_collapsed ? num_bins_total_compressed : num_bins_total);
    TMatrixT<double > U_use(n_t,n_t);
    if(pseudo_from_collapsed){
        CollapseModes(U, U_use);
    }else{
        U_use = U;
    }

    //Stats error are NOT added back in herebut treat them as Poisson later. Seems better
    //TMatrixDEigen eigen (U); // This was original, but caused a lot of "Error in <MakeSchurr>: too many iterations". Move to explicit symmetric matrix
    //Is this Really the best way to construct?!?
    // TMatrixDSym U_explicit_sym(n_t);
    // for(int i=0; i< n_t;i++){
    //   for(int j=i; j< n_t;j++){
    //     U_explicit_sym[i][j] = U(i,j);
    // }
    // }

    bool was_modified = false;

    //Get a Determinant, check things.
    double det = U_use.Determinant(); 
    std::cout<<"SBNchi::CholeskyDecomposition\t||\t Checking determinant of Covariance Matrix: U "<<det<<std::endl; 

    if(det < tol){
        std::cout<<"SBNchi::CholeskyDecomposition\t||\t This determinant is below tolerance of : "<<m_tolerance<<std::endl;
        std::cout<<"SBNchi::CholeskyDecomposition\t||\t Going to add this back to diagonal of full covariance matrix : "<<tol<<std::endl;
        for(int a =0; a<U_use.GetNcols(); a++){
            U_use(a,a) += tol;
            was_modified = true;
        }
        det = U_use.Determinant(); 
        std::cout<<"SBNchi::CholeskyDecomposition\t||\t The modified determinant is now: "<<det<<std::endl;
    }

    TMatrixDEigen eigen(U_use);
    TVectorD eigen_values = eigen.GetEigenValuesRe();
    TVectorD eigen_values_IM = eigen.GetEigenValuesIm();

    int n_zeros = 0;
    for(int i=0; i< eigen_values.GetNoElements(); i++){
        if(is_verbose)std::cout<<"SBNchi::CholeskyDecomposition\t||\t Eigenvalue "<<i<<" is Re: "<<eigen_values(i)<<" Im: "<<eigen_values_IM(i)<<std::endl;
        if(eigen_values(i)<tol){
            if(fabs(eigen_values(i)) > 0){
                if(is_verbose)std::cout<<"SBNchi::CholeskyDecomposition\t|| cov has a very small, < "<<tol<<" , negative eigenvalue. Adding it back to diagonal of : "<<eigen_values(i)<<std::endl;

                for(int a =0; a<U_use.GetNcols(); a++){
                    U_use(a,a) += eigen_values(i);
                    was_modified = true;
                }

            }else{
                std::cout<<"SBNchi::SampleCovariance\t|| 0 or negative eigenvalues! error: Value "<<eigen_values(i)<<" Tolerence "<<tol<<std::endl;
                //U_use.Print();
                //std::cout<<"Hmm"<<std::endl;
                //exit(EXIT_FAILURE);
            }
        }
    }


    //If everything is OK, lets pass this matrix back to SBNchi for use.
    if(was_modified){
        std::cout<<"We had to add on small diagonal terms to covariance matrix to Decompose it. Adding back to primary fractional covariance for consistency"<<std::endl;
        for(int i=0; i< n_t; i++){
            for(int j=0; j< n_t; j++){
                double mi = specin->full_vector.at(i)*specin->full_vector.at(j);
                double mcstat = 0;
                if(i==j) mcstat = pow(specin->full_err_vector[i],2);
                matrix_fractional_covariance(i,j) = ( mi==0   ?  0 :  (U_use(i,j) - mcstat  )/(mi)  );
            }
        }
        this->ReloadCoreSpectrum(specin);
    }


    //Seconndly attempt a Cholosky Decomposition
    TDecompChol * chol = new TDecompChol(U_use,tol);
    bool worked = chol->Decompose();

    if(!worked){
        std::cout<<"SBNchi::SampleCovariance\t|| Cholosky Decomposition Failed. Tolerance is set at "<<tol<<std::endl;
        //U_use.Print();
        std::cout<<"With eigens"<<std::endl;

        for(int i=0; i< eigen_values.GetNoElements(); i++){
            std::cout<<eigen_values(i)<<std::endl;
        }
        exit(EXIT_FAILURE);
    }

    TMatrixT<float> upper_trian(n_t,n_t);
    matrix_lower_triangular.ResizeTo(n_t,n_t);
    upper_trian = chol->GetU();
    matrix_lower_triangular = upper_trian;
    matrix_lower_triangular.T();

    vec_matrix_lower_triangular.resize(n_t, std::vector<float>(n_t));
    for(int i=0; i< n_t; i++){
        for(int j=0; j< n_t; j++){
            vec_matrix_lower_triangular[i][j] = matrix_lower_triangular[i][j];
            //            std::cout<<"Flormph "<<i<<" "<<j<<" "<<vec_matrix_lower_triangular[i][j]<<" "<<U_use(i,j)<<std::endl;
        }
    }


    //New Check, rebuild matrix and compare
    TMatrixDSym rebuild = chol->GetMatrix();
    TMatrixD rebuild2 = matrix_lower_triangular*upper_trian;

    for(int i=0; i< n_t; i++){
        for(int j=0; j< n_t; j++){
            double fd1 = fabs(rebuild(i,j)-U_use(i,j))/fabs(rebuild(i,j)+U_use(i,j));
            double fd2 = fabs(rebuild2(i,j)-U_use(i,j))/fabs(rebuild2(i,j)+U_use(i,j));
            if( fd1>m_tolerance || fd2>m_tolerance){
                //std::cout<<"ERROR the rebuilt matrix after Cholesky Decomp is not the same as the original."<<std::endl;
                //std::cout<<i<<" "<<j<<" Original: "<<U_use(i,j)<<" , Rebuild: "<<rebuild(i,j)<<" , Diff "<<fd1<<std::endl;
                //std::cout<<i<<" "<<j<<" Original: "<<U_use(i,j)<<" , Rebuild2: "<<rebuild2(i,j)<<" , Diff "<<fd2<<std::endl;
                //exit(EXIT_FAILURE);
            }
        }
    }
    std::cout<<"Rebuilt matrix is the same as input after Cholesky Decomp"<<std::endl;



    cholosky_performed = true;	
    delete chol; chol = nullptr;
    return 0;
}


TH1D SBNchi::SampleCovarianceVaryInput(SBNspec *specin, int num_MC, double maxchi){ 
    max_sample_chi_val = maxchi;
    std::vector<double>  tmp;
    return SampleCovarianceVaryInput(specin,num_MC,&tmp);
}

TH1D SBNchi::SampleCovarianceVaryInput(SBNspec *specin, int num_MC, std::vector<double> * chival){
    if(!cholosky_performed) this->PerformCholoskyDecomposition(specin); 

    float** a_vec_matrix_lower_triangular = new float*[num_bins_total];
    float** a_vec_matrix_inverted = new float*[num_bins_total_compressed];

    for(int i=0; i < num_bins_total; i++){
        a_vec_matrix_lower_triangular[i] = new float[num_bins_total];
    }

    for(int i=0; i < num_bins_total_compressed; i++){
        a_vec_matrix_inverted[i] = new float[num_bins_total_compressed];
    }

    for(int i=0; i < num_bins_total; i++){
        for(int j=0; j < num_bins_total; j++){
            a_vec_matrix_lower_triangular[i][j] = vec_matrix_lower_triangular[i][j]; 
        }
    }

    for(int i=0; i< num_bins_total_compressed; i++){
        for(int j=0; j< num_bins_total_compressed; j++){
            a_vec_matrix_inverted[i][j] = vec_matrix_inverted[i][j]; 
        }
    }

    float *a_specin = new float[num_bins_total];
    float *a_corein = new float[num_bins_total_compressed];

    for(int i=0; i< num_bins_total; i++){
        a_specin[i] = specin->full_vector[i];
    }

    for(int i=0; i< num_bins_total_compressed; i++) {
        a_corein[i] = core_spectrum.collapsed_vector[i];
    }

    TH1D ans("","",std::max(200,(int)max_sample_chi_val),0,max_sample_chi_val );
    //ans.GetXaxis()->SetCanExtend(kTRUE);
    is_verbose = false;

    std::vector<float> vec_chis (num_MC, 0.0);

    float* a_vec_chis  = (float*)vec_chis.data();
    int num_chival = chival->size();
    float* a_chival = new float[num_chival];

    int *nlower = new int[num_chival];
    for(int i=0; i< num_chival; i++){
        nlower[i]=0; 
        a_chival[i] = chival->at(i);
    }

    std::vector < float > gaus_sample_v(num_bins_total), sampled_fullvector_v(num_bins_total);
    std::vector<float> collapsed_v(num_bins_total_compressed, 0.0);

    float* gaus_sample = new float[num_bins_total];
    float* sampled_fullvector = new float[num_bins_total] ;
    float* collapsed = new float[num_bins_total_compressed];

    //  float gaus_sample[54];
    // float sampled_fullvector[54];
    //  float collapsed[38];

    //We will need a uniform dist and a Gaussian
    std::uniform_int_distribution<int> dist_int(0,pow(2,32));
    std::normal_distribution<float> dist_normal(0,1);

#ifdef USE_GPU
    unsigned long long seed[num_MC];
    unsigned long long seq = 0ULL;
    unsigned long long offset = 0ULL;
    curandState_t state;

    for(int i=0; i<num_MC; ++i) {
        seed[i] = dist_int(*rangen_twister);
    }
#endif

#ifdef USE_GPU
#pragma acc parallel loop  private(gaus_sample[:54],sampled_fullvector[:54],collapsed[:36],state) \
    copyin(this[0:1],							\
            a_specin[:num_bins_total],					\
            a_vec_matrix_lower_triangular[:num_bins_total][:num_bins_total],\ 
            a_corein[:num_bins_total_compressed],				\
            a_vec_matrix_inverted[:num_bins_total_compressed][:num_bins_total_compressed],	\
            seed[0:num_MC],						\
            a_chival[:num_chival],						\
            this->a_num_bins[:num_channels],				\
            this->a_num_subchannels[:num_channels])			\    
        copyout(a_vec_chis[:num_MC]) \
        copy(nlower[:num_chival])
#endif

        for(int i=0; i < num_MC;i++){

#ifdef USE_GPU
            unsigned long long seed_sd = seed[i];
            curand_init(seed_sd, seq, offset, &state);

            for(int a=0; a<num_bins_total; a++) {
                gaus_sample[a]= curand_normal(&state);
            }
#else
            for(int a=0; a<num_bins_total; a++) {
                gaus_sample[a]= dist_normal(*rangen_twister);
            }      
#endif

            for(int j = 0; j < num_bins_total; j++){
                sampled_fullvector[j] = a_specin[j];
                for(int k = 0; k < num_bins_total; k++){
                    sampled_fullvector[j] += a_vec_matrix_lower_triangular[j][k] * gaus_sample[k];
                }

                if(sampled_fullvector[j]<0) sampled_fullvector[j]=0.0;

                sampled_fullvector[j] = rangen->Poisson(sampled_fullvector[j]);
                //sampled_fullvector[j] = rangen->Poisson(a_specin[j]);
                //std::cout<<"P: "<<a_specin[j]<<" "<<sampled_fullvector[j]<<std::endl;
            }

            this->CollapseVectorStandAlone(sampled_fullvector, collapsed);

            a_vec_chis[i] = this->CalcChi(a_vec_matrix_inverted, a_corein, collapsed);
            //Just to get some pvalues that were asked for.

            for(int j=0; j< num_chival; j++){
#pragma acc atomic update
                if(a_vec_chis[i]>=a_chival[j]) nlower[j]++;
            }

        }

    is_verbose = false;


    for(int i=0; i<num_MC; i++){
        //       if (i<(int)1e3) 
        //         std::cout << "@i=" << a_vec_chis[i] << std::endl;
        ans.Fill(a_vec_chis[i]);
    }
    for(int n =0; n< num_chival; n++){
        chival->at(n) = nlower[n]/(double)num_MC;
    }



    delete[] a_corein;
    delete[] a_specin;
    delete[] nlower;

    for(int i=0; i < num_bins_total; i++){
        delete[] a_vec_matrix_lower_triangular[i];
    }

    for(int i=0; i < num_bins_total_compressed; i++){
        delete[] a_vec_matrix_inverted[i];  
    }

    delete[] a_vec_matrix_lower_triangular;
    delete[] a_vec_matrix_inverted;

    delete[] gaus_sample;
    delete[] sampled_fullvector;
    delete[] collapsed;

    return ans;
}


int SBNchi::CollapseVectorStandAlone(std::vector<double> * full_vector, std::vector<double> *collapsed_vector){
    for(int im = 0; im < num_modes; im++){
        for(int id =0; id < num_detectors; id++){
            int edge = id*num_bins_detector_block + num_bins_mode_block*im; // This is the starting index for this detector
            int out_edge = edge;
            int tmp_chan = 0;
            for(int ic = 0; ic < num_channels; ic++){
                int corner=edge;

                for(int j=0; j< num_bins[ic]; j++){

                    double tempval=0;

                    for(int sc = 0; sc < num_subchannels[ic]; sc++){
                        tempval += (*full_vector)[j+sc*num_bins[ic]+corner];
                        edge +=1;	//when your done with a channel, add on every bin you just summed
                    }
                    //we can size this vector beforehand and get rid of all push_back()

                    int collapsed_index = tmp_chan+out_edge;
                    (*collapsed_vector)[collapsed_index] = tempval;
                    tmp_chan++;
                }
            }
        }
    }


    return 0;
}

int SBNchi::CollapseVectorStandAlone(float* full_vector, float *collapsed_vector){

    //int tmp_num_bins[3] = {25,25,6};
    //int tmp_num_subchannels[3] = {2,1,1};

    int collapsed_index = 0;
    for(int im = 0; im < num_modes; im++){
        for(int id =0; id < num_detectors; id++){
            int edge = id*num_bins_detector_block + num_bins_mode_block*im; // This is the starting index for this detector
            int out_edge = edge;
            //            int chan = 0;
            for(int ic = 0; ic < num_channels; ic++){
                int corner=edge;

                for(int j=0; j< this->a_num_bins[ic]; j++){

                    float tempval=0;

                    for(int sc = 0; sc < this->a_num_subchannels[ic]; sc++){
                        tempval += full_vector[j+sc*this->a_num_bins[ic]+corner];
                        //tempval += full_vector[j+sc*tmp_num_bins[ic]+corner];
                        edge +=1;	//when your done with a channel, add on every bin you just summed
                    }
                    //we can size this vector beforehand and get rid of all push_back()

                    // int collapsed_index = chan+out_edge;
                    collapsed_vector[collapsed_index] = tempval;
                    collapsed_index++;
                }
            }
        }
    }


    return 0;
}


int SBNchi::CollapseVectorStandAlone(double* full_vector, double *collapsed_vector){

    //int tmp_num_bins[3] = {25,25,6};
    //int tmp_num_subchannels[3] = {2,1,1};

    for(int im = 0; im < num_modes; im++){
        for(int id =0; id < num_detectors; id++){
            int edge = id*num_bins_detector_block + num_bins_mode_block*im; // This is the starting index for this detector
            int out_edge = edge;
            int chan = 0;
            for(int ic = 0; ic < num_channels; ic++){
                int corner=edge;

                for(int j=0; j< this->a_num_bins[ic]; j++){

                    double tempval=0;

                    for(int sc = 0; sc < this->a_num_subchannels[ic]; sc++){
                        tempval += full_vector[j+sc*this->a_num_bins[ic]+corner];
                        //tempval += full_vector[j+sc*tmp_num_bins[ic]+corner];
                        edge +=1;	//when your done with a channel, add on every bin you just summed
                    }
                    //we can size this vector beforehand and get rid of all push_back()

                    int collapsed_index = chan+out_edge;
                    collapsed_vector[collapsed_index] = tempval;
                    chan++;
                }
            }
        }
    }
    return 0;
}

std::vector<float> SBNchi::GeneratePseudoExperiment(){

    core_spectrum.CollapseVector();
    if(!cholosky_performed || is_stat_only){
        std::cout<<"CHOL now"<<std::endl;
        PerformCholoskyDecomposition(&core_spectrum); 
    }

    int n_t =  (pseudo_from_collapsed ? num_bins_total_compressed : num_bins_total); 
    std::vector<float> sampled(n_t);
    is_verbose = false;

    std::vector<double> v_gaus(n_t, 0.0);
    for(int i=0; i< n_t; ++i){
        //v_gaus[i] = rangen->Gaus(0,1);
        //v_gaus[i] = rangen->Uniform(-1,1); 
        v_gaus[i] = (*m_dist_normal)(*rangen_twister); 
    }

    for(int i=0; i< n_t; ++i){
        sampled[i] = (pseudo_from_collapsed ? core_spectrum.collapsed_vector[i] : core_spectrum.full_vector[i]); 
        if(!is_stat_only){
            for(int j=0; j<n_t; ++j){
                sampled[i] += vec_matrix_lower_triangular[i][j]*v_gaus[j];
            }
        }
    }
    //Now poisson fluctuate the sampled spectrum
    for(int j=0; j<n_t; ++j){
        std::poisson_distribution<int> dist_pois(sampled[j]);
        sampled[j] = float(dist_pois(*rangen_twister));
    }

    std::vector<float> collapsed(num_bins_total_compressed,0.0);

    if(pseudo_from_collapsed){
        collapsed = sampled;
    }else{
        this->CollapseVectorStandAlone(&sampled[0], &collapsed[0]);
    }
    return collapsed;
}


std::vector<float> SBNchi::SampleCovariance(SBNspec *specin){
    if(!cholosky_performed) this->PerformCholoskyDecomposition(specin); 

    int n_t = specin->full_vector.size();

    TVectorT<float> u(n_t);
    for(int i=0; i<n_t; i++){
        u(i) = specin->full_vector[i];
    }

    std::normal_distribution<float> dist_normal(0,1);
    is_verbose = false;

    TVectorT<float> gaus_sample(n_t);
    TVectorT<float> multi_sample(n_t);
    for(int a=0; a<n_t; a++){
        gaus_sample(a) = dist_normal(*rangen_twister);	
    }

    multi_sample = u + matrix_lower_triangular*gaus_sample;

    std::vector<float> sampled_fullvector(n_t,0.0);
    for(int j=0; j<n_t; j++){
        sampled_fullvector[j] = multi_sample(j);
    }

    std::vector<float> collapsed(num_bins_total_compressed,0.0);
    this->CollapseVectorStandAlone(&sampled_fullvector[0], &collapsed[0]);

    return collapsed;
}


TH1D SBNchi::SamplePoissonVaryInput(SBNspec *specin, int num_MC, double maxchi){ 
    max_sample_chi_val = maxchi;
    std::vector<double>  tmp = {};
    return SamplePoissonVaryInput(specin,num_MC,&tmp);
}


//This one varies the input comparative spectrum, and as sucn has  only to calculate the matrix_systematics once
TH1D SBNchi::SamplePoissonVaryInput(SBNspec *specin, int num_MC, std::vector<double> *chival){

    float** a_vec_matrix_inverted = new float*[num_bins_total_compressed];

    for(int i=0; i < num_bins_total_compressed; i++){
        a_vec_matrix_inverted[i] = new float[num_bins_total_compressed];
    }
    for(int i=0; i< num_bins_total_compressed; i++){
        for(int j=0; j< num_bins_total_compressed; j++){
            a_vec_matrix_inverted[i][j] = vec_matrix_inverted[i][j]; 
        }
    }

    float *a_specin = new float[num_bins_total];
    float *a_corein = new float[num_bins_total_compressed];

    for(int i=0; i< num_bins_total; i++){
        a_specin[i] = specin->full_vector[i];
    }

    for(int i=0; i< num_bins_total_compressed; i++) {
        a_corein[i] = core_spectrum.collapsed_vector[i];
    }

    std::vector<float> vec_chis (num_MC, 0.0);

    float* a_vec_chis  = (float*)vec_chis.data();
    int num_chival = chival->size();
    float* a_chival = new float[num_chival];

    int *nlower = new int[num_chival];
    for(int i=0; i< num_chival; i++){
        nlower[i]=0; 
        a_chival[i]=chival->at(i); 
    }
    float* sampled_fullvector = new float[num_bins_total] ;
    float* collapsed = new float[num_bins_total_compressed];

    //So save the core one that we will sample for
    //ans.GetXaxis()->SetCanExtend(kTRUE);
    is_verbose = false;

    std::vector< std::poisson_distribution<int>> dist_pois;
    //std::vector< std::normal_distribution<float>> dist_pois;
    for(int j = 0; j < num_bins_total; j++){
        //for tesing purposes
        dist_pois.push_back(std::poisson_distribution<int>(a_specin[j])); 
        //dist_pois.push_back(std::normal_distribution<float>(a_specin[j],sqrt(a_specin[j])));
    }

    for(int i=0; i < num_MC;i++){

        for(int j = 0; j < num_bins_total; j++){

            //float p = dist_pois[j](*rangen_twister); 
            int p = dist_pois[j](*rangen_twister); 

            sampled_fullvector[j] =  (float)p;
            //std::cout<<"P: "<<a_specin[j]<<" "<<sampled_fullvector[j]<<" "<<p<<std::endl;
        }

        this->CollapseVectorStandAlone(sampled_fullvector, collapsed);
        a_vec_chis[i] = this->CalcChi(a_vec_matrix_inverted, a_corein, collapsed);

        for(int j=0; j< num_chival; j++){
            if(a_vec_chis[i]>=a_chival[j]) nlower[j]++;
        }

    }

    TH1D ans("","",std::max(200,(int)max_sample_chi_val),0,max_sample_chi_val);


    for(int i=0; i<num_MC; i++){
        ans.Fill(a_vec_chis[i]);
    }

    for(int n =0; n< num_chival; n++){
        chival->at(n) = nlower[n]/(double)num_MC;
    }

    is_verbose = true;

    delete[] a_corein;
    delete[] a_specin;
    delete[] nlower;

    for(int i=0; i < num_bins_total_compressed; i++){
        delete[] a_vec_matrix_inverted[i];  
    }

    delete[] a_vec_matrix_inverted;

    delete[] sampled_fullvector;
    delete[] collapsed;

    return ans;


}

TH1D SBNchi::SamplePoisson_NP(SBNspec *specin, SBNchi &chi_h0, SBNchi & chi_h1, int num_MC, double maxchi,int which_sample){ 
    max_sample_chi_val = maxchi;
    std::vector<double>  tmp = {};
    return SamplePoisson_NP(specin,chi_h0,chi_h1,num_MC,&tmp, which_sample);
}


TH1D SBNchi::SamplePoisson_NP(SBNspec *specin, SBNchi &chi_h0, SBNchi & chi_h1, int num_MC, std::vector<double> *chival, int which_sample){

    float** h0_vec_matrix_inverted = new float*[num_bins_total_compressed];
    float** h1_vec_matrix_inverted = new float*[num_bins_total_compressed];

    for(int i=0; i < num_bins_total_compressed; i++){
        h0_vec_matrix_inverted[i] = new float[num_bins_total_compressed];
        h1_vec_matrix_inverted[i] = new float[num_bins_total_compressed];
    }
    for(int i=0; i< num_bins_total_compressed; i++){
        for(int j=0; j< num_bins_total_compressed; j++){
            h0_vec_matrix_inverted[i][j] = chi_h0.vec_matrix_inverted[i][j]; 
            h1_vec_matrix_inverted[i][j] = chi_h1.vec_matrix_inverted[i][j]; 
        }
    }

    float *a_specin = new float[num_bins_total];

    float *h0_corein = new float[num_bins_total_compressed];
    float *h1_corein = new float[num_bins_total_compressed];


    for(int i=0; i< num_bins_total; i++){
        a_specin[i] = specin->full_vector[i];
    }

    for(int i=0; i< num_bins_total_compressed; i++) {
        h0_corein[i] = chi_h0.core_spectrum.collapsed_vector[i];
        h1_corein[i] = chi_h1.core_spectrum.collapsed_vector[i];
    }

    std::vector<float> vec_chis (num_MC, 0.0);

    float* a_vec_chis  = (float*)vec_chis.data();
    int num_chival = chival->size();
    float* a_chival = new float[num_chival];

    int *nlower = new int[num_chival];
    for(int i=0; i< num_chival; i++){
        nlower[i]=0; 
        a_chival[i]=chival->at(i); 
    }
    float* sampled_fullvector = new float[num_bins_total] ;
    float* collapsed = new float[num_bins_total_compressed];

    float min_delta_chi = 99999;

    //So save the core one that we will sample for
    //ans.GetXaxis()->SetCanExtend(kTRUE);
    is_verbose = false;

    std::vector< std::poisson_distribution<int>> dist_pois;
    //std::vector< std::normal_distribution<float>> dist_pois;
    for(int j = 0; j < num_bins_total; j++){
        //for tesing purposes
        dist_pois.push_back(std::poisson_distribution<int>(a_specin[j])); 
        //dist_pois.push_back(std::normal_distribution<float>(a_specin[j],sqrt(a_specin[j])));
    }

    std::cout<<otag<<" Starting to generate "<<num_MC<<" pseudo universes according to poisson distribution"<<std::endl;
    for(int i=0; i < num_MC;i++){

        if(which_sample==0){//Poisson Mode
            for(int j = 0; j < num_bins_total; j++){
                //float p = dist_pois[j](*rangen_twister); 
                int p = dist_pois[j](*rangen_twister); 
                sampled_fullvector[j] =  (float)p;
                //std::cout<<"P: "<<a_specin[j]<<" "<<sampled_fullvector[j]<<" "<<p<<std::endl;
            }

            this->CollapseVectorStandAlone(sampled_fullvector, collapsed);

        }else if(which_sample==1){//Covariance Sampling

            std::vector<float> exp  = this->GeneratePseudoExperiment();

            for(int j = 0; j < num_bins_total_compressed; j++){
                collapsed[j] = exp[j];
            }
        }

	//vector "collapsed" is varied h1 spectrum
        float val_chi_h0  = chi_h0.CalcChi(h0_vec_matrix_inverted, h0_corein, collapsed);
        float val_chi_h1  = chi_h1.CalcChi(h1_vec_matrix_inverted, h1_corein, collapsed);

        a_vec_chis[i] = val_chi_h0-val_chi_h1;

        //std::cout<<" "<<val_chi_h0<<" "<<val_chi_h1<<std::endl;

        /*
           std::cout<<" Pseudo Data "<<std::endl;
           for(int m=0; m<num_bins_total_compressed; m++) std::cout<<" "<<collapsed[m]<<" ";
           std::cout<<std::endl;
           std::cout<<" H0 "<<val_chi_h0<<std::endl;
           for(int m=0; m<num_bins_total_compressed; m++) std::cout<<" "<<h0_corein[m]<<" ";
           std::cout<<std::endl;
           std::cout<<" H1 "<<val_chi_h1<<std::endl;
           for(int m=0; m<num_bins_total_compressed; m++) std::cout<<" "<<h1_corein[m]<<" ";
           std::cout<<std::endl;
           std::cout<<" H0_ERR "<<std::endl;
           for(int m=0; m<num_bins_total_compressed; m++) std::cout<<" "<<h0_vec_matrix_inverted[m][m]<<" ";
           std::cout<<std::endl;
           std::cout<<" H1_ERR "<<std::endl;
           for(int m=0; m<num_bins_total_compressed; m++) std::cout<<" "<<h1_vec_matrix_inverted[m][m]<<" ";
           std::cout<<std::endl;
           std::cout<<"Delta H0-H1 "<<a_vec_chis[i]<<std::endl;
           */
        if(a_vec_chis[i] < min_delta_chi) min_delta_chi = a_vec_chis[i];

        for(int j=0; j< num_chival; j++){
            if(a_vec_chis[i]>=a_chival[j]) nlower[j]++;
        }

    }

    TH1D ans("","",std::max(200,(int)max_sample_chi_val),min_delta_chi,max_sample_chi_val);


    for(int i=0; i<num_MC; i++){
        ans.Fill(a_vec_chis[i]);
    }

    for(int n =0; n< num_chival; n++){
        chival->at(n) = nlower[n]/(double)num_MC;
    }

    is_verbose = true;

    delete[] h1_corein;
    delete[] h0_corein;
    delete[] a_specin;
    delete[] nlower;

    for(int i=0; i < num_bins_total_compressed; i++){
        delete[] h1_vec_matrix_inverted[i];  
        delete[] h0_vec_matrix_inverted[i];  
    }

    delete[] h1_vec_matrix_inverted;
    delete[] h0_vec_matrix_inverted;

    delete[] sampled_fullvector;
    delete[] collapsed;

    return ans;


}


std::vector<CLSresult> SBNchi::Mike_NP(SBNspec *specin, SBNchi &chi_h0, SBNchi & chi_h1, int num_MC, int which_sample, int id){

    std::vector<CLSresult> v_results(5);

    float** h0_vec_matrix_inverted = new float*[num_bins_total_compressed];
    float** h1_vec_matrix_inverted = new float*[num_bins_total_compressed];

    for(int i=0; i < num_bins_total_compressed; i++){
        h0_vec_matrix_inverted[i] = new float[num_bins_total_compressed];
        h1_vec_matrix_inverted[i] = new float[num_bins_total_compressed];
    }
    for(int i=0; i< num_bins_total_compressed; i++){
        for(int j=0; j< num_bins_total_compressed; j++){
            h0_vec_matrix_inverted[i][j] = chi_h0.vec_matrix_inverted[i][j]; 
            h1_vec_matrix_inverted[i][j] = chi_h1.vec_matrix_inverted[i][j]; 
        }
    }

    float *a_specin = new float[num_bins_total];

    float *h0_corein = new float[num_bins_total_compressed];
    float *h1_corein = new float[num_bins_total_compressed];


    for(int i=0; i< num_bins_total; i++){
        a_specin[i] = specin->full_vector[i];
    }

    for(int i=0; i< num_bins_total_compressed; i++) {
        h0_corein[i] = chi_h0.core_spectrum.collapsed_vector[i];
        h1_corein[i] = chi_h1.core_spectrum.collapsed_vector[i];
    }

    std::vector<float> vec_chis (num_MC, 0.0);
    std::vector<float> vec_pois (num_MC, 0.0);
    std::vector<float> vec_cnp (num_MC, 0.0);
    std::vector<float> vec_h0 (num_MC, 0.0);
    std::vector<float> vec_h1 (num_MC, 0.0);

    float* a_vec_chis  = (float*)vec_chis.data();
    float* a_vec_pois  = (float*)vec_pois.data();
    float* a_vec_cnp  = (float*)vec_cnp.data();
    float* a_vec_h0  = (float*)vec_h0.data();
    float* a_vec_h1  = (float*)vec_h1.data();


    float* sampled_fullvector = new float[num_bins_total] ;
    float* collapsed = new float[num_bins_total_compressed];
    float* cv_collapsed = new float[num_bins_total_compressed];

    is_verbose = false;
    this->CollapseVectorStandAlone(a_specin, cv_collapsed);


    std::cout<<otag<<" Starting to generate "<<num_MC<<" Pseduo-Experiments."<<std::endl;
    for(int i=0; i < num_MC;i++){

        //Generate our spectra
        if(which_sample==0){//Poisson Mode
            for(int j = 0; j < num_bins_total_compressed; j++){
                std::poisson_distribution<int> dist_pois(cv_collapsed[j]);
                collapsed[j] = float(dist_pois(*rangen_twister));
            }
        }else if(which_sample==1){//Covariance Sampling
            auto exp = this->GeneratePseudoExperiment();
            for(int j=0; j < num_bins_total_compressed; j++){
                collapsed[j] = exp[j];
            }
        }

        //Base Default Chi
        float val_chi_h0  = chi_h0.CalcChi_Pearson(h0_corein, collapsed);
        float val_chi_h1  = chi_h1.CalcChi_Pearson(h1_corein, collapsed);

        //Poisson Log Likli
        float val_pois_h0  = chi_h0.PoissonLogLiklihood(h0_corein, collapsed);
        float val_pois_h1  = chi_h1.PoissonLogLiklihood(h1_corein, collapsed);

        //CNP, going to need to recalculate and reinvert.
        float val_cnp_h0  = chi_h0.CalcChi_CNP(h0_corein, collapsed);
        float val_cnp_h1  = chi_h1.CalcChi_CNP(h1_corein, collapsed);

        a_vec_chis[i] = val_chi_h0 - val_chi_h1;
        a_vec_pois[i] = val_pois_h0 - val_pois_h1;
        a_vec_cnp[i]  = val_cnp_h0 - val_cnp_h1;
        a_vec_h0[i] = val_chi_h0;
        a_vec_h1[i] = val_chi_h1;

        if(i%1000==0) std::cout<<"Pseudo-Experiment: "<<i<<"/"<<num_MC<<" DeltaChi: "<<a_vec_chis[i]<<" PoisLogLiki: "<<a_vec_pois[i]<<" CNP_chi: "<<a_vec_cnp[i]<<std::endl;

        if(a_vec_chis[i] < v_results[0].m_min_value) v_results[0].m_min_value = a_vec_chis[i];
        if(a_vec_pois[i] < v_results[1].m_min_value) v_results[1].m_min_value = a_vec_pois[i];
        if(a_vec_cnp[i]  < v_results[2].m_min_value) v_results[2].m_min_value = a_vec_cnp[i];
        if(val_chi_h0  < v_results[3].m_min_value) v_results[3].m_min_value = val_chi_h0;
        if(val_chi_h1  < v_results[4].m_min_value) v_results[4].m_min_value = val_chi_h1;

        if(a_vec_chis[i] > v_results[0].m_max_value) v_results[0].m_max_value = a_vec_chis[i];
        if(a_vec_pois[i] > v_results[1].m_max_value) v_results[1].m_max_value = a_vec_pois[i];
        if(a_vec_cnp[i]  > v_results[2].m_max_value) v_results[2].m_max_value = a_vec_cnp[i];
        if(val_chi_h0  > v_results[3].m_max_value) v_results[3].m_max_value = val_chi_h0;
        if(val_chi_h1  > v_results[4].m_max_value) v_results[4].m_max_value = val_chi_h1;

    }

    for(int i=0; i<3;i++){
        std::cout<<"Res "<<i<<" "<<v_results[i].m_max_value<<" "<<v_results[i].m_min_value<<std::endl;
    }

    TH1D ans0(("0"+std::to_string(id)).c_str(),("0"+std::to_string(id)).c_str(),std::max(200,(int)v_results[0].m_max_value),v_results[0].m_min_value,v_results[0].m_max_value);
    TH1D ans1(("1"+std::to_string(id)).c_str(),("1"+std::to_string(id)).c_str(),std::max(200,(int)v_results[1].m_max_value),v_results[1].m_min_value,v_results[1].m_max_value);
    TH1D ans2(("2"+std::to_string(id)).c_str(),("2"+std::to_string(id)).c_str(),std::max(200,(int)v_results[2].m_max_value),v_results[2].m_min_value,v_results[2].m_max_value);
    TH1D ans3(("3"+std::to_string(id)).c_str(),("3"+std::to_string(id)).c_str(),std::max(200,(int)v_results[3].m_max_value),v_results[3].m_min_value,v_results[3].m_max_value);
    TH1D ans4(("4"+std::to_string(id)).c_str(),("4"+std::to_string(id)).c_str(),std::max(200,(int)v_results[4].m_max_value),v_results[4].m_min_value,v_results[4].m_max_value);

    for(int i=0; i<num_MC; i++){
        ans0.Fill(a_vec_chis[i]);
        ans1.Fill(a_vec_pois[i]);
        ans2.Fill(a_vec_cnp[i]);
        ans3.Fill(a_vec_h0[i]);
        ans4.Fill(a_vec_h1[i]);
    }
    v_results[0].m_pdf = ans0;
    v_results[1].m_pdf = ans1;
    v_results[2].m_pdf = ans2;
    v_results[3].m_pdf = ans3;
    v_results[4].m_pdf = ans4;

    v_results[0].m_values = vec_chis;
    v_results[1].m_values = vec_pois;
    v_results[2].m_values = vec_cnp;
    v_results[3].m_values = vec_h0;
    v_results[4].m_values = vec_h1;

    is_verbose = true;

    delete[] h1_corein;
    delete[] h0_corein;
    delete[] a_specin;

    for(int i=0; i < num_bins_total_compressed; i++){
        delete[] h1_vec_matrix_inverted[i];  
        delete[] h0_vec_matrix_inverted[i];  
    }

    delete[] h1_vec_matrix_inverted;
    delete[] h0_vec_matrix_inverted;

    delete[] sampled_fullvector;
    delete[] collapsed;

    return v_results;

}




/*
   std::vector<double> SBNchi::SampleCovarianceVaryInput_getpval(SBNspec *specin, int num_MC, std::vector<double> chival){
   if(!cholosky_performed) this->PerformCholoskyDecomposition(specin); 

   int n_t = specin->full_vector.size();
   std::vector<int> nlower(chival.size(),0);

   TVectorT<double> u(n_t);
   for(int i=0; i<n_t; i++){
   u(i) = specin->full_vector.at(i);
   }

   TRandom3 * rangen = new TRandom3(0);


   TH1D ans("","",100,0,100);
   ans.GetXaxis()->SetCanExtend(kTRUE);
   is_verbose = false;
   for(int i=0; i < num_MC;i++){

   TVectorT<double> gaus_sample(n_t);
   TVectorT<double> multi_sample(n_t);
   for(int a=0; a<n_t; a++){
   gaus_sample(a) = rangen->Gaus(0,1);	
   }

   multi_sample = u + matrix_lower_triangular*gaus_sample;

   std::vector<double> sampled_fullvector(n_t,0.0);
   for(int i=0; i<n_t; i++){
   sampled_fullvector.at(i) = multi_sample(i);
   }
   SBNspec sampled_spectra(sampled_fullvector, specin->xmlname ,false);

   sampled_spectra.CollapseVector(); //this line important isnt it!

   double thischi = this->CalcChi(&sampled_spectra);
   ans.Fill(thischi);

   for(int j=0; j< chival.size(); j++){
   if(thischi>=chival.at(j)) nlower.at(j)++;
   }


   if(i%1000==0) std::cout<<"SBNchi::SampleCovarianceVaryInput(SBNspec*, int) on MC :"<<i<<"/"<<num_MC<<". Ans: "<<thischi<<std::endl;
   }
   is_verbose = true;

   std::vector<double> pval;
   for(auto n: nlower){
   pval.push_back(n/(double)num_MC);

   }

   return pval;

   }


//This one varies the input comparative spectrum, and as sucn has  only to calculate the matrix_systematics once
std::vector<double> SBNchi::SamplePoissonVaryInput_getpval(SBNspec *specin, int num_MC, std::vector<double> chival){
std::vector<int> nlower(chival.size(),0);

TRandom3 *rangen = new TRandom3(0);

TH1D ans("","",100,0,100);
//So save the core one that we will sample for
ans.GetXaxis()->SetCanExtend(kTRUE);
is_verbose = false;
for(int i=0; i < num_MC;i++){

SBNspec tmp = *specin;
tmp.ScalePoisson(rangen);
tmp.CollapseVector(); //this line important isnt it!
//tmp.PrintFullVector();

double thischi = this->CalcChi(&tmp);
ans.Fill(thischi);

for(int j=0; j< chival.size(); j++){
    if(thischi>=chival.at(j)) nlower.at(j)++;
}

if(i%1000==0) std::cout<<"SBNchi::SamplePoissonVaryInput(SBNspec*, int) on MC :"<<i<<"/"<<num_MC<<". Ans: "<<thischi<<std::endl;
}
std::vector<double> pval;
for(auto n: nlower){
    pval.push_back(n/(double)num_MC);

}

is_verbose = true;
return pval;


}
*/

//This one varies the core spectrum, and as sucn has to recalculate the matrix_systematics each stem
TH1D SBNchi::SamplePoissonVaryCore(SBNspec *specin, int num_MC){
    double center = this->CalcChi(specin);
    int nlower=0;

    TRandom3 *rangen = new TRandom3(0);

    TH1D ans("MCans","MCans",100,center-100,center+200);
    //So save the core one that we will sample for
    SBNspec core  = core_spectrum;

    is_verbose = false;
    for(int i=0; i<num_MC;i++){

        SBNspec tmp = core;
        tmp.ScalePoisson(rangen);
        this->ReloadCoreSpectrum(&tmp);
        double thischi = this->CalcChi(specin);
        ans.Fill(thischi);
        if(thischi<=center)nlower++;

        if(i%1000==0) std::cout<<"SBNchi::SamplePoissonVaryCore(SBNspec*, int) on MC :"<<i<<"/"<<num_MC<<". Ans: "<<thischi<<std::endl;
    }
    std::cout<<"pval: "<<nlower/(double)num_MC<<std::endl;

    is_verbose = true;
    delete rangen; rangen = nullptr;
    return ans;
}


int SBNchi::DrawComparisonIndividualFracMatrix(SBNspec& sig, SBNspec& data, TMatrixT<double>& frac_matrix, std::string tag){
	return this->DrawComparisonIndividualFracMatrix(sig, data, frac_matrix, tag, false);
}


int SBNchi::DrawComparisonIndividualFracMatrix(SBNspec& sig, SBNspec& data, TMatrixT<double>& frac_matrix, std::string tag, bool width_norm){
	sig.CalcFullVector();
	sig.CalcErrorVector(); 

	TMatrixT<double> collapsed_syst_matrix = this->FillSystMatrix(frac_matrix, sig.full_vector, sig.full_err_vector, true);
	return this->DrawComparisonIndividual(sig, data, collapsed_syst_matrix, tag, width_norm);
}


int SBNchi::DrawComparisonIndividual(SBNspec& sig, SBNspec& data, std::string tag){
	return this->DrawComparisonIndividual(sig, data, tag, false);
}

int SBNchi::DrawComparisonIndividual(SBNspec& sig, SBNspec& data, std::string tag,bool width_norm){

	//calculate the MC intrinsic systematic matrix for sig.
	sig.CalcErrorVector();
	TMatrixT<double> collapsed_MCIntrinsic_matrix(sig.num_bins_total_compressed, sig.num_bins_total_compressed);
	collapsed_MCIntrinsic_matrix.Zero();
	for(int i=0; i<sig.num_bins_total_compressed; i++){
		collapsed_MCIntrinsic_matrix(i,i) = pow(sig.collapsed_err_vector.at(i), 2.0);
	}

	return this->DrawComparisonIndividual(sig, data, collapsed_MCIntrinsic_matrix, tag, width_norm);
}

int SBNchi::DrawComparisonIndividual(SBNspec& sig, SBNspec& data, TMatrixT<double>& inmatrix, std::string tag){
	return this->DrawComparisonIndividual(sig, data, inmatrix, tag, false);
}


//This one compares two SBNspec, and print the chi2, pvalue for each individual distribution
//basically almost the same as SBNspec::CompareSBNspecs()
// given collapsed systematic covariance matrix, generate the comparison plots
int SBNchi::DrawComparisonIndividual(SBNspec& sig, SBNspec& data, TMatrixT<double>& collapsed_syst_matrix, std::string tag, bool width_norm){
	// default color scheme
	std::vector<int> mycol = {kAzure -9, kRed-7, kGreen-3, kBlue-6, kMagenta-3, kYellow-7,  kOrange-3, kBlue, kBlue+2,  kGreen+1,kBlue-7, kPink, kViolet, kCyan,kMagenta,kAzure};
		
	if((collapsed_syst_matrix.GetNcols()!= collapsed_syst_matrix.GetNrows()) ||  (collapsed_syst_matrix.GetNcols() != this->num_bins_total_compressed)){
		std::cerr<< "ERROR: dimension of the matrix doesn't match the total compressed bin number, or it's not a symetric matrix"<< std::endl;
		exit(EXIT_FAILURE);
	}
	


	// single photon color scheme
	bool gLEE_plot= true;
	if(gLEE_plot){
                mycol.clear();
                
                std::map<std::string, int> color_channel_map;
                std::map<std::string, std::vector<double>> rgb_channel_map={
                        {"NCDelta", {255./255.,255./255.,153./255.}},
                        {"NCDeltaLEE", {0.97,0.75,0}},
                        {"NCPi0Coh", {255./255,189./255.,189./255.}},
                        {"NCPi0NotCoh", {1,0.4,0.4}},
                        {"NCMultiPi0", {0.9,0.9,1.0}},
                        {"CC1Pi0", {0.4,0.4,1.0}},
                        {"BNBOther", {0.6,0.8,1.0}},
                        {"Nue",{0.9,0.5,0.9}},
                        {"Dirt", {0.6,0.4,0.2}},
                        {"OTPCExtra", {0.2,0.5,0.2}},
                        {"BNBext", {0.2,0.8,0.2}}
                };
                
                std::map<std::string, std::vector<double>>::iterator iter;
                std::map<std::string, int>::iterator iter_int;
                TColor* t_col = NULL;
                for(iter = rgb_channel_map.begin(); iter!= rgb_channel_map.end(); ++iter){
                        int color_index = TColor::GetFreeColorIndex();
                        t_col = new TColor(color_index, iter->second.at(0),iter->second.at(1),iter->second.at(2));
                        color_channel_map.insert({iter->first, t_col->GetNumber()});
			delete t_col; t_col = nullptr;
                }
                
                for(int is = 0; is <subchannel_names[0].size(); is++){
                        std::string isubchannel_name = subchannel_names[0][is];
                        iter_int = color_channel_map.find(isubchannel_name);
                        if(iter_int == color_channel_map.end()){
                                std::cout << "Color of " << isubchannel_name << " is not defined, choose a random color" << std::endl;
                                mycol.push_back(is);
                        }
                        else{   
				 mycol.push_back(iter_int->second);
                        }
                }
         
        }

	sig.CollapseVector();
	sig.CalcErrorVector(); 
	data.CollapseVector();

	std::vector<TH1D> temp = sig.hist;
        std::vector<TH1D> temp_comp = data.hist;
	std::vector<double> sig_collapse = sig.collapsed_vector;
	std::vector<double> data_collapse = data.collapsed_vector;

	TFile* fcompare = new TFile(("SBNfit_compare_individual_plots_"+tag+".root").c_str(), "recreate");
	fcompare->cd();
	collapsed_syst_matrix.Write("collapsed_syst_matrix");

	int matrix_index= 0;  //to keep track of what part of collapsed matrix we should extract, for a certain channel

	for(int im = 0; im <mode_names.size(); im++){
                for(int id = 0; id <detector_names.size(); id++){
                        for(int ic = 0; ic <channel_names.size(); ic++){


				bool this_run = false;
                                bool this_run_comp = false;
				std::string canvas_name = mode_names.at(im)+"_"+detector_names.at(id)+"_"+channel_names.at(ic);
				TCanvas* Cstack= new TCanvas((tag+"_"+canvas_name).c_str(),(tag+" | "+canvas_name).c_str(),1450,1200);
				Cstack->cd();
                                THStack * hs = new THStack(canvas_name.c_str(),  canvas_name.c_str());
				TLegend legStack(0.11, 0.58, 0.89, 0.89);
                                legStack.SetNColumns(2);
                                legStack.SetLineWidth(0);
                                legStack.SetLineColor(kWhite);
                                int n=0;
                                int nc=0;
                                TH1D * hcomp;
                                TH1D *hsum;
                                double hcomp_sum=0;
                                double hsum_sum=0;

				for(auto &h : temp_comp){
                                        std::string test = h.GetName();
                                        if(test.find(canvas_name)!=std::string::npos){
                                                double total_events = h.GetSumOfWeights();
                                                hcomp_sum += total_events;
						if(width_norm) h.Scale(1,"width");
						h.SetLineColor(kBlack);
						if(!this_run_comp){
                                                        hcomp = (TH1D*)h.Clone(("comp_"+canvas_name).c_str());
                                                        hcomp->Reset();
                                                }

                                                std::ostringstream out;
                                                out << std::setprecision(2) << total_events;
                                                std::string hmm = "\t";
                                                std::string tmp = h.GetName() +hmm+ out.str();


                                                hcomp->Add(&h);
                                                nc++;

                                                this_run_comp=true;

                                        }
                                }

				std::vector<double> integral_sorter;
                                std::vector<TH1*> to_sort;
                                std::vector<std::string> l_to_sort;

                                for(auto &h : temp){
                                        std::string test = h.GetName();
                                        if(test.find(canvas_name)!=std::string::npos ){

                                                double total_events = h.GetSumOfWeights();
                                                hsum_sum += total_events;
						if(width_norm) h.Scale(1,"width");
						h.GetYaxis()->SetTitle("Events/GeV");
                                                h.SetMarkerStyle(20);
                                                h.SetMarkerColor(mycol[n]);
                                                h.SetFillColor(mycol[n]);
                                                if(gLEE_plot & (test.find("BNBext")!=std::string::npos)) h.SetFillStyle(3333);
                                                h.SetLineColor(kBlack);
                                                h.SetTitle(h.GetName());
						if(!this_run){
                                                        hsum = (TH1D*)h.Clone(("sum_"+canvas_name).c_str());
                                                        hsum->Reset();
                                                }
						std::ostringstream out;
                                                out <<std::fixed<< std::setprecision(2) << total_events;
                                                std::string hmm = " | ";
						std::string tmp_name = h.GetName();
                                                std::string tmp = map_subchannel_plotnames[tmp_name] +hmm+ out.str();
						hsum->Add(&h);
						n++;

                                                this_run=true;

                                                to_sort.push_back(&h);
                                                l_to_sort.push_back(tmp);
                                                integral_sorter.push_back(total_events);

                                                if(gLEE_plot){
                                                        hs->Add(&h, "HIST");
                                                        legStack.AddEntry(&h, tmp.c_str(),"f");
                                                }

                                        }
                                }

                                if(!gLEE_plot){
					for (int i: SortIndexes(integral_sorter)) {
                                        hs->Add(to_sort.at(i), "HIST");
                                        legStack.AddEntry(to_sort.at(i), l_to_sort.at(i).c_str(),"f");
                                }
                                }

				// ************* start dealing with error bars************************
				int step = num_bins[ic];   //num of bins perchannel
				TMatrixT<double> sub_syst_matrix = collapsed_syst_matrix.GetSub(matrix_index, matrix_index+step-1, matrix_index, matrix_index+step-1);
				TMatrixT<double> sub_all_err_matrix = sub_syst_matrix;
				double mc_intrinsic_err_square=0;  //overall stats error of the MC histogram integral
				double sig_bin_err, data_bin_err;
				double bin_width;
				for(int ib=0; ib< step; ib++){
					if(width_norm) bin_width = hsum->GetBinWidth(ib+1);
					else bin_width = 1.0;

					mc_intrinsic_err_square += pow(bin_width*hsum->GetBinError(ib+1), 2.0);				
					sig_bin_err = sqrt(sub_syst_matrix(ib,ib))/bin_width;
					data_bin_err = sqrt(hcomp->GetBinContent(ib+1)/bin_width);

					hsum->SetBinError(ib+1, sig_bin_err);
					hcomp->SetBinError(ib+1, data_bin_err);

					// ******to add stats error to the matrix**********
					//Pearson chi2
					//sub_all_err_matrix(ib,ib) += bin_width*hsum->GetBinContent(ib+1);
					//Neyman chi2
					//sub_all_err_matrix(ib,ib) += bin_width*hcomp->GetBinContent(ib+1);
					//CNP chi2
					sub_all_err_matrix(ib,ib) +=  bin_width*( hcomp->GetBinContent(ib+1) >0.001 ? 3.0/(1.0/hcomp->GetBinContent(ib+1) +  2.0/hsum->GetBinContent(ib+1))  : hsum->GetBinContent(ib+1)/2.0 );

				}	
				


				// *********** start calculate chi2 and Pvalue **************
				TMatrixT<double> sub_inverted_matrix = this->InvertMatrix(sub_all_err_matrix);
				double chi=0;
				for(int ir=0; ir < step ; ir++){
				   for(int ik=0; ik< step; ik++){
					chi += sub_inverted_matrix(ir, ik)*(sig_collapse[matrix_index+ir] - data_collapse[matrix_index+ir])*(sig_collapse[matrix_index+ik] - data_collapse[matrix_index+ik]);
				   }
				}
				double pvalue = TMath::Prob(chi, step);

				otag= "SBNchi::DrawComparisonIndividual\t||";
				std::cout<< otag<<canvas_name<<" has stats error "<< sqrt(hsum_sum)<<", overall instrinsic error: "<< sqrt(mc_intrinsic_err_square) << ", and a total systematic error: "<< sqrt(sub_syst_matrix.Sum())  <<", overall total error: " << sqrt(hsum_sum + sub_syst_matrix.Sum())<< std::endl;
				std::cout<<otag<<canvas_name<< " has chi2 of " << chi << ", ndf: " << step << ", and Pvalue: "<< pvalue<< std::endl;
				matrix_index += step; //update where to start for the next iteration
				


				// ************* start the drawing process *****************
				legStack.AddEntry(hsum, Form("MC Stack | %.2f", hsum_sum), "fl");
                                legStack.AddEntry(hcomp, Form("Data | %.2f", hcomp_sum), "flp");
				if(this_run && this_run_comp){
                                        double plot_pot=5e19;

                                        double title_size_ratio=0.11;
                                        double label_size_ratio=0.11;
                                        double title_offSet_ratioY = 0.45;
                                        double title_offSet_ratioX = 1.1;

                                        double title_size_upper=0.048;
                                        double label_size_upper=0.05;
                                        double title_offSet_upper = 0.85;

                                        Cstack->cd();
					TPad *pad0top = new TPad(("pad0top_"+canvas_name).c_str(), ("pad0top_"+canvas_name).c_str(), 0, 0.40, 1, 1.0);
                                        pad0top->SetBottomMargin(0); // Upper and lower plot are joined
					pad0top->Draw();             // Draw the upper pad: pad2top
					pad0top->cd();               // pad2top becomes the current pad
					hs->Draw();
  					gStyle->SetHatchesLineWidth(1);

					//draw error bar for 'sig' SBNspec
					hsum->SetFillColor(kBlack);
                                        hsum->SetFillStyle(3354);
                                        hsum->SetMarkerSize(0);
                                        hsum->SetLineWidth(1);
                                        hsum->Draw("E2 same");
					if(width_norm) hs->GetYaxis()->SetTitle(("Events/"+channel_units.at(ic)).c_str());
					else hs->GetYaxis()->SetTitle("Events");
                                        hs->GetYaxis()->SetTitleSize(title_size_upper);
                                        hs->GetYaxis()->SetLabelSize(label_size_upper);
                                        hs->GetYaxis()->SetTitleOffset(title_offSet_upper*1.2);

					gStyle->SetEndErrorSize(5);
                                        hcomp->SetMarkerStyle(20);
                                        hcomp->SetMarkerColor(kBlack);
                                        hcomp->SetMarkerSize(1.5);
                                        hcomp->SetLineWidth(2);
                                        hcomp->SetLineColor(kBlack);
                                        hcomp->Draw("E1P same");
					hs->SetMaximum(std::max(hs->GetMaximum(), hcomp->GetMaximum())*2.2);
                                        //hs->SetMaximum(std::max(hs->GetMaximum(), hcomp->GetMaximum())*1.1);
                                        hs->SetMinimum(0.001);

                                        Cstack->Update();
                                        legStack.Draw();
                                        TText t_text(0.7*(hs->GetXaxis()->GetXmax()), 1.4*std::max(hs->GetMaximum(), hcomp->GetMaximum()), (channel_names.at(ic)).c_str());
                                        t_text.SetTextSize(0.08);
                                        t_text.Draw();

                                        Cstack->cd();
					TPad *pad0bot = new TPad(("padbot_"+canvas_name).c_str(),("padbot_"+canvas_name).c_str(), 0, 0.1, 1, 0.40);
                                        pad0bot->SetTopMargin(0);
                                        pad0bot->SetBottomMargin(0.38);
                                        pad0bot->SetGridx(); // vertical grid
					pad0bot->Draw();
                                        pad0bot->cd();       // pad0bot becomes the current pad
					
			
					TH1* ratpre = (TH1*)hcomp->Clone(("ratio_"+canvas_name).c_str());
                                        ratpre->Divide(hsum);
                                        ratpre->SetStats(false);
                                        //to draw the 1 sigma error band on the ratio plot
                                        TH1* h_err = (TH1*)hsum->Clone("error_band");
                                        h_err->Divide(hsum);
                                        for(int i=0; i<h_err->GetNbinsX(); i++){
                                                h_err->SetBinError(i+1, hsum->GetBinError(i+1)/hsum->GetBinContent(i+1));
                                                ratpre->SetBinError(i+1, hcomp->GetBinError(i+1)/hsum->GetBinContent(i+1));
                                        }


					//ratpre->Draw("hist");
                                        ratpre->Draw("E1");
					ratpre->SetFillColor(kWhite);
					ratpre->SetLineWidth(2);
                                        h_err->Draw("E2same");
                                        Cstack->Update();
                                        gStyle->SetOptStat(0);
                                        TLine *line = new TLine(ratpre->GetXaxis()->GetXmin(),1.0,ratpre->GetXaxis()->GetXmax(),1.0 );
                                        line->Draw("same");
                                        ratpre->SetLineColor(kBlack);
                                        ratpre->SetTitle("");
                                        ratpre->GetYaxis()->SetTitle("Data/Prediction");
                                        ratpre->GetXaxis()->SetTitleOffset(title_offSet_ratioX);
                                        ratpre->GetYaxis()->SetTitleOffset(title_offSet_ratioY);
                                        ratpre->SetMinimum(std::min(0.5, ratpre->GetMinimum())*0.8);
                                        //ratpre->SetMinimum(ratpre->GetMinimum()*0.97);
                                        ratpre->SetMaximum(std::max(1.5, ratpre->GetMaximum())*1.2);
                                        //ratpre->SetMaximum(ratpre->GetMaximum()*1.03);
                                        ratpre->GetYaxis()->SetNdivisions(505, kTRUE);   //change the label division in y axis
					ratpre->GetYaxis()->SetTitleSize(title_size_ratio);
                                        ratpre->GetXaxis()->SetTitleSize(title_size_ratio);
                                        ratpre->GetYaxis()->SetLabelSize(label_size_ratio);
                                        ratpre->GetXaxis()->SetLabelSize(label_size_ratio);
					ratpre->GetXaxis()->SetTitle(channel_units.at(ic).c_str());
                                        Cstack->Update();
					
					Cstack->cd();
					TLatex* lat = new TLatex();
					lat->SetTextColor(kRed-9);
					lat->SetTextSize(0.03);
					lat->DrawLatex(0.1, 0.1, Form("Data/MC: %.2f", hcomp_sum/hsum_sum));
					lat->DrawLatex(0.4, 0.1, Form("(#chi^{2}/nDOF: %.2f/%d)", chi, step));
					lat->DrawLatex(0.7, 0.1, Form("(#chi^{2} P^{val}: %.3f)", pvalue));
					Cstack->Update();
                                        Cstack->SaveAs((tag+"_chi_"+canvas_name+".pdf").c_str(),"pdf");
                                        Cstack->Write(canvas_name.c_str() );
				}

			}
                }
        }

        fcompare->Close();


	return 0;
}
