#ifndef PROB_H_
#define PROB_H_

#include <cmath>
#include <vector>
#include <complex>
#include <iostream>
#include <TRandom3.h>
#include <TMatrixT.h>
#include <sstream>
#include <iomanip>

#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>


namespace sbn{

/*************************************************************
 *************************************************************
 *	ToDO:
 *	Actually:  add tau to NeutrinoModel. Simple and quick, but have to adjust constructors
 ************************************************************
 ************************************************************/

double SmearEnergyGaussian(double E, double percen, TRandom3 * rangen);


struct ComplexMatrix{
	ComplexMatrix(int dim);
	
	int dimension;
	TMatrixT<double> real;
	TMatrixT<double> imag;
	int mult(ComplexMatrix* in);
	int multI();
	int mult(double );
	int mult(double ,double);
	int add(ComplexMatrix* in );
	int conj();



	std::vector<double> MatrixExponent();
	std::vector<double> MatrixExponentTest(double L, std::vector<double> *, ComplexMatrix *);


	std::vector<double> GetEigenStuff(std::vector<double> *, ComplexMatrix *);

	int SetRotation(int,int, double);
	int SetComplexRotation(int,int,double,double);
	int SetDiagonal(std::vector<double> ms);
	int SetIdentity();

	int HermitianConjugate();
	int transpose();

	int Print();
};



struct NeutrinoModel{
	double mNu[3], Ue[3], Um[3], phi[3];
	double dm41Sq, dm51Sq, dm61Sq, dm54Sq, dm64Sq, dm65Sq;
	std::vector< std::vector < std::complex<double> > >  U; 
	//std::vector<std::vector<double>> dm;

	int numsterile;
	std::string mass_tag;

	//constructors!! Should overload these immensely for  3+1, 3+2, 3+3 and NULL
	NeutrinoModel();
	NeutrinoModel(double * mn, double * ue, double *um, double *ph);
	NeutrinoModel(double m4, double ue4, double um4);
	NeutrinoModel(double m4, double ue4, double um4, bool setMassTag);

	//
	NeutrinoModel(std::vector<double> mass_splittings, std::vector<double>angles, std::vector<double> phases);

	double UUem;
	double UUme;
	double UUmm;
	double UUee;

 	void Printall();

	void zero();
	void difference();

	//Might as well have the oscProb in the neutrioModel class as all the things are aready here.
	double oscProb(int init, int fin, double Ev, double L);
	double oscProb_dis(int a, double  Ev, double L);
	double oscProb_app(int a,int b, double  Ev, double L);

	double oscProbSin(double Ev, double L);
	double oscProbSinSq(double Ev, double L);


	double oscAmp(int a, int b, int which_dm, int sqornot);
	double oscAmp_dis(int a, int which_dm);
	double oscAmp_app(int a, int b, int which_dm, int sqornot);
};


// Simple 3+1 testing probabilities, completely redundant now.
double Pmue(double L, double E, double dm, double sin2);
double Pmm(double L, double Ev, double Dm, double sinSq2thmm);


//arbitraty U matrix version overload?
double oscProb(int a, int b, double Ev, double L,  std::vector< std::vector < std::complex<double> > >    U, std::vector<std::vector<double> > dm);



};
#endif
