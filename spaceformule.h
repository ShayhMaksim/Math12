#ifndef SPACEFORMULE_H
#define SPACEFORMULE_H
#include <cmath>
#include "MyLinearAlgebra.h"
using namespace MyLinearAlgebra;
namespace FSpace{
double JulianDate(int year, int month, int day, int hour, int minute, int second);
double SiderialTime(int year, int month, int day, int hour, int minute, int second);
TVector intoGSK(TVector R);
TMatrix intoA(double Fi,double S);
TMatrix intoA2(long double theta,long double i,long double u);
TVector intoSphere(const TVector & R);
TVector intoGeo(const TVector & R);
TVector KeplerElements(TVector R, TVector V,double mu);
TMatrix Mat_SK_iz_svyaz(double U, double Omega, double i) ;
}
#endif // SPACEFORMULE_H
