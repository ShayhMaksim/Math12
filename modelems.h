#ifndef MODELEMS_H
#define MODELEMS_H
#include "model.h"

class ModelEMS : public TModel
{
  public:
    long double G=6.67e-11;//гравитационная постоянная
    long double m[3]={7.349e22,5.792e24,1.989e30};
    long double mu[3];
    long double kappa[3]; //нормировка графитац параметров
    long double T;// безразмерное время
    const double a=1.495978707e11;
    const long double k=4.L * M_PI * M_PI;
    long double u;
    long double local_T=86400;
public:
    ModelEMS(long double t0, long double t1, long double SamplingIncrement);
    TVector getRight( const TVector& X, long double t) override;
};

#endif // MODELEMS_H
