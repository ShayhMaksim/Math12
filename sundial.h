#ifndef SUNDIAL_H
#define SUNDIAL_H

#include "model.h"

class Sundial : public TModel
{
    const double mu_s = 132712.43994e+6;
public:
    Sundial(long double t0, long double t1, long double SamplingIncrement,const TVector& X);
    TVector getRight( const TVector& X, long double t) override;
};

#endif // SUNDIAL_H
