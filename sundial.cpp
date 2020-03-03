#include "sundial.h"

Sundial::Sundial(long double t0, long double t1, long double SamplingIncrement,const TVector& X)
{
    this->t0=t0;
    this->t1=t1;
    this->SamplingIncrement=SamplingIncrement;
    X0=X;
}

TVector Sundial::getRight( const TVector& X, long double t)
{
    TVector Y(6);
    for (int i=0; i!=3; ++i)
        Y.SetItem(i, X[i+3]);
    TVector X_Coordin = X.Concat(0,2);

    TVector X_Dif = -mu_s*X_Coordin*(1/(pow(X_Coordin.length(),3)));
    for (int i=3; i!=6; ++i)
        Y.SetItem(i, X_Dif[i-3]);
   return Y;
}
