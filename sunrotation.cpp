#include "sunrotation.h"


SunRotation::SunRotation(long double t0, long double t1, long double SamplingIncrement): TModel()
{
    this->t0=t0;
    this->t1=t1;
    this->SamplingIncrement=SamplingIncrement;
}

TVector SunRotation::getRight(const TVector& X, long double t)
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


void SunRotation::setInitialPosition(long double theta, long double i, long double w, long double a, long double e, long double O)
{
    long double p=a*(1-e*e);
    TVector r_(3);
    r_[0]=p/(1+e*cos(O));
    r_[1]=0;
    r_[2]=0;
    long double u=w+O;
    TMatrix A=intoA2(theta,i,u);
    TVector r=A*r_;
    long double mu=398600.436e+9L;
    TVector V_(3);
    V_[0]=sqrt(mu/p)*e*sin(O);
    V_[1]=sqrt(mu/p)*(1+e*cos(O));
    V_[2]=0;
    TVector V=A*V_;
    X0=TVector(6);
    for(int i=0;i<3;i++)
    {
        X0[i]=r[i];
    }
    for(int i=3;i<6;i++)
    {
        X0[i]=V[i-3];
    }
}

