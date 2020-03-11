#include "adecorator.h"

ADecorator::ADecorator()
{

}

EMSDecorator::EMSDecorator(long double t0, long double t1, long double SamplingIncrement,
                           long double theta,long double i,long double w,long double a,long double e,long double O, std::string type)
{
    EMS=new ModelEMS(t0,t1,SamplingIncrement);

    if (type=="Sattelite")
    {
        Model=new Sattelite(t0,t1*EMS->T,SamplingIncrement*EMS->T,theta,i,w,a,e,O);
    }

    if(type=="NormalGPZ")
    {
        Model=new NormalGPZ(t0,t1*EMS->T,SamplingIncrement*EMS->T,theta,i,w,a,e,O);
    }
    this->t0=t0;
    this->t1=t1*EMS->T;
    this->SamplingIncrement=SamplingIncrement*EMS->T;

    this->GlobalT=t0;
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

void EMSDecorator::Run()
{
    TIntegrator * integrator=new TDormandPrince();
    integrator->setPrecision(1e-18);
    integrator->Run(EMS);
    delete integrator;

    for(int i=0;i<EMS->getResult().GetRowCount();i++)
    {
        TVector arg(9);
        for(int j=0;j<9;j++)
            arg[j]=EMS->getResult()(i,j+1);
        vec_dX_EMS.push_back(std::move(arg));
    }

}

TVector EMSDecorator::getRight( const TVector& X, long double t )
{
    TVector Y(6);
    TVector dX_EMS(9);
   if (t<GlobalT+SamplingIncrement){
       dX_EMS=vec_dX_EMS.at(GlobalIndex);
   } else {
       if(EMS->getResult().GetRowHigh()>GlobalIndex)
       {
           GlobalT+=SamplingIncrement;
           GlobalIndex++;
       }
   }
   for(int i=0;i<3;i++)
       Y[i]=X[i+3];


   TVector dx_Sattelite=Model->getRight(X,t);
   TVector RMoon(3);
   TVector RSun(3);


   RMoon=(dX_EMS.Concat(0,2)-dX_EMS.Concat(3,5))*EMS->a;
   for(int i=0;i<3;i++)
       RMoon[i]-=dx_Sattelite[i];

   RSun=(dX_EMS.Concat(6,8)-dX_EMS.Concat(3,5))*EMS->a;
   for(int i=0;i<3;i++)
       RSun[i]-=dx_Sattelite[i];

   Y[3]=dx_Sattelite[3]+EMS->mu[0]*EMS->G/pow(RMoon.length(),3)*RMoon[0]+EMS->mu[2]*EMS->G/pow(RSun.length(),3)*RSun[0];
   Y[4]=dx_Sattelite[4]+EMS->mu[0]*EMS->G/pow(RMoon.length(),3)*RMoon[1]+EMS->mu[2]*EMS->G/pow(RSun.length(),3)*RSun[1];
   Y[5]=dx_Sattelite[5]+EMS->mu[0]*EMS->G/pow(RMoon.length(),3)*RMoon[2]+EMS->mu[2]*EMS->G/pow(RSun.length(),3)*RSun[2];

   return Y;
}




