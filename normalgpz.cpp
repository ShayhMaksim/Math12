#include "normalgpz.h"

NormalGPZ::NormalGPZ(long double t0, long double t1, long double SamplingIncrement,
                     long double theta,long double i,long double w,long double a,long double e,long double O):
    Sattelite (t0,t1,SamplingIncrement,theta,i,w,a,e,O)
{

}

TVector  NormalGPZ::getRight( const TVector& X, long double t)
{
    // L,Fi, ro
    // x,y,z
    TVector X_Coordin = X.Concat(0,2);
    TVector Sph=intoSphere(X_Coordin);

    TVector GSph(3);
    long double Circle1=0;
    long double Circle2=0;
    for(int n=2;n<=8;n=n+2)
    {
        Circle1+=(n+1)*powl(a/Sph[2],n)*C0N0(n)*Pnm(n,0,Sph[1]);
        Circle2+=powl(a/Sph[2],n)*C0N0(n)*Pn_(n,Sph[1]);
    }

    GSph[0]=-mu_s/(Sph[2]*Sph[2])-mu_s/(Sph[2]*Sph[2])*Circle1;//+omega*omega*Sph[2]*cos(Sph[1])*cos(Sph[1]);//ro
    GSph[1]=mu_s/(Sph[2]*Sph[2])*Circle2;//-omega*omega*Sph[2]*Sph[2]/2*sin(2*Sph[1]);//Fi
    GSph[2]=0;//L

    TVector GGeo(3);
    TMatrix A(3,3);
    long double r0=sqrtl(X_Coordin[0]*X_Coordin[0]+X_Coordin[1]*X_Coordin[1]);

    A(0,0)=X_Coordin[0]/Sph[2];
    A(0,1)=-X_Coordin[0]*X_Coordin[2]/(Sph[2]*r0);
    A(0,2)=-X_Coordin[1]/r0;
    A(1,0)=X_Coordin[1]/Sph[2];
    A(1,1)=-X_Coordin[1]*X_Coordin[2]/(Sph[2]*r0);
    A(1,2)=X_Coordin[0]/r0;
    A(2,0)=X_Coordin[2]/Sph[2];
    A(2,1)=r0/Sph[2];
    A(2,2)=0;

    GGeo=A*GSph;


    TVector Y(6);
    for (int i=0; i!=3; ++i)
        Y.SetItem(i, X[i+3]);
    for (int i=3; i!=6; ++i)
        Y.SetItem(i, GGeo[i-3]);

    return Y;

}



 long double NormalGPZ::J0(int n)
 {
     long double J0=0;
     if (n==2) J0=1082.62575E-6L;
     if (n==4) J0=-2.37089E-6L;
     if (n==6) J0=6.08E-9L;
     if (n==8) J0=-1.4E-11L;
     return J0;
 }

 long double NormalGPZ::C0N0(int n)
 {
     long double C0N0=-J0(n)/sqrtl(2*n+1);
     return C0N0;
 }

 long double NormalGPZ::delta(int k)
 {
     if (k == 0) { return 0.5L; }
              else return 1;
 }

 long double NormalGPZ::Pnm(int n,int m, long double Fi)
 {
     long double Result=0;
     if ((n==m) && (n!=0) && (m!=0)) {
         Result=Pnm(n-1,m-1,Fi)*cos(Fi)*sqrt((2*n+1)/(2*n*delta(m-1)));
     }
     if (n>m)
     {
         Result=Pnm(n-1,m,Fi)*sin(Fi)*sqrtl((4*n*n-1)/(n*n-m*m))-Pnm(n-2,m,Fi)*sqrtl(((n-1)*(n-1)-m*m)*(2*n+1)/((n*n-m*m)*(2*n-3)));
     }
     if (n<m) Result=0;
     if ((n==m) && (m==0)) Result=1;
     return Result;
 }

 long double NormalGPZ::Pn_(int n,long double Fi)
 {
     long double Result = sqrtl(0.5L*n*(n+1))*Pnm(n,1,Fi);
     return Result;
 }

 void NormalGPZ::toQVector()
 {
     for (int i=0;i<Result.GetRowCount();i++)
             {
                     for(int j=1;j<Result.GetColCount()-3;j++)
                     {
                         arguments[j-1].push_back(Result(i,j));
                     }
             }
 }
