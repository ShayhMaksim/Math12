#include "spaceformule.h"

namespace FSpace {

double JulianDate(int year, int month, int day, int hour, int minute, int second) {
    int a = (14 - month) / 12;
    int Y = year + 4800 - a;
    int M = month + 12 * a - 3;
    int JDN = day + (153 * M + 2) / 5 + 365 * Y+ Y / 4 - Y / 100 + Y / 400 -32045;
    double JD = JDN + (hour - 12)/24.0 + minute/1440.0 + second/86400.0;
    return JD;
}

double SiderialTime(int year, int month, int day, int hour, int minute, int second){
    int d = (int)(JulianDate(year, month, day, hour, minute, second) - 2451544.5);
    double t = d / 36525.0;
    double S_G_0 = 24110.584841 + 8640184.812866 * t + 0.093104 * pow(t, 2.) - 6.2e-6 * pow(t, 3.);
    return S_G_0;
}

TVector intoGSK(TVector R)
{
     long double RE=6371.3;
     TVector X(3);
     X[0]=(RE)*cos(R[1])*cos(R[0]);
     X[1]=(RE)*cos(R[1])*sin(R[0]);
     X[2]=(RE)*sin(R[1]);
     return X;
}

TMatrix intoA(double Fi,double S)
{
    TMatrix A(3,3);
    A(0,0)=-sin(Fi)*cos(S);
    A(0,1)=-sin(Fi)*sin(S);
    A(0,2)=cos(Fi);
    A(1,0)=cos(Fi)*cos(S);
    A(1,1)=cos(Fi)*sin(S);
    A(1,2)=-sin(Fi);
    A(2,0)=-sin(S);
    A(2,1)=cos(S);
    A(2,2)=0;
    return A;
}

TMatrix intoA2(long double theta,long double i,long double u)
{
    TMatrix A(3,3);
    A(0,0)=cos(u)*cos(theta)-sin(u)*sin(theta)*cos(i);
    A(0,1)=-sin(u)*cos(theta)-cos(u)*sin(theta)*cos(i);
    A(0,2)=sin(i)*sin(theta);
    A(1,0)=cos(u)*sin(theta)+sin(u)*cos(theta)*cos(i);
    A(1,1)=-sin(u)*sin(theta)+cos(u)*cos(theta)*cos(i);
    A(1,2)=-sin(i)*cos(theta);
    A(2,0)=sin(u)*sin(i);
    A(2,1)=cos(u)*sin(i);
    A(2,2)=cos(i);
    return A;
}

TVector intoGeo(const TVector & R)
{
   TVector G(3);
   G[0]=R[2]*cos(R[1])*cos(R[0]);
   G[1]=R[2]*cos(R[1])*sin(R[0]);
   G[2]=R[2]*sin(R[1]);
   return G;
}

TVector intoSphere(const TVector & R)
{
   TVector G(3);
   G[0]=atan2(R[1],R[0]);//y/x
   G[1]=atan2(R[2],sqrt(R[0]*R[0]+R[1]*R[1]));// z/(x^2+y^2)
   G[2]=sqrt(R[0]*R[0]+R[1]*R[1]+R[2]*R[2]);// ro=x^2+y^2+z^2
   return G;
}

TVector KeplerElements(TVector R, TVector V,double mu)
{
    TVector c = R^V;
    TVector f = (V^c) - (R.norm()*mu);
    double p = c.length()*c.length()/mu;
    double e = f.length()/mu;
    double a = p/(1-e*e);
    TVector e_x(3);
    e_x.SetItem(0,1);
    TVector e_z(3);
    e_z.SetItem(2,1);
    TVector e_c = c.norm();
    TVector e_r = R.norm();
    TVector e_f = f.norm();
    TVector e_omega = (e_z^e_c).norm();
    double i = acos(e_z*e_c);
    double TVR = atan(((e_x^e_omega)*e_z)*(1/(e_x*e_omega)));
    double omega = atan(-((e_omega^e_f)*e_c)*(1/(e_omega*e_f)));
    double V_kep = atan(((e_f^e_r)*e_c)*(1/(e_f*e_r)));
    TVector Kepler(6);
    Kepler.SetItem(0, i); // i
    Kepler.SetItem(1, TVR); // TVR
    Kepler.SetItem(2, omega); // omega
    Kepler.SetItem(3, a); // a
    Kepler.SetItem(4, e); // e
    Kepler.SetItem(5, V_kep); // V
    return Kepler;
}

TMatrix Mat_SK_iz_svyaz(double U, double Omega, double i)
{
    TMatrix A(3,3);
    A(0,0) = cos(U) * cos(Omega) - sin(U) * sin(Omega) * cos(i);
    A(0,1) = cos(U) * sin(Omega) + sin(U) * cos(Omega) * cos(i);
    A(0,2) = sin(U) * sin(i);

    A(1,0) = -sin(U) * cos(Omega) - cos(U) * sin(Omega) * cos(i);
    A(1,1) = -sin(U) * sin(Omega) + cos(U) * cos(Omega) * cos(i);
    A(1,2) = cos(U) * sin(i);

    A(2,0) = sin(Omega) * sin(i);
    A(2,1) = -sin(i) * cos(Omega);
    A(2,2) = cos(i);

    return A;
}

}
