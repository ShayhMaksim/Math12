#include "modelems.h"

ModelEMS::ModelEMS(long double t0, long double t1, long double SamplingIncrement):TModel()
{
    this->t0=t0;
    this->t1=t1;
    this->SamplingIncrement=SamplingIncrement;

    for(int i=0;i<3;i++){
        mu[i]=G*m[i];
    }
       for(int i=0;i<3;i++){
            kappa[i]=mu[i]/mu[2];
       }

    T=2*M_PI*a*sqrt(a/mu[2]);

    //Расчет начального положения
    TVector xi10(new long double[3]{5.77103476e-01,-8.32119380e-01,  -4.85579076e-05},3);
    TVector xi20(new long double[3]{5.75566367e-01,-8.29881892e-01,  -5.36699450e-05},3);
    TVector xi30=-kappa[0]*xi10-kappa[1]*xi20;

    TVector VE(new long double[3]{1.388633512282171E-02,9.678934168415631E-03,3.429889230737491E-07},3);
    TVector VM(new long double[3]{1.434571674368357E-02,9.997686898668805E-03,-5.149408819470315E-05},3);

    u=sqrt(mu[2]/a)/2./M_PI;

    TVector VE0(3);
    TVector VM0(3);

    for(int i=0;i<3;i++){
        VE0[i]=VE[i] * a / (local_T*u);
        VM0[i]=VM[i] * a / (local_T*u);
    }

    TVector VS0= - kappa[0] * VM0 - kappa[1] * VE0;



    X0=TVector(new long double [18]{
            xi10[0],xi10[1], xi10[2],
            xi20[0],xi20[1], xi20[2],
            xi30[0],xi30[1], xi30[2],
            VM0[0], VM0[1],   VM0[2],
            VE0[0], VE0[1],   VE0[2],
            VS0[0], VS0[1],   VS0[2],
},18);



}

TVector ModelEMS::getRight( const TVector& X, long double t)
{
    TVector r12=X.Concat(3,5)-X.Concat(0,2);
    TVector r13=X.Concat(6,8)-X.Concat(0,2);
    TVector r23=X.Concat(6,8)-X.Concat(3,5);

//    TVector X1(3),X2(3),X3(3);

//    X1[0]=X[0];
//    X1[1]=X[1];
//    X1[2]=X[2];

//    X2[0]=X[3];
//    X2[1]=X[4];
//    X2[2]=X[5];

//    X3[0]=X[6];
//    X3[1]=X[7];
//    X3[2]=X[8];

//    TVector r12=X2-X1;
//    TVector r13=X3-X1;
//    TVector r23=X3-X2;

    TVector Y(18);
    for(int i=0;i<9;i++){
        Y[i]=X[i+9];
    }

    TVector E0,E1,E2;
    E0=k*kappa[1]/pow(r12.length(),3)*r12+k*kappa[2]/pow(r13.length(),3)*r13;
    E1=-k*kappa[0]/pow(r12.length(),3)*r12+k*kappa[2]/pow(r23.length(),3)*r23;
    E2=-k*kappa[0]/pow(r13.length(),3)*r13-k*kappa[1]/pow(r23.length(),3)*r23;

    Y[9]=E0[0];
    Y[10]=E0[1];
    Y[11]=E0[2];

    Y[12]=E1[0];
    Y[13]=E1[1];
    Y[14]=E1[2];

    Y[15]=E2[0];
    Y[16]=E2[1];
    Y[17]=E2[2];

    return Y;
}
