#include "abnormalgpz.h"

AbnormalGPZ::AbnormalGPZ(long double t0, long double t1, long double SamplingIncrement,
                         long double theta,long double i,long double w,long double a,long double e,long double O):
    NormalGPZ (t0, t1, SamplingIncrement,
                theta, i, w, a, e, O)
{
    std::string line;
    std::ifstream in("InitParam.txt");
    if (in.is_open())
    {
        int j=0;
        while(getline(in,line))
        {
            int index=0;
            std::string arg="";
            for(decltype (line.size()) i=0;i<=line.size();i++)
            {
                char a=line[i];
                if (a==',') a='.';
                if ((a!=' ')||(a!='\r')||(a!='\0')) arg+=a;
                if ((a==' ') || (a=='\r')||(a=='\0'))
                {
                        if (index==0) this->arg[j].e=std::atof(arg.c_str())*1e-10L;
                        if (index==1) this->arg[j].x=std::atof(arg.c_str())*1000L;
                        if (index==2) this->arg[j].y=std::atof(arg.c_str())*1000L;
                        if (index==3) this->arg[j].z=std::atof(arg.c_str())*1000L;
                        arg=" ";
                        index++;
                }
            }
            j++;
        }
    }
    in.close();
}

TVector  AbnormalGPZ::getRight( const TVector& X, long double t)
{
    // L,Fi, ro
    // x,y,z
    TVector X_Coordin=X.Concat(0,2);
    TVector dX=NormalGPZ::getRight(X,t);

    TVector deltaGGeo(3);

    long double sum[3]={0,0,0};

    for(int i=0;i<60;i++)
    {
        long double ri=sqrtl(powl(X_Coordin[0]-arg[i].x,2)+powl(X_Coordin[1]-arg[i].y,2)+powl(X_Coordin[2]-arg[i].z,2));
        sum[0]+=arg[i].e*(X_Coordin[0]-arg[i].x)/powl(ri,3);
        sum[1]+=arg[i].e*(X_Coordin[1]-arg[i].y)/powl(ri,3);
        sum[2]+=arg[i].e*(X_Coordin[2]-arg[i].z)/powl(ri,3);
    }

    for(int i=0;i<3;i++)
    deltaGGeo[i]=-mu_s*sum[i];

    TVector Y(6);
    for (int i=0; i!=3; ++i)
        Y.SetItem(i, dX[i]);
    for (int i=3; i!=6; ++i)
        Y.SetItem(i, dX[i]+deltaGGeo[i-3]);
    return Y;
}
