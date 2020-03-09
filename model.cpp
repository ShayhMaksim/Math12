//---------------------------------------------------------------------------

#include "model.h"



//---------------------------------------------------------------------------
void TModel::addResult( const TVector& X, long double& t )
{

	if (N == Result.GetRowCount())
		Result.Resize(N + 1, getOrder() + 1);

	Result(N, 0) = t;
	for (int i = X.GetSize(); i>0; i--)
		Result(N, i) = X[i-1];
	N++;
}

void TModel::clearResult() 
{ 

	Result.Resize(0, 0);
	N = 0; 
}

void TModel::prepareResult() 
{ 

	Result.Resize((int)((t1 - t0)/SamplingIncrement) + 1, getOrder() + 1);

	N = 0; 
}



//---------------------------------------------------------------------------


//---------------------------------------------------------------------------

TDormandPrince::TDormandPrince():TIntegrator()
{
	long double v=1;
        		while ((1+v)>1)
        		{
        			u=v;
        			v=v/2;
        		}
}



void TDormandPrince::Run(TModel* Model)
{
    long double
				t = Model->getT0(),

				t_out = t,

				t1 = Model->getT1(),

				h,

				h_new = Model->getSamplingIncrement(),

			    e = 0;

    TVector
            X = Model->getInitialConditions(),

            X1( X.GetSize() ),

            X2( X.GetSize() ),

            Xout ( X.GetSize() );

	Model->prepareResult();

	int N = 0;

	while ( t < t1 )
    {

        h = h_new;

        K[0] = Model->getRight( X, t );
		K[1] = Model->getRight(X+K[0]*h*a[1][0] , t+c[1]*h);
		K[2] = Model->getRight(X+K[0]*h*a[2][0]+K[1]*a[2][1]*h, t+c[2]*h);
		K[3] = Model->getRight(X+K[0]*h*a[3][0]+K[1]*a[3][1]*h+K[2]*a[3][2]*h, t+c[3]*h);
		K[4] = Model->getRight(X+K[0]*h*a[4][0]+K[1]*a[4][1]*h+K[2]*h*a[4][2]+K[3]*h*a[4][3], t+c[4]*h);
		K[5] = Model->getRight(X+K[0]*h*a[5][0]+K[1]*a[5][1]*h+K[2]*h*a[5][2]+K[3]*h*a[5][3]+K[4]*h*a[5][4], t+c[5]*h);
		K[6] = Model->getRight(X+K[0]*h*a[6][0]+K[1]*a[6][1]*h+K[2]*h*a[6][2]+K[3]*h*a[6][3]+K[4]*h*a[6][4]+K[5]*h*a[6][5], t+c[6]*h);


        X1 = X+K[0]*b1[0]*h+K[1]*b1[1]*h+K[2]*b1[2]*h+K[3]*b1[3]*h+K[4]*b1[4]*h+K[5]*b1[5]*h+K[6]*b1[6]*h;
        X2 = X+K[0]*b2[0]*h+K[1]*b2[1]*h+K[2]*b2[2]*h+K[3]*b2[3]*h+K[4]*b2[4]*h+K[5]*b2[5]*h+K[6]*b2[6]*h;


        long double max[X.GetSize()];
        double e1=0;
        for(int i=0;i<X.GetSize();i++)
        {
            long double max1[4]{(0.00001) , abs(X1[i]), abs(X[i]) , 0.5*u/Eps};
        	long double a=max1[0];

        	for(int j=0;j<4;j++)
        			if (a<=max1[j]) a=max1[j];

        	max[i]=a;
        	e1=e1+pow(h*(X1.GetItem(i)-X2.GetItem(i))/max[i],2);
        }
        e=pow(e1/X.GetSize(), 0.5);


        long double min;
        long double max2;


        if (5>pow((e/Eps), 0.2)/0.9) { min=pow(e/Eps, 0.2)/0.9; } else min=5;

        if (0.1>=min) {max2=0.1 ;} else max2=min;

        h_new = h/max2;


        if ( e > Eps )
            continue;


		while ( (t_out < t + h) && (t_out <= t1) )
        {
            long double theta = (t_out - t)/h,
                        b[6];

            b[0] =  theta*(1+theta*(-1337./480+theta*(1039./360+theta*(-1163./1152))));
            b[1] = 0;
            b[2] = 100.*theta*theta*(1054./9275+theta*(-4682./27825+theta*(379./5565)))/3.;
            b[3] = -5.*theta*theta*(27./40+theta*(-9./5+theta*(83./96)))/2.;
            b[4] = 18225.*theta*theta*(-3./250+theta*(22./375+theta*(-37./600)))/848.;
            b[5] = -22.*theta*theta*(-3./10+theta*(29./30+theta*(-17./24)))/7.;

        	Xout=X;
        		for(int j=0;j<6;j++)
        		{
        			Xout = Xout+K[j]*(b[j]*h);
        		}

            Model->addResult( Xout, t_out );

            t_out += Model->getSamplingIncrement();
        }

        X = X2;
        t += h;

        N++;
    }
}



