//---------------------------------------------------------------------------

#ifndef modelH
#define modelH


#include "MyLinearAlgebra.h"
#include <math.h>

//---------------------------------------------------------------------------

using namespace MyLinearAlgebra;

class TModel
{
    protected:

	    TVector X0;

        long double SamplingIncrement;

        long double t0, t1;

        TMatrix Result;

		int N;
		
    public:

        TModel()
            : SamplingIncrement( 1e-1 )
            , t0( 0 )
            , t1( 1 )
			, N( 0 )
        {}

        virtual TVector getRight( const TVector& X, long double t ) = 0;

        inline TVector getInitialConditions() const { return X0; }

        inline int getOrder() const { return X0.GetSize(); }

        inline long double getSamplingIncrement() const { return SamplingIncrement; }

        inline long double getT0() const { return t0; }
        inline long double getT1() const { return t1; }
       

		inline TMatrix getResult() const { return Result; }

        virtual void addResult( const TVector& X, long double &t );

		virtual void clearResult();

		virtual void prepareResult();
		virtual ~TModel(){};
};

//---------------------------------------------------------------------------
class TIntegrator
{
    protected:

        long double Eps;
    public:

        TIntegrator() : Eps( 1e-8 ) {}
        inline void setPrecision( long double Eps ) { this->Eps = Eps; }
        inline long double getPrecision() const { return Eps; }

        virtual void Run(TModel* Model)=0;
        virtual ~TIntegrator(){}
};

//---------------------------------------------------------------------------

class TDormandPrince : public TIntegrator
{
    private:

    const long double c[7] = { 0, 1./5, 3./10, 4./5, 8./9, 1., 1. };
	const long double a[7][6] = {
	    { 0. },
	    { 1./5 },
	    { 3./40, 9./40 },
	    { 44./45, -56./15, 32./9 },
	    { 19372./6561, -25360./2187, 64448./6561, -212./729 },
	    { 9017./3168, -355./33, 46732./5247, 49./176, -5103./18656 },
	    { 35./384, 0., 500./1113, 125./192, -2187./6784, 11./84 }
	};
	const long double b1[7] = { 35./384, 0., 500./1113, 125./192, -2187./6784, 11./84, 0 };
	const long double b2[7] = { 5179./57600, 0., 7571./16695, 393./640, -92097./339200, 187./2100, 1./40 };


        TVector K[7];

		long double u;
    public:
        TDormandPrince();

        virtual void Run(TModel* Model);

};


#endif
 
