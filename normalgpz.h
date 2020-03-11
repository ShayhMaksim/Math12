#ifndef NORMALGPZ_H
#define NORMALGPZ_H
#include "Sattelite.h"
#include <QVector>
class NormalGPZ : public Sattelite
{
protected:

    const long double omega = 7.292115E-5L;
    const long double a = 6378136;

public:
    QVector<double> arguments[3];
    NormalGPZ(long double t0, long double t1, long double SamplingIncrement,
              long double theta,long double i,long double w,long double a,long double e,long double O);
    TVector getRight( const TVector& X, long double t);
    long double J0(int n);
    long double  C0N0(int n);
    long double  delta(int k);
    long double  Pnm( int  n,int  m,long double Fi);
    long double Pn_( int n,long double Fi);
    void toQVector();
};

#endif // NORMALGPZ_H
