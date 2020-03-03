#ifndef ABNORMALGZ_H
#define ABNORMALGZ_H
#include "normalgpz.h"

struct Element
{
    long double e;
    long double x;
    long double y;
    long double z;
};

class AbnormalGPZ : public NormalGPZ
{
    Element arg[60];
public:
    AbnormalGPZ(long double t0, long double t1, long double SamplingIncrement0);
    TVector getRight( const TVector& X, long double t);
};

#endif // ABNORMALGZ_H
