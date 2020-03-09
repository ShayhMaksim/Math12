#ifndef SUNROTATION_H
#define SUNROTATION_H
#include "model.h"
#include "spaceformule.h"
#include <QVector>
using namespace FSpace;
class Sattelite : public TModel
{
protected:
     const double mu_s = 398600.436e+9;
public:
    Sattelite(long double t0, long double t1, long double SamplingIncrement);
    TVector getRight( const TVector& X, long double t) override;
    void setInitialPosition(long double theta,long double i,long double w,long double a,long double e,long double O);
};



#endif // SUNROTATION_H
