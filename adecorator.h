#ifndef ADECORATOR_H
#define ADECORATOR_H
#include "model.h"
#include "modelems.h"
#include "Sattelite.h"
#include "normalgpz.h"
#include <thread>
class ADecorator:public TModel
{
public:
   TModel * Model;
    ADecorator();
    virtual void Run()=0;
};

class EMSDecorator:public ADecorator
{
public:
    EMSDecorator(long double t0, long double t1, long double SamplingIncrement,
                 long double theta,long double i,long double w,long double a,long double e,long double O,std::string type);
    virtual TVector getRight( const TVector& X, long double t ) override;
    void Run() override;
    ModelEMS * EMS;
    inline void setEMS( ModelEMS * EMS){
        this->EMS=EMS;
    }
    long double GlobalT;
    int GlobalIndex=0;
    std::vector<TVector> vec_dX_EMS;
};



#endif // ADECORATOR_H
