#include "mainwindow.h"
#include <QApplication>
#include "adecorator.h"
int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    MainWindow w;
    w.show();
//    EMSDecorator *ar=new EMSDecoratorNGPZ(0,0.01,0.00001,0,M_PI_4,0,10000000,0,0);
//    TIntegrator * integrator=new TDormandPrince();
//    integrator->setPrecision(1e-18);
//    ar->Run();
//    integrator->Run(ar);
//    delete integrator;
//    ar->getResult().toText("go");
    return a.exec();
}
