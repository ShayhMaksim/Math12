#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QLogValueAxis>
#include <QLineSeries>
#include <QValueAxis>
#include <QChart>
#include <QChartView>
#include <thread>
#include "sundial.h"
#include "Sattelite.h"
#include "modelems.h"
#include "normalgpz.h"
#include "adecorator.h"
QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE


using namespace QtCharts;
class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private slots:
    void on_pushButton_2_clicked();

    void on_MainWindow_iconSizeChanged(const QSize &iconSize);

private:
    double theta=0, i=M_PI_4, w=0, a=10000000,  e=0, O=0;
    Ui::MainWindow *ui;
    ModelEMS * E_S;
    NormalGPZ * E_M;
    ADecorator *Decorator;
    ADecorator *DecoratorNGPZ;
    void setSeries(QLineSeries *series,TMatrix Result,int first,int second,int j,int k);
    void setDialog(long double theta,long double i,long double w,long double a,long double e,long double O);
};

#endif // MAINWINDOW_H
