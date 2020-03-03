#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QLineSeries>
MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);


    E_S=new ModelEMS(0,1800,0.018);

    E_M=new SunRotation(0,50000,1);

    E_M->getInitialPosition(0,0,0,10000000,0,0);
    ui->progressBar->setRange(0,100);
    ui->progressBar->setValue((float)E_M->getT0()/E_M->getT1()*100);
    ui->progressBar_2->setRange(0,100);
    ui->progressBar_2->setValue((float)E_S->getT0()/E_S->getT1()*100);


    std::thread thread_1([&](){
        TIntegrator * integrator=new TDormandPrince();
        integrator->setPrecision(1e-16);
        integrator->Run(E_M);
        delete integrator;
    });

    std::thread thread_2([&](){
        TIntegrator * integrator2=new TDormandPrince();
        integrator2->setPrecision(1e-16);
        integrator2->Run(E_S);
        delete integrator2;
    });


    thread_1.join();
    thread_2.join();


    ui->graphicsView->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);

    ui->progressBar->setValue(100);
    ui->progressBar_2->setValue(100);

    auto res=E_S->getResult();
    res.toText("hello.txt");
}

MainWindow::~MainWindow()
{
    delete ui;
}


void MainWindow::on_pushButton_2_clicked()
{
    QLineSeries *series=new QLineSeries();


    TMatrix Result;
    if(ui->comboBox_3->currentText()=="Earth-Moon")
    {
        Result=E_M->getResult();
    } else if(ui->comboBox_3->currentText()=="Sun-Earth")
    {
        Result=E_S->getResult();
    }

  //t,x,y,z

    QString arg[]={"t","x","y","z"};
    int first=0,second=0;

    for(int i=0;i<4;i++)
    {
     if (ui->comboBox->currentText()==arg[i])
     {
         first=i;
     }

     if (ui->comboBox_2->currentText()==arg[i])
     {
         second=i;
     }
    }


    for (int i=0;i<Result.GetRowCount();i++)
            {
                //   *series << QPointF(  Result(i,first),Result(i,second) );
                *series<<QPointF((Result(i,2)-Result(i,8))*E_S->a,(Result(i,3)-Result(i,9))*E_S->a);
            }


    //создаем график и добавляем в него синусоиду
    QChart *chart=new QChart();
    chart->addSeries(series);
    chart->legend()->hide();
    chart->setTitle("Graphics");

    //настройка осей графика
    QValueAxis *axisX=new QValueAxis();
    axisX->setTitleText(arg[first]);
    axisX->setLabelFormat("%i");
    chart->addAxis(axisX,Qt::AlignBottom);
    series->attachAxis(axisX);
    //настройка осей графика


    QValueAxis *axisY=new QValueAxis();
    axisY->setTitleText(arg[second]);
    axisY->setLabelFormat("%i");
    chart->addAxis(axisY,Qt::AlignLeft);
    series->attachAxis(axisY);

    //создаем представление графика


    ui->graphicsView->setChart(chart);


}

void MainWindow::on_MainWindow_iconSizeChanged(const QSize &iconSize)
{
   // ui->graphicsView->setSizePolicy(&iconSize);
     //   QChart *chart=new QChart();
    ui->graphicsView->setSizeIncrement(iconSize);

}
