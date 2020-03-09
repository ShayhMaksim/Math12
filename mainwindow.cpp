#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QLineSeries>
MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);


    E_S=new ModelEMS(0,1,0.0018);

    E_M=new SunRotation(0,50000,1);

    E_M->getInitialPosition(0,0,0,10000000,0,0);
    ui->progressBar->setRange(0,100);
    ui->progressBar->setValue((float)E_M->getT0()/E_M->getT1()*100);


    auto thread_1=std::async (std::launch::async, ([&](){
        TIntegrator * integrator=new TDormandPrince();
        integrator->setPrecision(1e-18);
        integrator->Run(E_M);
        delete integrator;
    }));

    auto thread_2=std::async(std::launch::async,([&](){
        TIntegrator * integrator2=new TDormandPrince();
        integrator2->setPrecision(1e-16);
        integrator2->Run(E_S);
        delete integrator2;
    }));


    thread_1.get();
    thread_2.get();


    ui->graphicsView->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);


}

MainWindow::~MainWindow()
{
    delete ui;
}


void MainWindow::on_pushButton_2_clicked()
{
    std::future<void> thread_1;
    std::future<void> thread_2;

    if(ui->checkBox->isChecked())
    {
        if(ui->comboBox_3->currentText()=="Earth-Sattelite")
        {
            delete E_M;
             E_M=new SunRotation(0,ui->lineEdit->text().toDouble()*86400,1);
             E_M->getInitialPosition(0,0,0,10000000,0,0);

            thread_1=std::async (std::launch::async, ([&](){
                TIntegrator * integrator=new TDormandPrince();
                integrator->setPrecision(1e-18);
                integrator->Run(E_M);
                delete integrator;
            }));
        }
        else if(ui->comboBox_3->currentText()=="Earth-Moon")
        {
            delete E_S;
            E_S=new ModelEMS(0,ui->lineEdit->text().toDouble(),0.0018);
            thread_2=std::async(std::launch::async,([&](){
                TIntegrator * integrator2=new TDormandPrince();
                integrator2->setPrecision(1e-16);
                integrator2->Run(E_S);
                delete integrator2;
            }));
        }
    }

    QLineSeries *series=new QLineSeries();

    ui->progressBar->setValue(0);
  //t,x,y,zzx
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
    ui->progressBar->setValue(10);

    TMatrix Result;

    if(ui->checkBox->isChecked())
    {
        if(ui->comboBox_3->currentText()=="Earth-Sattelite")
        {
              thread_1.get();
        }
        else if(ui->comboBox_3->currentText()=="Earth-Moon")
        {
              thread_2.get();
        }
    }

    ui->progressBar->setValue(50);

    if(ui->comboBox_3->currentText()=="Earth-Sattelite")
    {
        Result=E_M->getResult();
        for (int i=0;i<Result.GetRowCount();i++)
                {
                       *series << QPointF(  Result(i,first),Result(i,second) );
                       ui->progressBar->setValue(50+30*(double)i/Result.GetRowCount());
                }

    }
    else if(ui->comboBox_3->currentText()=="Earth-Moon")
    {
        Result=E_S->getResult();
        for (int i=0;i<Result.GetRowCount();i++)
                {
                    if (first==0){
                      *series<<QPointF((Result(i,first))*E_S->a,(Result(i,second)-Result(i,second+3))*E_S->a);
                       ui->progressBar->setValue(10+70*(double)i/Result.GetRowCount());
                    }
                    else if (second==0)
                    {
                          *series<<QPointF((Result(i,first)-Result(i,first+3))*E_S->a,(Result(i,second))*E_S->a);
                            ui->progressBar->setValue(10+70*(double)i/Result.GetRowCount());
                    }
                    else {
                    *series<<QPointF((Result(i,first)-Result(i,first+3))*E_S->a,(Result(i,second)-Result(i,second+3))*E_S->a);
                            ui->progressBar->setValue(10+70*(double)i/Result.GetRowCount());
                    }
                }
    }
    ui->progressBar->setValue(80);
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
    ui->progressBar->setValue(90);

    ui->graphicsView->setChart(chart);

    ui->progressBar->setValue(100);
}

void MainWindow::on_MainWindow_iconSizeChanged(const QSize &iconSize)
{
   // ui->graphicsView->setSizePolicy(&iconSize);
     //   QChart *chart=new QChart();
    ui->graphicsView->setSizeIncrement(iconSize);

}
