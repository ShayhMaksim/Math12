#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QLineSeries>
MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    E_S=new ModelEMS(0,0.1,0.00001);
    E_M=new NormalGPZ(0,E_S->T*0.003,0.00001*E_S->T,theta,i,w,a,e,O);
    Decorator=new EMSDecorator(0,0.003,0.00001,theta,i,w,a,e,O,"Sattelite");
    DecoratorNGPZ=new EMSDecorator(0,0.003,0.00001,theta,i,w,a,e,O,"NormalGPZ");

   // E_M->setInitialPosition(theta,i,w,a,e,O);
    ui->progressBar->setRange(0,100);
    ui->progressBar->setValue((float)E_M->getT0()/E_M->getT1()*100);

    auto thread_1=std::async (std::launch::async, ([&](){
        TIntegrator * integrator=new TDormandPrince();
        integrator->setPrecision(1e-12);
        integrator->Run(E_M);
        delete integrator;
    }));

    auto thread_2=std::async(std::launch::async,([&](){
        TIntegrator * integrator2=new TDormandPrince();
        integrator2->setPrecision(1e-12);
        integrator2->Run(E_S);
        delete integrator2;
    }));

    thread_1.get();

    Decorator->Run();
    DecoratorNGPZ->Run();

    auto thread_3=std::async(std::launch::async,([&](){
        TIntegrator * integrator3=new TDormandPrince();
        integrator3->setPrecision(1e-12);
        integrator3->Run(Decorator);
        delete integrator3;
    }));


        TIntegrator * integrator4=new TDormandPrince();
        integrator4->setPrecision(1e-12);
        integrator4->Run(DecoratorNGPZ);
        delete integrator4;


    thread_2.get();
    thread_3.get();



    setDialog(theta*180/M_PI,i*180/M_PI,w*180/M_PI,a,e,O*180/M_PI);
}


void MainWindow::setDialog(long double theta, long double i, long double w, long double a, long double e, long double O)
{
    std::string string[6];
    string[0]="theta: "+std::to_string(theta);
    string[1]="i: "+std::to_string(i);
    string[2]="w: "+std::to_string(w);
    string[3]="a: "+std::to_string(a);
    string[4]="e: "+std::to_string(e);
    string[5]="O: "+std::to_string(O);

    std::string main_text=string[0]+"\n"+
                          string[1]+"\n"+
            string[2]+"\n"+
            string[3]+"\n"+
            string[4]+"\n"+
            string[5]+"\n";
    ui->textBrowser->setText(QString::fromStdString(main_text));

}


MainWindow::~MainWindow()
{
    delete ui;
}


void MainWindow::setSeries(QLineSeries *series,TMatrix Result,int first,int second,int j,int k)
{
    for (int i=0;i<Result.GetRowCount();i++)
            {
                if (first==0){
                  *series<<QPointF((Result(i,first)),((Result(i,second+j)-Result(i,second+k))*E_S->a*1e-8));
                   ui->progressBar->setValue(10+70*(double)i/Result.GetRowCount());
                }
                else if (second==0)
                {
                      *series<<QPointF(((Result(i,first+j)-Result(i,first+k))*E_S->a)*1e-8,(Result(i,second)));
                        ui->progressBar->setValue(10+70*(double)i/Result.GetRowCount());
                }
                else {
                    long double res1=((Result(i,first+j)-Result(i,first+k))*E_S->a*1e-8);
                    long double res2=((Result(i,second+j)-Result(i,second+k))*E_S->a*1e-8);
                    *series<<QPointF(res1,res2);
                    ui->progressBar->setValue(10+70*(double)i/Result.GetRowCount());
                }
            }
}


void MainWindow::on_pushButton_2_clicked()
{
    ui->graphicsView->clearFocus();
    std::future<void> thread;

    ui->progressBar->setValue(0);

    if((ui->checkBox->isChecked())&&(ui->lineEdit->text().isEmpty())) return;


    if((ui->checkBox->isChecked()))
    {

        if(!ui->lineEdit_2->text().isEmpty()){
            if ((ui->lineEdit_2->text().toDouble()>=0)&&(ui->lineEdit_2->text().toDouble()<360))
                theta=ui->lineEdit_2->text().toDouble()*M_PI/180;
        }
        if(!ui->lineEdit_3->text().isEmpty()){
            if ((ui->lineEdit_3->text().toDouble()>=0)&&(ui->lineEdit_3->text().toDouble()<=90))
                i=ui->lineEdit_3->text().toDouble()*M_PI/180;
        }
        if(!ui->lineEdit_4->text().isEmpty()){
            if ((ui->lineEdit_4->text().toDouble()>=-90)&&(ui->lineEdit_4->text().toDouble()<=90))
                w=ui->lineEdit_4->text().toDouble()*M_PI/180;
        }
        if(!ui->lineEdit_5->text().isEmpty()){
            if ((ui->lineEdit_5->text().toDouble()>=6500000)&&(ui->lineEdit_5->text().toDouble()<=50000000))
                a=ui->lineEdit_5->text().toDouble();
        }
        if(!ui->lineEdit_6->text().isEmpty()){
            if ((ui->lineEdit_6->text().toDouble()>=0)&&(ui->lineEdit_6->text().toDouble()<=0.5))
                e=ui->lineEdit_6->text().toDouble();
        }
        if(!ui->lineEdit_7->text().isEmpty()){
            if ((ui->lineEdit_7->text().toDouble()>=0)&&(ui->lineEdit_7->text().toDouble()<360))
                O=ui->lineEdit_7->text().toDouble()*M_PI/180;
        }

        setDialog(theta*180/M_PI,i*180/M_PI,w*180/M_PI,a,e,O*180/M_PI);

        if(ui->comboBox_3->currentText()=="Earth-Sattelite")
        {
            delete E_M;
             E_M=new NormalGPZ(0,ui->lineEdit->text().toDouble()*E_S->T,0.001*ui->lineEdit->text().toDouble()*E_S->T,theta,i,w,a,e,O);
        //     E_M->setInitialPosition(0,M_PI_4,0,10000000,0,0);

            thread=std::async (std::launch::async, ([&](){
                TIntegrator * integrator=new TDormandPrince();
                integrator->setPrecision(1e-12);
                integrator->Run(E_M);
                delete integrator;
            }));
        }
        else if((ui->comboBox_3->currentText()=="Earth-Moon")
                ||(ui->comboBox_3->currentText()=="Sun-Earth")
                ||(ui->comboBox_3->currentText()=="Sun-Moon"))
        {
            delete E_S;
            E_S=new ModelEMS(0,ui->lineEdit->text().toDouble(),0.001*ui->lineEdit->text().toDouble());
            thread=std::async(std::launch::async,([&](){
                TIntegrator * integrator2=new TDormandPrince();
                integrator2->setPrecision(1e-12);
                integrator2->Run(E_S);
                delete integrator2;
            }));
        } else if (ui->comboBox_3->currentText()=="Earth-Sattelite(CGP Decorator)")
        {        
            thread=std::async(std::launch::async,([&](){
             //   delete Decorator;
                Decorator=new EMSDecorator(0,ui->lineEdit->text().toDouble(),0.001*ui->lineEdit->text().toDouble(),theta,i,w,a,e,O,"Sattelite");
                TIntegrator * integrator3=new TDormandPrince();
                integrator3->setPrecision(1e-12);
                Decorator->Run();
                integrator3->Run(Decorator);
                delete integrator3;
            }));
        } else if (ui->comboBox_3->currentText()=="Earth-Sattelite(NGP Decorator)")
        {
            thread=std::async(std::launch::async,([&](){
               // delete DecoratorNGPZ;
                DecoratorNGPZ=new EMSDecorator(0,ui->lineEdit->text().toDouble(),0.001*ui->lineEdit->text().toDouble(),theta,i,w,a,e,O,"NormalGPZ");
                TIntegrator * integrator3=new TDormandPrince();
                integrator3->setPrecision(1e-12);
                DecoratorNGPZ->Run();
                integrator3->Run(DecoratorNGPZ);
                delete integrator3;
            }));
        }
    }


    QLineSeries *series=new QLineSeries();

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
    ui->progressBar->setValue(10);

    TMatrix Result;


    if(ui->checkBox->isChecked()){
        thread.get();
    }

    ui->progressBar->setValue(50);

    if(ui->comboBox_3->currentText()=="Earth-Sattelite")
    {
        Result=E_M->getResult();
        for (int i=0;i<Result.GetRowCount();i++)
                {
                       *series << QPointF(  (qreal)Result(i,first),(qreal)Result(i,second) );
                       ui->progressBar->setValue(50+30*(double)i/Result.GetRowCount());
                }
    }
    else if(ui->comboBox_3->currentText()=="Earth-Moon")
    {
             setSeries(series,E_S->getResult(),first,second,0,3);
    } else if (ui->comboBox_3->currentText()=="Sun-Earth"){
        setSeries(series,E_S->getResult(),first,second,3,6);
    } else if(ui->comboBox_3->currentText()=="Sun-Moon"){
             setSeries(series,E_S->getResult(),first,second,0,6);
    } else if(ui->comboBox_3->currentText()=="Earth-Sattelite(CGP Decorator)")
    {
        Result=Decorator->getResult();
        for (int i=0;i<Result.GetRowCount();i++)
                {
                       *series << QPointF(  (qreal)Result(i,first),(qreal)Result(i,second) );
                       ui->progressBar->setValue(50+30*(double)i/Result.GetRowCount());
                }
    } else if(ui->comboBox_3->currentText()=="Earth-Sattelite(NGP Decorator)")
    {
        Result=DecoratorNGPZ->getResult();
        for (int i=0;i<Result.GetRowCount();i++)
                {
                       *series << QPointF(  (qreal)Result(i,first),(qreal)Result(i,second) );
                       ui->progressBar->setValue(50+30*(double)i/Result.GetRowCount());
                }
    }


    ui->progressBar->setValue(80);



    //создаем график и добавляем в него синусоиду
    QChart *chart=new QChart();
    chart->addSeries(series);
    chart->legend()->hide();
    if ((ui->comboBox_3->currentText()=="Sun-Earth")
            ||(ui->comboBox_3->currentText()=="Sun-Moon")
            ||(ui->comboBox_3->currentText()=="Earth-Moon"))
        chart->setTitle("График: масштаб 1e+8");
    else
        chart->setTitle("График");

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

    ui->graphicsView->setRubberBand(QChartView::RectangleRubberBand);

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
