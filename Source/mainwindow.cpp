#include "mainwindow.h"
#include "ui_mainwindow.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_buttonNewSet_clicked()
{
    ui->OGL->resetVariables();
    ui->OGL->generateSet();
    ui->OGL->update();
}

void MainWindow::on_buttonShowCCW_clicked()
{
    ui->OGL->showCCW = !ui->OGL->showCCW;
    ui->OGL->update();
}

void MainWindow::on_buttonShowCH_clicked()
{
    ui->OGL->showCH = !ui->OGL->showCH;
    ui->OGL->update();
}

void MainWindow::on_buttonShowHOrder_clicked()
{
    ui->OGL->showHOrder = !ui->OGL->showHOrder;
    ui->OGL->update();
}

void MainWindow::on_buttonShowTri_clicked()
{
    ui->OGL->showTri = !ui->OGL->showTri;
    ui->OGL->update();
}

void MainWindow::on_buttonShowIntEdges_clicked()
{
    ui->OGL->showIntEdges = !ui->OGL->showIntEdges;
    ui->OGL->update();
}

void MainWindow::on_buttonShowDTri_clicked()
{
    ui->OGL->showDTri = !ui->OGL->showDTri;
    ui->OGL->update();
}

void MainWindow::on_buttonShowEdgeFlip_clicked()
{
    ui->OGL->showEdgeFlip = !ui->OGL->showEdgeFlip;
    ui->OGL->update();
}

void MainWindow::on_buttonShowCircumcircles_clicked()
{
    ui->OGL->showCircumcircles = !ui->OGL->showCircumcircles;
    ui->OGL->update();
}

void MainWindow::on_buttonShowCandidate_clicked()
{
    ui->OGL->showCandidate = !ui->OGL->showCandidate;
    ui->OGL->update();
}

void MainWindow::on_buttonShowLEC_clicked()
{
    ui->OGL->showLEC = !ui->OGL->showLEC;
    ui->OGL->update();
}
