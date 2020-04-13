#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private slots:
    void on_buttonNewSet_clicked();
    void on_buttonShowCCW_clicked();
    void on_buttonShowCH_clicked();
    void on_buttonShowHOrder_clicked();
    void on_buttonShowTri_clicked();
    void on_buttonShowIntEdges_clicked();
    void on_buttonShowDTri_clicked();
    void on_buttonShowEdgeFlip_clicked();
    void on_buttonShowCircumcircles_clicked();

    void on_buttonShowCandidate_clicked();

    void on_buttonShowLEC_clicked();

private:
    Ui::MainWindow *ui;
};

#endif // MAINWINDOW_H
