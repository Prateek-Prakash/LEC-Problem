#ifndef GLWIDGET_H
#define GLWIDGET_H

#define SET 10
#define POINT_SIZE 3.0

#include <QOpenGLWidget>
#include <ctime>
#include <QMatrix4x4>
#include <QtMath>
#include <QDebug>

struct VERTEX
{
  double x, y, r;
};

struct TRIANGLE
{
    int p1, p2, p3;
};

struct EDGE
{
    int p1, p2;
    bool isD = false;
};

class GLWidget : public QOpenGLWidget
{
    Q_OBJECT

public:
    explicit GLWidget(QWidget *parent = 0);

    void initializeGL();
    void paintGL();
    void resizeGL(int w, int h);

    VERTEX PT[SET];

    bool showCCW;
    bool showCH;
    bool showHOrder;
    bool showTri;
    bool showIntEdges;
    bool showDTri;
    bool showEdgeFlip;
    bool showCircumcircles;
    bool showCandidate;
    bool showLEC;

    int indexCCW;
    int iCCW[SET];

    int indexCH;
    int iCH[SET];

    int iHO[SET];

    int indexTri;
    TRIANGLE TRI[SET * SET];
    TRIANGLE DTRI[SET * SET];

    int indexEdge;
    EDGE ED[SET * SET];
    EDGE DED[SET * SET];

    void resetVariables();
    void generateSet();
    void plotSet();
    void generateCH();
    void sortSetCCW();
    int findMinY();
    void plotCCW();
    double isLeft(double xP1, double xP2, double xP3, double yP1, double yP2, double yP3);
    void plotCH();
    void generateDTri();
    void findIntEdges();
    void plotIntEdges();
    void addEdge(int v1, int v2);
    void generateRandTri();
    void generateHOrder();
    void plotHOrder();
    void plotTri();
    void plotDTri();
    void plotEdgeFlip();

    void flipEdges();
    void getEdgeTriangles(EDGE E, int &t1, int &t2, int &C, int &D);

    int getCHIndex(int v);

    void printInfo();

    bool onCH(int V);
    bool inTriangle(double xP1, double xP2, double xP3, double yP1, double yP2, double yP3, double X, double Y);
    bool isDelaunayEdge(double Ax, double Ay, double Bx, double By, double Cx, double Cy, double Dx, double Dy);

    VERTEX CC[SET * SET];
    void generateCC();
    void plotCircumcircles();
    void drawCircle(float cx, float cy, float r, int num_segments);
    void plotCandidate();

    VERTEX CAND[SET * SET];
    int indexCand;
    void generateCand();

    int indexLEC;
    VERTEX LEC[SET];
    void determineLEC();
    bool isEmptyCircle(double x, double y, double radius);
    void plotLEC();
};

#endif // GLWIDGET_H
