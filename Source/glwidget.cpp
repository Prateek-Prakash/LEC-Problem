#include "glwidget.h"

GLWidget::GLWidget(QWidget *parent) : QOpenGLWidget(parent)
{
    resetVariables();
}

void GLWidget::initializeGL()
{
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_LINE_SMOOTH);
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);

    glClearColor(0.0, 0.0, 0.0, 1.0);
    glPointSize(POINT_SIZE);

    generateSet();
}

void GLWidget::paintGL()
{
    glClear(GL_COLOR_BUFFER_BIT);

    if(showCCW)
    {
        plotCCW();
    }

    if(showHOrder)
    {
        plotHOrder();
    }

    if(showCH)
    {
        plotCH();
    }

    if(showTri)
    {
        plotTri();
    }

    if(showDTri)
    {
        plotDTri();
    }

    if(showIntEdges)
    {
        plotIntEdges();
    }

    if(showEdgeFlip)
    {
        plotEdgeFlip();
    }

    if(showCircumcircles)
    {
        plotCircumcircles();
    }

    if(showCandidate)
    {
        plotCandidate();
    }

    if(showLEC)
    {
        plotLEC();
    }

    plotSet();
}

void GLWidget::resizeGL(int w, int h)
{
    (void) w, h;
}

void GLWidget::resetVariables()
{

    showCCW = false;
    showCH = false;
    showHOrder = false;
    showTri = false;
    showIntEdges = false;
    showDTri = false;
    showEdgeFlip = false;
    showCircumcircles = false;
    showCandidate = false;
    showLEC = false;
    indexCCW = 0;
    indexCH = 0;
    indexTri = 0;
    indexEdge = 0;
    indexCand = 0;
    indexLEC = 0;
}

void GLWidget::generateSet()
{
    srand(time(NULL));
    for(int i = 0; i < SET; i++)
    {
        PT[i].x = (-90 + (rand() % (90 - (-90)))) / 100.0;
        PT[i].y = (-90 + (rand() % (90 - (-90)))) / 100.0;
    }

    /*
    PT[0].x = 0.170000; PT[0].y = 0.780000;
    PT[1].x = -0.650000; PT[1].y = -0.080000;
    PT[2].x = -0.290000; PT[2].y = 0.140000;
    PT[3].x = 0.540000; PT[3].y = -0.180000;
    PT[4].x = 0.680000; PT[4].y = -0.860000;
    PT[5].x = -0.400000; PT[5].y = 0.810000;
    */

    generateCH();
}

void GLWidget::plotSet()
{
    glColor3d(1.0, 1.0, 1.0);
    glBegin(GL_POINTS);
        for(int i = 0; i < SET; i++)
        {
            glVertex2d(PT[i].x, PT[i].y);
        }
    glEnd();
}

void GLWidget::generateCH()
{
    sortSetCCW();

    int tempIndex = 0;

    iCH[indexCH++] = iCCW[tempIndex++];
    iCH[indexCH++] = iCCW[tempIndex++];

    while(tempIndex < SET)
    {
        if(isLeft(PT[iCH[indexCH - 2]].x, PT[iCH[indexCH - 1]].x, PT[iCCW[tempIndex]].x, PT[iCH[indexCH - 2]].y, PT[iCH[indexCH - 1]].y, PT[iCCW[tempIndex]].y) >= 0)
        {
            iCH[indexCH++] = iCCW[tempIndex++];
        }
        else
        {
            iCH[indexCH--] = iCCW[tempIndex];
        }
    }

    generateDTri();
}

void GLWidget::sortSetCCW()
{
    // Point 0
    iCCW[0] = findMinY();
    indexCCW++;

    double previousAngle = 0.0;
    double lowestAngle = 360.0;
    int tempIndex = 0.0;

    // Point N
    for(int j = 1; j < SET; j++)
    {
        lowestAngle = 360.0;
        tempIndex = 0;
        for(int i = 0; i < SET; i++)
        {
            double tempAngle = atan2(PT[iCCW[0]].y - PT[i].y, PT[iCCW[0]].x - PT[i].x);
            if(j == 1 && tempAngle <= lowestAngle)
            {
                lowestAngle = tempAngle;
                tempIndex = i;
            }
            else if(tempAngle <= lowestAngle && tempAngle > previousAngle)
            {
                lowestAngle = tempAngle;
                tempIndex = i;
            }
        }
        previousAngle = lowestAngle;
        iCCW[j] = tempIndex;
        indexCCW++;
    }
}

int GLWidget::findMinY()
{
    // Point 0
    double tempY = 100.0;
    double tempX = -100.0;
    int index = 0;
    for(int i = 0; i < SET; i++)
    {
        if(PT[i].y <= tempY)
        {
            if(PT[i].y == tempY)
            {
                if(PT[i].x > tempX)
                {
                    index = i;
                    tempY = PT[i].y;
                    tempX = PT[i].x;
                }
            }
            else
            {
                index = i;
                tempY = PT[i].y;
                tempX = PT[i].x;
            }
        }
    }
    return index;
}

void GLWidget::plotCCW()
{
    bool color = true;
    double x = 0.0;
    for(int i = 1; i < SET; i++)
    {
        if(color)
        {
            x = 1.0;
        }
        else
        {
            x = 0.0;
        }
        color = !color;
        glColor3d(1.0, x, 0.0);
        glBegin(GL_LINES);
            glVertex2d(PT[iCCW[0]].x, PT[iCCW[0]].y);
            glVertex2d(PT[iCCW[i]].x, PT[iCCW[i]].y);
        glEnd();
    }
}

void GLWidget::plotEdgeFlip()
{
    glColor3d(1.0, 0.5, 0.0);
    for(int i = 0; i < indexEdge; i++)
    {
        glBegin(GL_LINES);
            glVertex2d(PT[DED[i].p1].x, PT[DED[i].p1].y);
            glVertex2d(PT[DED[i].p2].x, PT[DED[i].p2].y);
        glEnd();
    }
}

void GLWidget::plotDTri()
{
    glColor3d(1.0, 0.0, 1.0);
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    glBegin(GL_TRIANGLES);
    for(int i = 0; i < indexTri; i++)
    {
        glVertex2d(PT[DTRI[i].p1].x, PT[DTRI[i].p1].y);
        glVertex2d(PT[DTRI[i].p2].x, PT[DTRI[i].p2].y);
        glVertex2d(PT[DTRI[i].p3].x, PT[DTRI[i].p3].y);
    }
    glEnd();
}

double GLWidget::isLeft(double xP1, double xP2, double xP3, double yP1, double yP2, double yP3)
{
    return ((xP2 - xP1) * (yP3 - yP1)) - ((xP3 - xP1) * (yP2 - yP1));
}

void GLWidget::plotCH()
{
    glColor3d(0.0, 1.0, 0.0);
    glBegin(GL_LINE_LOOP);
    for(int i = 0; i < indexCH; i++)
    {
        glVertex2d(PT[iCH[i]].x, PT[iCH[i]].y);
    }
    glEnd();
}

void GLWidget::generateDTri()
{
    generateHOrder();
    generateRandTri();
    findIntEdges();
    flipEdges();

    generateCC();
}

void GLWidget::printInfo()
{
    qDebug("--\nPrinting All Information\n--\nPT\n--");
    for(int i = 0; i < SET; i++)
    {
        qDebug("PT %d: %f, %f", i, PT[i].x, PT[i].y);
    }
    qDebug("--\nTRI\n--");
    for(int i = 0; i < indexTri; i++)
    {
        qDebug("TRI %d: %d, %d, %d", i, TRI[i].p1, TRI[i].p2, TRI[i].p3);
    }
    qDebug("--\nED\n--");
    for(int i = 0; i < indexTri; i++)
    {
        qDebug("ED %d: %d, %d", i, ED[i].p1, ED[i].p2);
    }

}

void GLWidget::flipEdges()
{
    //printInfo();

    for(int i = 0; i < indexEdge; i++)
    {
        int tri1, tri2, C, D;
        if(DED[i].isD == false)
        {
            getEdgeTriangles(DED[i], tri1, tri2, C, D);
            int A = DED[i].p1;
            int B = DED[i].p2;

            //qDebug("--\nDTRI 1: %d, %d, %d", DTRI[tri1].p1, DTRI[tri1].p2, DTRI[tri1].p3);
            //qDebug("DTRI 2: %d, %d, %d", DTRI[tri2].p1, DTRI[tri2].p2, DTRI[tri2].p3);
            //qDebug("DED %d: %d, %d", i, DED[i].p1, DED[i].p2);
            //qDebug("%d, %d, %d, %d", A, B, C, D);

            if(!isDelaunayEdge(PT[A].x, PT[A].y, PT[B].x, PT[B].y, PT[C].x, PT[C].y, PT[D].x, PT[D].y))
            {
                //qDebug("Flipped!");
                DED[i].p1 = C;
                DED[i].p2 = D;
                DTRI[tri1].p1 = A;
                DTRI[tri1].p2 = D;
                DTRI[tri1].p3 = C;
                DTRI[tri2].p1 = B;
                DTRI[tri2].p2 = C;
                DTRI[tri2].p3 = D;
                DED[i].isD = true;
                i = 0;
                //qDebug("New DTRI 1: %d, %d, %d", DTRI[tri1].p1, DTRI[tri1].p2, DTRI[tri1].p3);
                //qDebug("New DTRI 2: %d, %d, %d", DTRI[tri2].p1, DTRI[tri2].p2, DTRI[tri2].p3);
                //qDebug("New DED %d: %d, %d", i, DED[i].p1, DED[i].p2);
            }
        }

    }
}

void GLWidget::getEdgeTriangles(EDGE E, int &t1, int &t2, int &C, int &D)
{
    int A = E.p1;
    int B = E.p2;
    int count = 0;
    for(int i = 0; i < indexTri; i++)
    {
        if((A == DTRI[i].p1 || A == DTRI[i].p2 || A == DTRI[i].p3)
                && (B == DTRI[i].p1 || B == DTRI[i].p2 || B == DTRI[i].p3))
        {
            if(count++ == 0)
            {
                t1 = i;
            }
            else
            {
                t2 = i;
            }
        }
    }

    // Find C
    // TRI 1
    if(A == DTRI[t1].p1 && B == DTRI[t1].p2)
    {
        C = DTRI[t1].p3;
    }
    else if(A == DTRI[t1].p2 && B == DTRI[t1].p3)
    {
        C = DTRI[t1].p1;
    }
    else if(A == DTRI[t1].p3 && B == DTRI[t1].p1)
    {
        C = DTRI[t1].p2;
    }
    // TRI 2
    else if(A == DTRI[t2].p1 && B == DTRI[t2].p2)
    {
        C = DTRI[t2].p3;
    }
    else if(A == DTRI[t2].p2 && B == DTRI[t2].p3)
    {
        C = DTRI[t2].p1;
    }
    else if(A == DTRI[t2].p3 && B == DTRI[t2].p1)
    {
        C = DTRI[t2].p2;
    }

    // Find D
    // TRI 1
    if(DTRI[t1].p1 != A && DTRI[t1].p1 != B && DTRI[t1].p1 != C)
    {
        D = DTRI[t1].p1;
    }
    else if(DTRI[t1].p2 != A && DTRI[t1].p2 != B && DTRI[t1].p2 != C)
    {
        D = DTRI[t1].p2;
    }
    else if(DTRI[t1].p3 != A && DTRI[t1].p3 != B && DTRI[t1].p3 != C)
    {
        D = DTRI[t1].p3;
    }
    // TRI 2
    else if(DTRI[t2].p1 != A && DTRI[t2].p1 != B && DTRI[t2].p1 != C)
    {
        D = DTRI[t2].p1;
    }
    else if(DTRI[t2].p2 != A && DTRI[t2].p2 != B && DTRI[t2].p2 != C)
    {
        D = DTRI[t2].p2;
    }
    else if(DTRI[t2].p3 != A && DTRI[t2].p3 != B && DTRI[t2].p3 != C)
    {
        D = DTRI[t2].p3;
    }
}

void GLWidget::findIntEdges()
{
    for(int i = 0; i < indexTri; i++)
    {
        if(onCH(TRI[i].p1) && onCH(TRI[i].p2))
        {
            if(qAbs(getCHIndex(TRI[i].p1) - getCHIndex(TRI[i].p2)) < (indexCH - 1)
                    && qAbs(getCHIndex(TRI[i].p1) - getCHIndex(TRI[i].p2)) > 1)
            {
                addEdge(TRI[i].p1, TRI[i].p2);
            }
        }
        else
        {
            addEdge(TRI[i].p1, TRI[i].p2);
        }

        if(onCH(TRI[i].p2) && onCH(TRI[i].p3))
        {
            if(qAbs(getCHIndex(TRI[i].p2) - getCHIndex(TRI[i].p3)) < (indexCH - 1)
                    && qAbs(getCHIndex(TRI[i].p2) - getCHIndex(TRI[i].p3)) > 1)
            {
                addEdge(TRI[i].p2, TRI[i].p3);
            }
        }
        else
        {
            addEdge(TRI[i].p2, TRI[i].p3);
        }

        if(onCH(TRI[i].p3) && onCH(TRI[i].p1))
        {
            if(qAbs(getCHIndex(TRI[i].p3) - getCHIndex(TRI[i].p1)) < (indexCH - 1)
                    && qAbs(getCHIndex(TRI[i].p3) - getCHIndex(TRI[i].p1)) > 1)
            {
                addEdge(TRI[i].p3, TRI[i].p1);
            }
        }
        else
        {
            addEdge(TRI[i].p3, TRI[i].p1);
        }
    }

    // Make Copy
    for(int i = 0; i < indexEdge; i++)
    {
        DED[i].p1 = ED[i].p1;
        DED[i].p2 = ED[i].p2;
        DED[i].isD = ED[i].isD;
    }
}

int GLWidget::getCHIndex(int v)
{
    for(int i = 0; i < indexCH; i++)
    {
        if(v == iCH[i])
            return i;
    }
    return 100;
}

void GLWidget::addEdge(int v1, int v2)
{
    if(indexEdge == 0)
    {
        ED[indexEdge].p1 = v1;
        ED[indexEdge++].p2 = v2;
        return;
    }
    bool edgeExists = false;
    for(int i = 0; i < indexEdge; i ++)
    {
        if((v1 == ED[i].p1 && v2 == ED[i].p2) || (v2 == ED[i].p1 && v1 == ED[i].p2))
        {
            edgeExists = true;
        }
    }
    if(!edgeExists)
    {
        ED[indexEdge].p1 = v1;
        ED[indexEdge++].p2 = v2;
    }
}

void GLWidget::plotIntEdges()
{
    glColor3d(0.0, 0.5, 1.0);
    for(int i = 0; i < indexEdge; i++)
    {
        glBegin(GL_LINES);
            glVertex2d(PT[ED[i].p1].x, PT[ED[i].p1].y);
            glVertex2d(PT[ED[i].p2].x, PT[ED[i].p2].y);
        glEnd();
    }
}

void GLWidget::generateRandTri()
{
    for(int i = 0; i < indexCH - 2; i++)
    {
        TRI[indexTri].p1 = iCH[0];
        TRI[indexTri].p2 = iCH[i + 1];
        TRI[indexTri++].p3 = iCH[i + 2];
    }
    for(int i = 0; i < SET; i++)
    {
        if(onCH(iHO[i]))
        {
            // Convex Hull Point
        }
        else
        {
            for(int j = 0; j < indexTri; j++)
            {
                if(inTriangle(PT[TRI[j].p1].x, PT[TRI[j].p1].y, PT[TRI[j].p2].x, PT[TRI[j].p2].y, PT[TRI[j].p3].x, PT[TRI[j].p3].y, PT[iHO[i]].x, PT[iHO[i]].y))
                {
                    int tempIndex = indexTri - 1;
                    indexTri += 2;

                    // Divide Outer Triangle
                    int p1 = TRI[j].p1;
                    int p2 = TRI[j].p2;
                    int p3 = TRI[j].p3;
                    int p = iHO[i];

                    while(tempIndex > j)
                    {
                        TRI[tempIndex + 2] = TRI[tempIndex--];
                    }

                    TRI[j].p3 = p;

                    TRI[j + 1].p1 = p2;
                    TRI[j + 1].p2 = p3;
                    TRI[j + 1].p3 = p;

                    TRI[j + 2].p1 = p3;
                    TRI[j + 2].p2 = p1;
                    TRI[j + 2].p3 = p;

                    // Exit For Loop
                    break;
                }
            }
        }
    }

    // Make Copy
    for(int i = 0; i < indexTri; i++)
    {
        DTRI[i].p1 = TRI[i].p1;
        DTRI[i].p2 = TRI[i].p2;
        DTRI[i].p3 = TRI[i].p3;
    }
}

void GLWidget::generateHOrder()
{
    double previousX = -100;
    int lowestIndex = 0;

    for(int j = 0; j < SET; j++)
    {
        double lowestNum = 100;
        for(int i = 0; i < SET; i++)
        {
            if(PT[i].x < lowestNum && PT[i].x > previousX)
            {
                lowestNum = PT[i].x;
                lowestIndex = i;
            }
            else if (PT[i].x < lowestNum && PT[i].x == previousX)
            {
                // Handle Same X Values
                int oldPoint = false;
                for(int k = 0; k < j; k++)
                {
                    if(PT[i].x == PT[iHO[k]].x && PT[i].y == PT[iHO[k]].y)
                    {
                        oldPoint = true;
                    }
                }
                if(!oldPoint)
                {
                    lowestNum = PT[i].x;
                    lowestIndex = i;
                }
            }
        }

        previousX = lowestNum;
        iHO[j] = lowestIndex;
        iHO[j] = lowestIndex;
    }
}

void GLWidget::plotHOrder()
{
    glColor3d(0.0, 0.0, 1.0);
    glBegin(GL_LINE_STRIP);
    for(int i = 0; i < SET; i++)
    {
        glVertex2d(PT[iHO[i]].x, PT[iHO[i]].y);
    }
    glEnd();
}

void GLWidget::plotTri()
{
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    glColor3d(0.5, 0.5, 0.5);
    glBegin(GL_TRIANGLES);
    for(int i = 0; i < indexTri; i++)
    {
        glVertex2d(PT[TRI[i].p1].x, PT[TRI[i].p1].y);
        glVertex2d(PT[TRI[i].p2].x, PT[TRI[i].p2].y);
        glVertex2d(PT[TRI[i].p3].x, PT[TRI[i].p3].y);
    }
    glEnd();
}

bool GLWidget::onCH(int V)
{
    for(int i = 0; i < indexCH; i++)
    {
        if(V == iCH[i])
            return true;
    }
    return false;
}

bool GLWidget::inTriangle(double xP1, double yP1, double xP2, double yP2, double xP3, double yP3, double X, double Y)
{
    // Barycentric Coordinates
    float alpha = ((yP2 - yP3)*(X - xP3) + (xP3 - xP2)*(Y - yP3)) /
            ((yP2 - yP3)*(xP1 - xP3) + (xP3 - xP2)*(yP1 - yP3));
    float beta = ((yP3 - yP1)*(X - xP3) + (xP1 - xP3)*(Y - yP3)) /
           ((yP2 - yP3)*(xP1 - xP3) + (xP3 - xP2)*(yP1 - yP3));
    float gamma = 1.0f - alpha - beta;

    if(alpha > 0 && beta > 0 && gamma > 0)
    {
        return true;
    }
    else
    {
        return false;
    }
}

bool GLWidget::isDelaunayEdge(double Ax, double Ay, double Bx, double By, double Cx, double Cy, double Dx, double Dy)
{
    QMatrix4x4 tempMatrix
            (
                Ax, Ay, (Ax*Ax) + (Ay*Ay), 1,
                Bx, By, (Bx*Bx) + (By*By), 1,
                Cx, Cy, (Cx*Cx) + (Cy*Cy), 1,
                Dx, Dy, (Dx*Dx) + (Dy*Dy), 1
            );
    if(tempMatrix.determinant() > 0)
    {
        return false;
    }
    else
    {
        return true;
    }
}

void GLWidget::generateCC()
{
    for(int i = 0; i < indexTri; i++)
    {
        int A = DTRI[i].p1;
        int B = DTRI[i].p2;
        int C = DTRI[i].p3;

        // First Bisector
        double AxB = (PT[A].x + PT[B].x) / 2.0;
        double AyB = (PT[A].y + PT[B].y) / 2.0;

        double nAB = (PT[B].y - PT[A].y);
        double dAB = (PT[B].x - PT[A].x);
        double mAB = -(dAB / nAB);

        double bAB = AyB - (mAB * AxB);

        // Second Bisector
        double BxC = (PT[B].x + PT[C].x) / 2.0;
        double ByC = (PT[B].y + PT[C].y) / 2.0;

        double nBC = (PT[C].y - PT[B].y);
        double dBC = (PT[C].x - PT[B].x);
        double mBC = -(dBC / nBC);

        double bBC = ByC - (mBC * BxC);

        // X
        double x = (bBC - bAB) / (mAB - mBC);
        double y = (mAB * x) + bAB;

        CC[i].x = x;
        CC[i].y = y;
        CC[i].r = qSqrt(((CC[i].x - PT[DTRI[i].p1].x)*(CC[i].x - PT[DTRI[i].p1].x)) + ((CC[i].y - PT[DTRI[i].p1].y)*(CC[i].y - PT[DTRI[i].p1].y)));
        // qDebug("CC %d: %f, %f", i, x, y);
    }

    generateCand();
}

void GLWidget::plotCircumcircles()
{
    glColor3d(0.5, 0.5, 0.0);
    for(int i = 0; i < indexTri; i++)
    {
        glBegin(GL_POINTS);
            glVertex2d(CC[i].x, CC[i].y);
        glEnd();
        drawCircle(CC[i].x, CC[i].y, CC[i].r, 100);
    }
}

void GLWidget::drawCircle(float cx, float cy, float r, int num_segments)
{
    glBegin(GL_LINE_LOOP);
    for (int ii = 0; ii < num_segments; ii++)
    {
        float theta = 2.0f * 3.1415926f * float(ii) / float(num_segments);
        float x = r * cosf(theta);
        float y = r * sinf(theta);
        glVertex2f(x + cx, y + cy);
    }
    glEnd();
}

void GLWidget::generateCand()
{
    // Check Circumcenters
    for(int i = 0; i < indexTri; i++)
    {
        bool isCandidate = false;
        for(int j = 0; j < indexTri; j++)
        {
            if(inTriangle(PT[DTRI[j].p1].x, PT[DTRI[j].p1].y, PT[DTRI[j].p2].x, PT[DTRI[j].p2].y, PT[DTRI[j].p3].x, PT[DTRI[j].p3].y, CC[i].x, CC[i].y))
            {
                isCandidate = true;
            }
        }
        if(isCandidate)
        {
            CAND[indexCand].x = CC[i].x;
            CAND[indexCand].y = CC[i].y;
            CAND[indexCand++].r = CC[i].r;
        }
    }

    // Check CH + Voronoi Intersect
    for(int i = 0; i < indexCH; i++)
    {
        int A = iCH[i];
        int B;
        if(i < indexCH - 1)
        {
            B = iCH[i + 1];
        }
        else
        {
            B = iCH[0];
        }

        double x = (PT[A].x + PT[B].x) / 2.0;
        double y = (PT[A].y + PT[B].y) / 2.0;

        CAND[indexCand].x = x;
        CAND[indexCand].y = y;
        CAND[indexCand++].r = qSqrt(((x - PT[A].x)*(x - PT[A].x)) + ((y - PT[A].y)*(y - PT[A].y)));
    }

    // Determine LEC
    determineLEC();
}

void GLWidget::plotCandidate()
{
    glColor3d(0.5, 0.1, 0.0);
    glBegin(GL_POINTS);
    for(int i = 0; i <indexCand; i++)
    {
        glBegin(GL_POINTS);
            glVertex2d(CAND[i].x, CAND[i].y);
        glEnd();
        drawCircle(CAND[i].x, CAND[i].y, CAND[i].r, 100);
    }
    glEnd();
}

void GLWidget::determineLEC()
{
    LEC[indexLEC].x = 0.0;
    LEC[indexLEC].y = 0.0;
    LEC[indexLEC++].r = 0.0;

    for(int i = 0; i < indexCand; i++)
    {
        if(LEC[indexLEC - 1].r < CAND[i].r)
        {
            // qDebug("%d: %f < %f", i, LEC[indexLEC - 1].r, CAND[i].r);
            if(isEmptyCircle(CAND[i].x, CAND[i].y, CAND[i].r))
            {
                indexLEC = 0;
                LEC[indexLEC].x = CAND[i].x;
                LEC[indexLEC].y = CAND[i].y;
                LEC[indexLEC++].r = CAND[i].r;
            }
            else
            {
                // qDebug("Not Empty Circle!");
            }
        }
        else if(LEC[indexLEC - 1].r == CAND[i].r)
        {
            // qDebug("%d: %f = %f", i, LEC[indexLEC - 1].r, CAND[i].r);
            if(isEmptyCircle(CAND[i].x, CAND[i].y, CAND[i].r))
            {
                LEC[indexLEC].x = CAND[i].x;
                LEC[indexLEC].y = CAND[i].y;
                LEC[indexLEC++].r = CAND[i].r;
            }
            else
            {
                // qDebug("Not Empty Circle!");
            }

        }
        else
        {
            // qDebug("%d: %f > %f", i, LEC[indexLEC - 1].r, CAND[i].r);
        }
    }
}

bool GLWidget::isEmptyCircle(double x, double y, double radius)
{
    int count = 0;
    for(int i = 0; i < SET; i++)
    {
        double tempRadius = qSqrt(((x - PT[i].x)*(x - PT[i].x)) + ((y - PT[i].y)*(y - PT[i].y)));
        if(tempRadius  < radius)
        {
            count++;
        }
    }

    if(count > 0)
    {
        return false;
    }

    return true;
}

void GLWidget::plotLEC()
{
    glColor3d(0.45, 0.85, 0.75);
    for(int i = 0; i < indexLEC; i++)
    {
        glBegin(GL_POINTS);
            glVertex2d(LEC[i].x, LEC[i].y);
        glEnd();
        drawCircle(LEC[i].x, LEC[i].y, LEC[i].r, 100);
    }
}
