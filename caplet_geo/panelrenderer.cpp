/*
CREATED : Jan 31, 2013
MODIFIED: Feb 15, 2013
AUTHOR  : Yu-Chung Hsiao
EMAIL   : yuchsiao@mit.edu

This file is part of CAPLET.

CAPLET is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

CAPLET is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the Lesser GNU General Public License
along with CAPLET.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "panelrenderer.h"

#include "debug.h"

#include <QtOpenGL>
#include <QRadialGradient>

#include <iostream>
#include <cmath>
#include <algorithm>
#include <limits>


#define DEBUG


using namespace std;
enum {X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, X4, Y4, Z4};
enum BOUNDS {XL,XU, YL,YU, ZL, ZU};


PanelRenderer::PanelRenderer( QWidget* parent)
    : QGLWidget( parent ),
      colorScheme( BYLAYER ),
      isLoaded( false ),
      rotationX (0.0),
      rotationY (0.0),
      rotationZ (0.0),
      translationX (0.0),
      translationY (0.0),
      translationZ (-10.0),
      mouseMode (0),
      scale (1.0f),
      outlineWidth (2.0f),
      outline(true),
      colorIndex (0)
{

    //* set up colors
    int comp_saturation = 105;
    int saturation      = 230;
    int comp_saturation2 = 190;
    int saturation2     = 255;

    nColor = 10;
    color = new QColor*[nColor];
    color[0] = new QColor(comp_saturation, comp_saturation, saturation);
    color[1] = new QColor(comp_saturation, saturation, comp_saturation);
    color[2] = new QColor(saturation, comp_saturation, comp_saturation);
    color[3] = new QColor(comp_saturation, saturation, saturation);
    color[4] = new QColor(saturation, comp_saturation, saturation);
    color[5] = new QColor(saturation, saturation, comp_saturation);

    color[6] = new QColor(comp_saturation2, saturation2, saturation2);
    color[7] = new QColor(saturation2, comp_saturation2, saturation2);
    color[8] = new QColor(comp_saturation2, comp_saturation2, saturation2);
    color[9] = new QColor(saturation2, comp_saturation2, comp_saturation2);


}

PanelRenderer::~PanelRenderer(){
    for (int i=0; i<nColor; i++){
        delete color[i];
    }
    delete color;
}

void PanelRenderer::plotOutline(bool outlineFlag){
    outline = outlineFlag;
    updateGL();
}

//**
//* PanelRenderer::loadPanels
//* - int n              :
//* - float (*panels)[12]:
//* - float (*normals)[3]:
//* - float *axisCoord   : two 3D points as the boundary of the show area
bool x1Less(const Rectangle &p1, const Rectangle &p2){
    return p1.x1 < p2.x1;
}
bool x2Less(const Rectangle &p1, const Rectangle &p2){
    return p1.x2 < p2.x2;
}
bool y1Less(const Rectangle &p1, const Rectangle &p2){
    return p1.y1 < p2.y1;
}
bool y2Less(const Rectangle &p1, const Rectangle &p2){
    return p1.y2 < p2.y2;
}
bool z1Less(const Rectangle &p1, const Rectangle &p2){
    return p1.z1 < p2.z1;
}
bool z2Less(const Rectangle &p1, const Rectangle &p2){
    return p1.z2 < p2.z2;
}

void PanelRenderer::setGLColorScheme(int scheme)
{
    //* if no update, return to avoid inf loops
    if ( colorScheme == scheme ){
        return;
    }
    colorScheme = scheme;

    if ( isLoaded == true ){
        conductorFPList2RectGLForDisplay();
        updateGL();
    }
}


void PanelRenderer::loadGLRects(const ConductorFPList *condFPListPtr)
{
    this->condFPListPtr = condFPListPtr;
    conductorFPList2RectGLForDisplay();
    computeGLDisplayBoundary();
    updateGL();

    isLoaded = true;
}


void PanelRenderer::setupGLRenderProperty()
{
    renderPropertyVec.resize(rectGLForDisplay.size());

    //* set up renderproperties
    for ( unsigned int i=0; i<rectGLForDisplay.size(); ++i ){
        renderPropertyVec[i].colorIndex = i % nColor;
    }
}



void PanelRenderer::computeGLDisplayBoundary()
{
    //* set up bounds
    float xmin = numeric_limits<float>::max();
    float xmax = numeric_limits<float>::min();
    float ymin = numeric_limits<float>::max();
    float ymax = numeric_limits<float>::min();
    float zmin = numeric_limits<float>::max();
    float zmax = numeric_limits<float>::min();

    for ( unsigned int i=0; i<rectGLForDisplay.size(); ++i ){
        for ( RectangleGLList::const_iterator eachRectIt = rectGLForDisplay[i].begin();
              eachRectIt!=rectGLForDisplay[i].end(); ++eachRectIt )
        {
            if ( (*eachRectIt).x1 < xmin ){
                xmin = (*eachRectIt).x1;
            }
            if ( (*eachRectIt).x2 > xmax ){
                xmax = (*eachRectIt).x2;
            }
            if ( (*eachRectIt).y1 < ymin ){
                ymin = (*eachRectIt).y1;
            }
            if ( (*eachRectIt).y2 > ymax ){
                ymax = (*eachRectIt).y2;
            }
            if ( (*eachRectIt).z1 < zmin ){
                zmin = (*eachRectIt).z1;
            }
            if ( (*eachRectIt).z2 > zmax ){
                zmax = (*eachRectIt).z2;
            }
        }
    }

    axis[XL] = xmin;
    axis[XU] = xmax;
    axis[YL] = ymin;
    axis[YU] = ymax;
    axis[ZL] = zmin;
    axis[ZU] = zmax;

    xc = (axis[XU]+axis[XL])/2;
    yc = (axis[YU]+axis[YL])/2;
    zc = (axis[ZU]+axis[ZL])/2;
    xr = (axis[XU]-axis[XL])/2;
    yr = (axis[YU]-axis[YL])/2;
    zr = (axis[ZU]-axis[ZL])/2;

    maxr = max(xr, max(yr, zr));
}



void PanelRenderer::initView(){
	translationX = xc;//0.0f;
	translationY = yc;//0.0f;
	translationZ = zc;

    scale = 0.9f;

	rotationX = 0.0f;
	rotationY = 0.0f;
	rotationZ = 0.0f;

    resizeGL(this->width(), this->height());
    updateGL();
}

void PanelRenderer::initializeGL(){
    glEnable(GL_DEPTH_TEST);

#ifdef WHITE_BACKGROUND
    glClearColor(1.0f, 1.0f, 1.0f, 1.0f); //Set the background color
#endif
#ifndef WHITE_BACKGROUND
    glClearColor(0.08f, 0.11f, 0.18f, 1.0f); //Set the background color
#endif
    glEnable(GL_LIGHTING);	//Enable lighting
    glEnable(GL_LIGHT0);	//Enable light #0
    glEnable(GL_LIGHT1);	//Enable light #1
    glEnable(GL_NORMALIZE); //Have OpenGL automatically normalize our normals
    glShadeModel(GL_FLAT);//Enable smooth shading
}

/*_____________________________________________________________________________
 *
 * PanelRenderer:resizeGL
 * : to maintain fixed aspect ratio when the widget size is modified
 *_____________________________________________________________________________
 */
void PanelRenderer::resizeGL(int width, int height){

    float hwRatio = (float)height/width;

    glViewport(0,0,width,height);
    glMatrixMode(GL_PROJECTION);

    glLoadIdentity();

    float zDepth = 3.0;

    if ( hwRatio > yr/xr ){
        glOrtho(-xr,+xr,-xr*hwRatio,+xr*hwRatio,
                (-xr*zDepth),(+xr*zDepth));
    }else{
        glOrtho(-yr/hwRatio,+yr/hwRatio,-yr,+yr,
                (-yr*zDepth),(+yr*zDepth));
    }

    glMatrixMode(GL_MODELVIEW);
}


void PanelRenderer::paintGL()
{
    //* clean up screens
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    //* process rotation centered at display center
    glRotatef(rotationX, 1.0, 0.0, 0.0);
    glRotatef(rotationY, 0.0, 1.0, 0.0);
    glRotatef(rotationZ, 0.0, 0.0, 1.0);

    //* scale objects
    glScalef(scale, scale, scale);
    glTranslatef(-translationX, -translationY, -translationZ);

    //* set up lights
    GLfloat ambientColor[] = {0.7f, 0.7f, 0.7f, 1.0f}; //Color (0.2, 0.2, 0.2)
    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambientColor);
    //* light0
    GLfloat lightColor0[] = {1.0f, 1.0f, 1.0f, 1.0f};
    GLfloat lightPos0[] = {5.0f, 5.0f, 10.0f, 0.0f};
    glLightfv(GL_LIGHT0, GL_SPECULAR, lightColor0);
    glLightfv(GL_LIGHT0, GL_POSITION, lightPos0);
    //* light1
    GLfloat lightColor1[] = {0.5f, 0.5f, 0.5f, 1.0f};
    GLfloat lightPos1[] = {-5.0f, -5.0f, -10.0f, 0.0f};
    glLightfv(GL_LIGHT1, GL_DIFFUSE, lightColor1);
    glLightfv(GL_LIGHT1, GL_POSITION, lightPos1);

    const float shift = 10e-9;
    const float shift2  = 20e-9;


    if ( isLoaded ){
        for ( unsigned int i=0; i<rectGLForDisplay.size();++i ){
            glMaterialfv(GL_FRONT, GL_DIFFUSE, getColor(renderPropertyVec[i].colorIndex));

            glBegin(GL_QUADS);
            for (RectangleGLList::iterator rectIt= rectGLForDisplay[i].begin();
                 rectIt!= rectGLForDisplay[i].end(); ++rectIt)
            {
                glNormal3f(rectIt->xn, rectIt->yn, rectIt->zn);

                if (rectIt->xn != 0){
                    if (rectIt->shapeNormalDistance!=0 && rectIt->xn > 0){
                        glVertex3f( rectIt->x1+shift, rectIt->y1, rectIt->z1 );
                        glVertex3f( rectIt->x1+shift, rectIt->y2, rectIt->z1 );
                        glVertex3f( rectIt->x1+shift, rectIt->y2, rectIt->z2 );
                        glVertex3f( rectIt->x1+shift, rectIt->y1, rectIt->z2 );
                    }
                    else if (rectIt->shapeNormalDistance!=0 && rectIt->xn < 0){
                        glVertex3f( rectIt->x1-shift, rectIt->y1, rectIt->z1 );
                        glVertex3f( rectIt->x1-shift, rectIt->y2, rectIt->z1 );
                        glVertex3f( rectIt->x1-shift, rectIt->y2, rectIt->z2 );
                        glVertex3f( rectIt->x1-shift, rectIt->y1, rectIt->z2 );
                    }
                    else if (rectIt->shapeShift!=0 && rectIt->xn > 0){
                        glVertex3f( rectIt->x1+shift2, rectIt->y1, rectIt->z1 );
                        glVertex3f( rectIt->x1+shift2, rectIt->y2, rectIt->z1 );
                        glVertex3f( rectIt->x1+shift2, rectIt->y2, rectIt->z2 );
                        glVertex3f( rectIt->x1+shift2, rectIt->y1, rectIt->z2 );
                    }
                    else if (rectIt->shapeShift!=0 && rectIt->xn < 0){
                        glVertex3f( rectIt->x1-shift2, rectIt->y1, rectIt->z1 );
                        glVertex3f( rectIt->x1-shift2, rectIt->y2, rectIt->z1 );
                        glVertex3f( rectIt->x1-shift2, rectIt->y2, rectIt->z2 );
                        glVertex3f( rectIt->x1-shift2, rectIt->y1, rectIt->z2 );
                    }
                    glVertex3f( rectIt->x1, rectIt->y1, rectIt->z1 );
                    glVertex3f( rectIt->x1, rectIt->y2, rectIt->z1 );
                    glVertex3f( rectIt->x1, rectIt->y2, rectIt->z2 );
                    glVertex3f( rectIt->x1, rectIt->y1, rectIt->z2 );
                }
                else if (rectIt->yn != 0){
                    if (rectIt->shapeNormalDistance!=0 && rectIt->yn > 0){
                        glVertex3f( rectIt->x1, rectIt->y1+shift, rectIt->z1 );
                        glVertex3f( rectIt->x2, rectIt->y1+shift, rectIt->z1 );
                        glVertex3f( rectIt->x2, rectIt->y1+shift, rectIt->z2 );
                        glVertex3f( rectIt->x1, rectIt->y1+shift, rectIt->z2 );
                    }
                    else if (rectIt->shapeNormalDistance!=0 && rectIt->yn < 0){
                        glVertex3f( rectIt->x1, rectIt->y1-shift, rectIt->z1 );
                        glVertex3f( rectIt->x2, rectIt->y1-shift, rectIt->z1 );
                        glVertex3f( rectIt->x2, rectIt->y1-shift, rectIt->z2 );
                        glVertex3f( rectIt->x1, rectIt->y1-shift, rectIt->z2 );
                    }
                    else if (rectIt->shapeShift!=0 && rectIt->yn > 0){
                        glVertex3f( rectIt->x1, rectIt->y1+shift2, rectIt->z1 );
                        glVertex3f( rectIt->x2, rectIt->y1+shift2, rectIt->z1 );
                        glVertex3f( rectIt->x2, rectIt->y1+shift2, rectIt->z2 );
                        glVertex3f( rectIt->x1, rectIt->y1+shift2, rectIt->z2 );
                    }
                    else if (rectIt->shapeShift!=0 && rectIt->yn < 0){
                        glVertex3f( rectIt->x1, rectIt->y1-shift2, rectIt->z1 );
                        glVertex3f( rectIt->x2, rectIt->y1-shift2, rectIt->z1 );
                        glVertex3f( rectIt->x2, rectIt->y1-shift2, rectIt->z2 );
                        glVertex3f( rectIt->x1, rectIt->y1-shift2, rectIt->z2 );
                    }
                    glVertex3f( rectIt->x1, rectIt->y1, rectIt->z1 );
                    glVertex3f( rectIt->x2, rectIt->y1, rectIt->z1 );
                    glVertex3f( rectIt->x2, rectIt->y1, rectIt->z2 );
                    glVertex3f( rectIt->x1, rectIt->y1, rectIt->z2 );
                }
                else if (rectIt->zn != 0){
                    //* include both Z and FLAT cases
                    if (rectIt->shapeNormalDistance!=0 && rectIt->zn > 0){
                        glVertex3f( rectIt->x1, rectIt->y1, rectIt->z1+shift );
                        glVertex3f( rectIt->x2, rectIt->y1, rectIt->z1+shift );
                        glVertex3f( rectIt->x2, rectIt->y2, rectIt->z1+shift );
                        glVertex3f( rectIt->x1, rectIt->y2, rectIt->z1+shift );
                    }
                    else if (rectIt->shapeNormalDistance!=0 && rectIt->zn < 0){
                        glVertex3f( rectIt->x1, rectIt->y1, rectIt->z1-shift );
                        glVertex3f( rectIt->x2, rectIt->y1, rectIt->z1-shift );
                        glVertex3f( rectIt->x2, rectIt->y2, rectIt->z1-shift );
                        glVertex3f( rectIt->x1, rectIt->y2, rectIt->z1-shift );
                    }
                    else if (rectIt->shapeShift!=0 && rectIt->zn > 0){
                        glVertex3f( rectIt->x1, rectIt->y1, rectIt->z1+shift2 );
                        glVertex3f( rectIt->x2, rectIt->y1, rectIt->z1+shift2 );
                        glVertex3f( rectIt->x2, rectIt->y2, rectIt->z1+shift2 );
                        glVertex3f( rectIt->x1, rectIt->y2, rectIt->z1+shift2 );
                    }
                    else if (rectIt->shapeShift!=0 && rectIt->zn < 0){
                        glVertex3f( rectIt->x1, rectIt->y1, rectIt->z1-shift2 );
                        glVertex3f( rectIt->x2, rectIt->y1, rectIt->z1-shift2 );
                        glVertex3f( rectIt->x2, rectIt->y2, rectIt->z1-shift2 );
                        glVertex3f( rectIt->x1, rectIt->y2, rectIt->z1-shift2 );
                    }
                    glVertex3f( rectIt->x1, rectIt->y1, rectIt->z1 );
                    glVertex3f( rectIt->x2, rectIt->y1, rectIt->z1 );
                    glVertex3f( rectIt->x2, rectIt->y2, rectIt->z1 );
                    glVertex3f( rectIt->x1, rectIt->y2, rectIt->z1 );
                }
                else{
                    cerr << "ERROR: rect has zero normal" << endl;
                }
            }
            glEnd();
        }


        if ( outline == true ){
            GLfloat blackColor[] = { 0.0, 0.0, 0.0, 1.0 };
            glMaterialfv(GL_FRONT, GL_DIFFUSE, blackColor);
            glLineWidth(outlineWidth);

            for ( unsigned int i=0; i<rectGLForDisplay.size();++i ){
                for (RectangleGLList::iterator rectIt= rectGLForDisplay[i].begin();
                     rectIt!= rectGLForDisplay[i].end(); ++rectIt)
                {
                    glBegin(GL_LINE_LOOP);
                    glNormal3f(rectIt->xn, rectIt->yn, rectIt->zn);

                    if (rectIt->xn != 0){
                        if (rectIt->shapeNormalDistance!=0 && rectIt->xn > 0){
                            glVertex3f( rectIt->x1+shift, rectIt->y1, rectIt->z1 );
                            glVertex3f( rectIt->x1+shift, rectIt->y2, rectIt->z1 );
                            glVertex3f( rectIt->x1+shift, rectIt->y2, rectIt->z2 );
                            glVertex3f( rectIt->x1+shift, rectIt->y1, rectIt->z2 );
                        }
                        else if (rectIt->shapeNormalDistance!=0 && rectIt->xn < 0){
                            glVertex3f( rectIt->x1-shift, rectIt->y1, rectIt->z1 );
                            glVertex3f( rectIt->x1-shift, rectIt->y2, rectIt->z1 );
                            glVertex3f( rectIt->x1-shift, rectIt->y2, rectIt->z2 );
                            glVertex3f( rectIt->x1-shift, rectIt->y1, rectIt->z2 );
                        }
                        else if (rectIt->shapeShift!=0 && rectIt->xn > 0){
                            glVertex3f( rectIt->x1+shift2, rectIt->y1, rectIt->z1 );
                            glVertex3f( rectIt->x1+shift2, rectIt->y2, rectIt->z1 );
                            glVertex3f( rectIt->x1+shift2, rectIt->y2, rectIt->z2 );
                            glVertex3f( rectIt->x1+shift2, rectIt->y1, rectIt->z2 );
                        }
                        else if (rectIt->shapeShift!=0 && rectIt->xn < 0){
                            glVertex3f( rectIt->x1-shift2, rectIt->y1, rectIt->z1 );
                            glVertex3f( rectIt->x1-shift2, rectIt->y2, rectIt->z1 );
                            glVertex3f( rectIt->x1-shift2, rectIt->y2, rectIt->z2 );
                            glVertex3f( rectIt->x1-shift2, rectIt->y1, rectIt->z2 );
                        }
                        else{
                        glVertex3f( rectIt->x1, rectIt->y1, rectIt->z1 );
                        glVertex3f( rectIt->x1, rectIt->y2, rectIt->z1 );
                        glVertex3f( rectIt->x1, rectIt->y2, rectIt->z2 );
                        glVertex3f( rectIt->x1, rectIt->y1, rectIt->z2 );
                        }
                    }
                    else if (rectIt->yn != 0){
                        if (rectIt->shapeNormalDistance!=0 && rectIt->yn > 0){
                            glVertex3f( rectIt->x1, rectIt->y1+shift, rectIt->z1 );
                            glVertex3f( rectIt->x2, rectIt->y1+shift, rectIt->z1 );
                            glVertex3f( rectIt->x2, rectIt->y1+shift, rectIt->z2 );
                            glVertex3f( rectIt->x1, rectIt->y1+shift, rectIt->z2 );
                        }
                        else if (rectIt->shapeNormalDistance!=0 && rectIt->yn < 0){
                            glVertex3f( rectIt->x1, rectIt->y1-shift, rectIt->z1 );
                            glVertex3f( rectIt->x2, rectIt->y1-shift, rectIt->z1 );
                            glVertex3f( rectIt->x2, rectIt->y1-shift, rectIt->z2 );
                            glVertex3f( rectIt->x1, rectIt->y1-shift, rectIt->z2 );
                        }
                        else if (rectIt->shapeShift!=0 && rectIt->yn > 0){
                            glVertex3f( rectIt->x1, rectIt->y1+shift2, rectIt->z1 );
                            glVertex3f( rectIt->x2, rectIt->y1+shift2, rectIt->z1 );
                            glVertex3f( rectIt->x2, rectIt->y1+shift2, rectIt->z2 );
                            glVertex3f( rectIt->x1, rectIt->y1+shift2, rectIt->z2 );
                        }
                        else if (rectIt->shapeShift!=0 && rectIt->yn < 0){
                            glVertex3f( rectIt->x1, rectIt->y1-shift2, rectIt->z1 );
                            glVertex3f( rectIt->x2, rectIt->y1-shift2, rectIt->z1 );
                            glVertex3f( rectIt->x2, rectIt->y1-shift2, rectIt->z2 );
                            glVertex3f( rectIt->x1, rectIt->y1-shift2, rectIt->z2 );
                        }
                        else{
                        glVertex3f( rectIt->x1, rectIt->y1, rectIt->z1 );
                        glVertex3f( rectIt->x2, rectIt->y1, rectIt->z1 );
                        glVertex3f( rectIt->x2, rectIt->y1, rectIt->z2 );
                        glVertex3f( rectIt->x1, rectIt->y1, rectIt->z2 );
                        }
                    }
                    else if (rectIt->zn != 0){
                        //* include both Z and FLAT cases
                        if (rectIt->shapeNormalDistance!=0 && rectIt->zn > 0){
                            glVertex3f( rectIt->x1, rectIt->y1, rectIt->z1+shift );
                            glVertex3f( rectIt->x2, rectIt->y1, rectIt->z1+shift );
                            glVertex3f( rectIt->x2, rectIt->y2, rectIt->z1+shift );
                            glVertex3f( rectIt->x1, rectIt->y2, rectIt->z1+shift );
                        }
                        else if (rectIt->shapeNormalDistance!=0 && rectIt->zn < 0){
                            glVertex3f( rectIt->x1, rectIt->y1, rectIt->z1-shift );
                            glVertex3f( rectIt->x2, rectIt->y1, rectIt->z1-shift );
                            glVertex3f( rectIt->x2, rectIt->y2, rectIt->z1-shift );
                            glVertex3f( rectIt->x1, rectIt->y2, rectIt->z1-shift );
                        }
                        else if (rectIt->shapeShift!=0 && rectIt->zn > 0){
                            glVertex3f( rectIt->x1, rectIt->y1, rectIt->z1+shift2 );
                            glVertex3f( rectIt->x2, rectIt->y1, rectIt->z1+shift2 );
                            glVertex3f( rectIt->x2, rectIt->y2, rectIt->z1+shift2 );
                            glVertex3f( rectIt->x1, rectIt->y2, rectIt->z1+shift2 );
                        }
                        else if (rectIt->shapeShift!=0 && rectIt->zn < 0){
                            glVertex3f( rectIt->x1, rectIt->y1, rectIt->z1-shift2 );
                            glVertex3f( rectIt->x2, rectIt->y1, rectIt->z1-shift2 );
                            glVertex3f( rectIt->x2, rectIt->y2, rectIt->z1-shift2 );
                            glVertex3f( rectIt->x1, rectIt->y2, rectIt->z1-shift2 );
                        }
                        else{
                        glVertex3f( rectIt->x1, rectIt->y1, rectIt->z1 );
                        glVertex3f( rectIt->x2, rectIt->y1, rectIt->z1 );
                        glVertex3f( rectIt->x2, rectIt->y2, rectIt->z1 );
                        glVertex3f( rectIt->x1, rectIt->y2, rectIt->z1 );
                        }
                    }
                    else{
                        cerr << "ERROR: rect has zero normal" << endl;
                    }
                    glEnd();
                }
            }

        }

    }

}

void PanelRenderer::conductorFPList2RectGLForDisplay()
{
    int nLayer = condFPListPtr->front().layer.size();
    rectGLForDisplay.clear();

    if ( colorScheme == BYLAYER )
    {//*group rects by layer
        rectGLForDisplay.resize( nLayer );
        for ( list<ConductorFP>::const_iterator eachCondIt = condFPListPtr->begin();
              eachCondIt != condFPListPtr->end(); ++eachCondIt)
        {
            for ( int i=0; i<nLayer; ++i )
            {
                for ( int j=0; j<nDir; ++j )
                {

                    rectGLForDisplay[i].insert(rectGLForDisplay[i].end(),
                                              eachCondIt->layer[i][j].begin(),
                                              eachCondIt->layer[i][j].end());
                }
            }
        }
    }
    else if ( colorScheme == BYCONDUCTOR )
    {//*group rects by conductor
        rectGLForDisplay.resize( condFPListPtr->size() );
        list<ConductorFP>::const_iterator eachCondIt;
        unsigned int i;

        for ( i=0, eachCondIt = condFPListPtr->begin();
              i<rectGLForDisplay.size(); ++i, ++eachCondIt )
        {
            for ( int j=0; j<nLayer; ++j )
            {
                for ( int k=0; k<nDir; ++k )
                {

                    rectGLForDisplay[i].insert(rectGLForDisplay[i].end(),
                                              eachCondIt->layer[j][k].begin(),
                                              eachCondIt->layer[j][k].end());

                }
            }
        }
    }
    setupGLRenderProperty();
}



void PanelRenderer::mousePressEvent(QMouseEvent *event)
{
	lastPos = event->pos();

	enum MOUSE_MODE { ROTATE, TRANSLATE, ZOOM };
	enum MOUSE_AXIS { AXIS_XY, AXIS_Z };

	int	lx = lastPos.x();
	int	ly = lastPos.y();
	QPoint center = rect().center();

	int xCenter = center.x();
	int yCenter = center.y();
	int maxRadius2 = (width()>height())? height()/2 : width()/2;
	int lastRadius2 = static_cast<int>(std::sqrt((lx-xCenter)*(lx-xCenter) + (ly-yCenter)*(ly-yCenter)));

	mouseAxis = ( lastRadius2 > maxRadius2 * 0.9 ) ? AXIS_Z : AXIS_XY;


	switch( event->modifiers() ){
	case Qt::CTRL: // zoom
		mouseMode = ZOOM;
		break;
	case Qt::SHIFT: // translate
		mouseMode = TRANSLATE;
		break;
	default:
		mouseMode = ROTATE;
	}
}


void PanelRenderer::mouseMoveEvent(QMouseEvent *event)
{
	enum MOUSE_MODE { ROTATE, TRANSLATE, ZOOM };
	enum MOUSE_AXIS { AXIS_XY, AXIS_Z };

	int	lx = lastPos.x();
	int	ly = lastPos.y();
	int	ex = event->x();
	int ey = event->y();


	GLfloat dx = GLfloat(ex - lx) / width();
	GLfloat dy = GLfloat(ey - ly) / height();

	switch ( mouseMode ){
	case ROTATE:
		if ( mouseAxis == AXIS_Z ){
			int xCenter = rect().center().x();
			int yCenter = rect().center().y();
			if ( lx > xCenter ){
				rotationZ -= 180 * dy;
			}else{
				rotationZ += 180 * dy;
			}
			if ( ly > yCenter ){
				rotationZ += 180 * dx;
			}else{
				rotationZ -= 180 * dx;
			}
		}else{
			rotationX += 180 * dy;
			rotationY += 180 * dx;
		}
		break;
	case TRANSLATE:
		if ( mouseAxis == AXIS_Z ){
            translationZ += dy*scale*maxr;
		}else{
            translationX -= dx*scale*maxr;
            translationY += dy*scale*maxr;
		}
		break;
	case ZOOM:
		scale *= 1+dx*2-dy*2;

	}

	updateGL();
	lastPos = event->pos();
}

void PanelRenderer::keyPressEvent(QKeyEvent *event){

	switch(event->key()){
	case Qt::Key_F:
		initView();
		switch( event->modifiers() ){
		case Qt::SHIFT: // xz view
			rotationX = 90;
			rotationY = 180;
			break;
		case Qt::CTRL:  // yz view
			rotationX = -90;
			rotationZ = -90;
		}
		updateGL();
		break;
	case Qt::Key_D:
		initView();
		switch( event->modifiers() ){
		case Qt::SHIFT: // -xz view
			rotationX = -90;
			break;
		case Qt::CTRL:  // -yz view
			rotationX = -90;
			rotationZ =  90;
			break;
		default:
			rotationY = 180;
		}
		updateGL();
        break;
    case Qt::Key_R:
        initView();
        rotationX = -60;
        rotationY = 0;
        rotationZ = 60;
        updateGL();
        break;
	}
}


