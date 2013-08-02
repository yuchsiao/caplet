/*
CREATED : Jan 31, 2013
AUTHOR  : Yu-Chung Hsiao
EMAIL   : project.caplet@gmail.com

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

#ifndef PANELRENDERER_H
#define PANELRENDERER_H

#include "geoloader.h"

#include <QGLWidget>
#include <vector>



class QRadialGradient;
class QColor;

class RenderProperty{
public:
    int colorIndex;
};

class PanelRenderer : public QGLWidget
{
    Q_OBJECT

public:
	PanelRenderer( QWidget* parent = 0);
    ~PanelRenderer();
	void keyPressEvent(QKeyEvent * event);
    void loadGLRects( const ConductorFPList *condFPListPtr);
    void initView();


    enum COLORSCHEME { BYLAYER, BYCONDUCTOR };

    void plotOutline(bool outlineFlag);
    void clear();

    //* color
    QColor** color;
    int      nColor;

public slots:
    void setGLColorScheme( int scheme );


protected:
	void initializeGL();
	void resizeGL(int width, int height);
	void paintGL();

	void mousePressEvent(QMouseEvent *event);
	void mouseMoveEvent(QMouseEvent *event);


private:
    void conductorFPList2RectGLForDisplay();
    void setupGLRenderProperty();
    void computeGLDisplayBoundary();

    std::vector<RenderProperty>     renderPropertyVec;

    const ConductorFPList           *condFPListPtr;

    std::vector<RectangleGLList>    rectGLForDisplay;

    int     colorScheme;


    float	(**panelArray)[12];
    float	(**normalArray)[3];
	int		nPanels;
	bool	isLoaded;
	float	axis[6];
	float	xc;
	float	yc;
	float	zc;
	float	xr;
	float	yr;
	float	zr;
	float	maxr;

	QPoint lastPos;

	GLfloat rotationX;
	GLfloat rotationY;
	GLfloat rotationZ;
	GLfloat translationX;
	GLfloat translationY;
    GLfloat translationZ;
	QRadialGradient gradient;
	int		mouseMode;
	int		mouseAxis;
	float	scale;
    float	outlineWidth;

    bool    outline;

    //* color setup
    GLfloat  currentColor[4];
    int      colorIndex;
    inline GLfloat* getColor(int colorIndex){
        currentColor[0] = color[colorIndex]->redF();
        currentColor[1] = color[colorIndex]->greenF();
        currentColor[2] = color[colorIndex]->blueF();
        currentColor[3] = 1.0f;
        return currentColor;
    }

    static const int nDir = 6;
};


#endif // PANELRENDERER_H
