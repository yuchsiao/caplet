/*
CREATED : Jan 31, 2013
MODIFIED: Feb 15, 2013
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

#ifndef COLORPALETTE_H
#define COLORPALETTE_H

#include <QDialog>

class QColor;
class QLabel;
class QVBoxLayout;

namespace Ui {
class ColorPalette;
}

class ColorPalette : public QDialog
{
    Q_OBJECT
    
public:
    explicit ColorPalette(QColor** color, int nColor, QWidget *parent = 0);
    ~ColorPalette();
    
private:
    QLabel** label;
    int      nLabel;
    Ui::ColorPalette *ui;
};

#endif // COLORPALETTE_H
