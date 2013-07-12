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

#include "colorpalette.h"
#include "ui_colorpalette.h"

#include <QColor>
#include <QLabel>
#include <QVBoxLayout>
#include <QDebug>

#include <iostream>

ColorPalette::ColorPalette(QColor** color, int nColor, QWidget *parent) :
    QDialog(parent), nLabel(nColor),
    ui(new Ui::ColorPalette)
{
    ui->setupUi(this);
    int r;
    int g;
    int b;

    label = new QLabel*[nLabel];
    for (int i=0; i<nLabel; ++i){
        label[i] = new QLabel(QString::number(i+1));
        label[i]->setVisible(true);
        ui->formLayout->addRow(label[i]);
        r = color[i]->red();
        g = color[i]->green();
        b = color[i]->blue();
        QString style = "QLabel { background-color : rgb("
                + QString::number(r) + ", "
                + QString::number(g) + ", "
                + QString::number(b) + ")}";
        label[i]->setStyleSheet(style);
        label[i]->setAlignment(Qt::AlignHCenter);
    }
}

ColorPalette::~ColorPalette()
{
    for (int i=0; i<nLabel; ++i){
        delete label[i];
    }
    delete label;

    delete ui;
}
