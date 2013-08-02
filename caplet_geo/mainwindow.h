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

#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include "geoloader.h"
#include "panelrenderer.h"
#include "colorpalette.h"

#include <QMainWindow>
#include <QString>
#include <QActionGroup>
#include <QDoubleValidator>



class PanelRenderer;
class GeoLoader;


namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT
    
public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();
    void loadFile(float unit = 1e-9); // default: 1nm grid for GDS
    void clear();

    void displayGeometry();
    void displayPWCBasisFunction();
    void displayInstantiableBasisFunction();

public slots:
    void byLayerClicked();
    void byConductorClicked();

protected:
    void keyPressEvent(QKeyEvent *event);

private slots:
    //* file
    void on_actionOpenFile_triggered();
    void on_actionExit_triggered();

    //**
    //* The buttons and menu actions are bundled together
    //*   through _triggered() and _clicked().
    //* The functions on_action*_triggered() do self-checking
    //*   and perform the actions, whereas the functions on_*_clicked()
    //*   are the shortcut to connect clicked() signals to slots
    //*   and call on_action*_triggered() for the actual actions.
    //* Better implementation can be possibly done but
    //*   for now this way serves the purpose though tedious.

    //* color scheme
    void on_actionColorSchemeLayer_triggered();
    void on_actionColorSchemeConductor_triggered();
    void on_colorSchemeLayerRadio_clicked();
    void on_colorSchemeConductorRadio_clicked();
    void on_actionOutline_triggered(bool checked);
    void on_outlineBox_clicked(bool checked);

    //* solver type
    void on_actionNone_triggered();
    void on_actionCAPLET_triggered();
    void on_actionFASTCAP_triggered();
    void on_actionStandardBEM_triggered();
    void on_solverNoneRadio_clicked();
    void on_solverCapletRadio_clicked();
    void on_solverFastcapRadio_clicked();
    void on_solverStandardRadio_clicked();

    //* parameter setting event
    void on_pwcSizeLineEdit_returnPressed();
    void on_extractionButton_clicked();
    void on_referenceButton_clicked();
    void on_flatCheckBox_clicked(bool checked);
    void on_archCheckBox_clicked(bool checked);
    void on_archLengthLineEdit_returnPressed();

    //* misc event
    void on_actionPalette_triggered();
    void on_actionAbout_triggered();

private:
    bool                isLoaded;
    float               unit;
    PanelRenderer*      panelRenderer;
    GeoLoader*          geoLoader;
    ColorPalette*       colorpalette;

    QString             fileName;
    QString             fileBaseName;
    QString             fileSuffix;
    QString             canonicalPath;

    QActionGroup*       colorSchemeGroup;
    QActionGroup*       solverGroup;
    QDoubleValidator*   pwcSizeValidator;
    QDoubleValidator*   archLengthValidator;
    QIntValidator*      coreNumValidator;

    QDoubleValidator*   initPWCSizeValidator;
    QDoubleValidator*   alphaValidator;
    QDoubleValidator*   epsilonValidator;

    void logTime(const QString &msg, bool flagBlue=true);
    void log(const QString &msg);

private:
    Ui::MainWindow *ui;
};

#endif // MAINWINDOW_H
