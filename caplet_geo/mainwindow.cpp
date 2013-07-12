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

#include "mainwindow.h"
#include "ui_mainwindow.h"

#include <QtGui>
#include <QDebug>
#include <QErrorMessage>
#include <QFileInfo>

#include "panelrenderer.h"
#include "geoloader.h"

#include <iostream>
#include <fstream>
#include <string>
#include <list>
#include <vector>
#include <limits>

using namespace std;

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    panelRenderer(0),
    geoLoader(0),
    fileName(""),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    setWindowTitle(QString("Caplet Geo"));

    geoLoader     = new GeoLoader();
    panelRenderer = new PanelRenderer();

    setCentralWidget(panelRenderer);

    //* menuBar actions
    colorSchemeGroup = new QActionGroup(this);
    ui->actionColorSchemeLayer->setActionGroup(colorSchemeGroup);
    ui->actionColorSchemeConductor->setActionGroup(colorSchemeGroup);

    solverGroup = new QActionGroup(this);
    ui->actionNone->setActionGroup(solverGroup);
    ui->actionCAPLET->setActionGroup(solverGroup);
    ui->actionFASTCAP->setActionGroup(solverGroup);
    ui->actionStandardBEM->setActionGroup(solverGroup);

    pwcSizeValidator = new QDoubleValidator(this);
    pwcSizeValidator->setBottom(0);

    archLengthValidator = new QDoubleValidator(this);
    archLengthValidator->setBottom(0);

    coreNumValidator = new QIntValidator(this);
    coreNumValidator->setBottom(1);

    initPWCSizeValidator = new QDoubleValidator(this);
    initPWCSizeValidator->setBottom(0);

    alphaValidator   = new QDoubleValidator(this);
    alphaValidator  ->setBottom(1);

    epsilonValidator = new QDoubleValidator(this);
    epsilonValidator->setTop(100);
    epsilonValidator->setBottom(0);

    //* validator
    ui->pwcSizeLineEdit    ->setValidator(pwcSizeValidator);
    ui->archLengthLineEdit ->setValidator(archLengthValidator);
    ui->coreNumLineEdit    ->setValidator(coreNumValidator);
    ui->initPWCSizeLineEdit->setValidator(initPWCSizeValidator);
    ui->alphaLineEdit      ->setValidator(alphaValidator);
    ui->epsilonLineEdit    ->setValidator(epsilonValidator);

    //* log text
    ui->logText->setReadOnly(true);
    ui->logText->ensureCursorVisible();
}


MainWindow::~MainWindow(){
    delete epsilonValidator;
    delete alphaValidator;
    delete initPWCSizeValidator;
    delete coreNumValidator;
    delete archLengthValidator;
    delete pwcSizeValidator;
    delete ui;
    clear();
}

void MainWindow::clear(){
    delete geoLoader;
    delete panelRenderer;
}

void MainWindow::displayGeometry()
{
    if (isLoaded == false) {
        return;
    }
    panelRenderer->loadGLRects(&(geoLoader->getGeometryConductorList(unit)));
}

void MainWindow::displayPWCBasisFunction()
{
    if (isLoaded == false) {
        return;
    }
    float suggestedPanelSize = ui->pwcSizeLineEdit->text().toFloat() * unit;
    const ConductorFPList &condList = geoLoader->getPWCBasisFunction(unit, suggestedPanelSize);
    panelRenderer->loadGLRects(&condList);
}

void MainWindow::displayInstantiableBasisFunction()
{
    if (isLoaded == false) {
        return;
    }
    const float archUnit = 1e-6;
    float archLength= ui->archLengthLineEdit->text().toFloat() * archUnit;

    if (ui->flatCheckBox->isChecked()==false){
        //* no flat and no arch
        archLength = -1;
    }
    else if (ui->archCheckBox->isChecked()==false){
        //* no arch
        archLength = 0;
    }

    const ConductorFPList &condList = geoLoader->getInstantiableBasisFunction(unit, archLength);
    panelRenderer->loadGLRects(&condList);
}

//**
//* loadFile
//* - geoLoader: load geometry
//* - panelRenderer: load rects
void MainWindow::loadFile(float unit)
{
    //* set unit (int to float ratio)
    this->unit = unit;

    //* initial empty file name guard
    if (fileName.isEmpty()){
        return;
    }

    try{
        geoLoader->loadGeo( fileName.toUtf8().data() );
        isLoaded = true;
    }
    catch (FileNotFoundError &e){
        QMessageBox::critical(this, "File Not Found", fileName + " is not found.");
        return;
    }
    catch (GeometryNotManhattanError &e){
        QMessageBox::critical(this, "Invliad Geometries", e.what());
        fileName = "";
        return;
    }
    on_actionNone_triggered();
}

//**
//* keyPressEvent
//* - Receive top view and two side views
//* - ESC key is received by other mechanism
void MainWindow::keyPressEvent(QKeyEvent *e)
{
    switch ( e->key() ){
    case Qt::Key_F:
    case Qt::Key_D:
    case Qt::Key_R:
        panelRenderer->keyPressEvent(e);
    default:
        QWidget::keyPressEvent(e);
    }
}


//****
//*
//* button and menu action functions
//*
//****
void MainWindow::byLayerClicked(){
    panelRenderer->setGLColorScheme(PanelRenderer::BYLAYER);
}

void MainWindow::byConductorClicked(){
    panelRenderer->setGLColorScheme(PanelRenderer::BYCONDUCTOR);
}

void MainWindow::on_actionOpenFile_triggered(){
    QString inputFileName = QFileDialog::getOpenFileName(
            this,
            tr("Open File"),
            QDir::currentPath(),
                tr("Caplet Geometry files [*.geo] (*.geo);;Fastcap files [*.qui] (*.qui);;All files (*.*)") );

    if ( inputFileName.isEmpty()==true ){
        return;
    }

    //* check file exsistence
    ifstream fin(inputFileName.toUtf8().data());
    if ( !fin.is_open() ){
        QMessageBox::critical(this, "File Not Found", inputFileName + " is not found.");
        fin.close();
        return;
    }
    fin.close();

    //* reset
    ui->pwcSizeLineEdit->setEnabled(true);

    QString suffix = QFileInfo(inputFileName).suffix();
    if (suffix.compare("geo")==0){
        if ( fileName.compare(inputFileName)!=0 ){

            //* inputFileName exists and can be opened
            fileName        = inputFileName;
            QFileInfo fileInfo(inputFileName);
            fileBaseName    = fileInfo.completeBaseName();
            fileSuffix      = fileInfo.suffix();
            canonicalPath   = fileInfo.canonicalPath();

            loadFile();
            panelRenderer->initView();
        }
    }
    //* not quite working
    //* only for display purpose
    else if (suffix.compare("qui")==0){
        fileName = inputFileName;
        geoLoader->loadQui(fileName.toUtf8().data());
        isLoaded = true;
        panelRenderer->loadGLRects(&(geoLoader->getPWCBasisFunction()));
        panelRenderer->initView();
        ui->pwcSizeLineEdit->setDisabled(true);
    }
    else if (suffix.compare("gds")==0){

    }
}

void MainWindow::on_actionExit_triggered(){
    close();
}

void MainWindow::on_actionColorSchemeLayer_triggered(){
    if (ui->colorSchemeLayerRadio->isChecked()==false){
        ui->colorSchemeLayerRadio->click();
    }
    byLayerClicked();
}

void MainWindow::on_actionColorSchemeConductor_triggered(){
    if (ui->colorSchemeConductorRadio->isChecked()==false){
        ui->colorSchemeConductorRadio->click();
    }
    byConductorClicked();
}

void MainWindow::on_colorSchemeLayerRadio_clicked(){
    if (ui->actionColorSchemeLayer->isChecked()==false){
        ui->actionColorSchemeLayer->trigger();
    }
}

void MainWindow::on_colorSchemeConductorRadio_clicked(){
    if (ui->actionColorSchemeConductor->isChecked()==false){
        ui->actionColorSchemeConductor->trigger();
    }
}

void MainWindow::on_actionNone_triggered()
{
    if ( ui->solverNoneRadio->isChecked()==false ){
        ui->solverNoneRadio->click();
    }
    displayGeometry();
}

void MainWindow::on_actionCAPLET_triggered(){
    if ( ui->solverCapletRadio->isChecked()==false ){
        ui->solverCapletRadio->click();
    }
    displayInstantiableBasisFunction();
}

void MainWindow::on_actionFASTCAP_triggered(){
    if ( ui->solverFastcapRadio->isChecked()==false ){
        ui->solverFastcapRadio->click();
    }
    displayPWCBasisFunction();
}

void MainWindow::on_actionStandardBEM_triggered(){
    if ( ui->solverStandardRadio->isChecked()==false ){
        ui->solverStandardRadio->click();
    }
    displayPWCBasisFunction();
}

void MainWindow::on_actionOutline_triggered(bool checked){
    if ( ui->outlineBox->isChecked()!=ui->actionOutline->isChecked()){
        ui->outlineBox->click();
    }
    panelRenderer->plotOutline(checked);
}


void MainWindow::on_solverNoneRadio_clicked()
{
    if ( ui->actionNone->isChecked()==false){
        ui->actionNone->trigger();
    }
}

void MainWindow::on_solverCapletRadio_clicked(){
    if ( ui->actionCAPLET->isChecked()==false ){
        ui->actionCAPLET->trigger();
    }
}

void MainWindow::on_solverFastcapRadio_clicked(){
    if ( ui->actionFASTCAP->isChecked()==false ){
        ui->actionFASTCAP->trigger();
    }
}

void MainWindow::on_solverStandardRadio_clicked(){
    if ( ui->actionStandardBEM->isChecked()==false ){
        ui->actionStandardBEM->trigger();
    }
}

void MainWindow::on_outlineBox_clicked(bool checked){
    if ( checked != ui->actionOutline->isChecked()){
        ui->actionOutline->trigger();
    }
}


void MainWindow::on_pwcSizeLineEdit_returnPressed()
{
    if (ui->solverStandardRadio->isChecked() == false &&
        ui->solverFastcapRadio->isChecked() == false){

        on_actionFASTCAP_triggered();
    }
    else if (ui->solverStandardRadio->isChecked()){
        on_actionStandardBEM_triggered();
    }
    else if (ui->solverFastcapRadio->isChecked()){
        on_actionFASTCAP_triggered();
    }
}

void MainWindow::on_extractionButton_clicked()
{
    const float percent = 0.01;

    QString pathFileBaseName = canonicalPath+"/"+fileBaseName;

    if (ui->solverFastcapRadio->isChecked()){
        pathFileBaseName += "_"+ui->pwcSizeLineEdit->text();
        geoLoader->runFastcap(pathFileBaseName.toUtf8().data());
        const ExtractionInfo result = geoLoader->getLastResult();
        QString resultText = result.toString().c_str();
        float error = result.compare(geoLoader->getReferenceResult());
        vector<float> errorVector = result.compareDiagonal(geoLoader->getReferenceResult());
        QString errorVectorString;
        for (unsigned i=0; i<errorVector.size(); ++i){
            errorVectorString += QString::number(errorVector[i]/percent) + " ";
        }

        QString msg = QString("# (") + QDateTime().currentDateTime().toString() + QString("): Fastcap extraction is done!");
        logTime(msg);
        log(pathFileBaseName);
        log(resultText);
        log(QString("Diagonal Error(%) = ")+errorVectorString);
        log(QString("Error(%)          = ")+QString::number(error/percent));
    }
    else if (ui->solverCapletRadio->isChecked()){
        int coreNum = ui->coreNumLineEdit->text().toInt();
        pathFileBaseName += "_np" + QString::number(coreNum);
        geoLoader->runCaplet(pathFileBaseName.toUtf8().data(), coreNum);
        const ExtractionInfo result = geoLoader->getLastResult();
        QString resultText = result.toString().c_str();
        float error = result.compare(geoLoader->getReferenceResult());
        vector<float> errorVector = result.compareDiagonal(geoLoader->getReferenceResult());
        QString errorVectorString;
        for (unsigned i=0; i<errorVector.size(); ++i){
            errorVectorString += QString::number(errorVector[i]/percent) + " ";
        }

        QString msg = QString("# (") + QDateTime().currentDateTime().toString() + QString("): Caplet extraction is done!");
        logTime(msg);
        log(pathFileBaseName);
        log(resultText);
        log(QString("Diagonal Error(%) = ")+errorVectorString);
        log(QString("Error(%) = ")+QString::number(error/percent));
    }
    else if (ui->solverStandardRadio->isChecked()){
        pathFileBaseName += "_"+ui->pwcSizeLineEdit->text();
        geoLoader->runCapletQui(pathFileBaseName.toUtf8().data());
        const ExtractionInfo result = geoLoader->getLastResult();
        QString resultText = result.toString().c_str();
        float error = result.compare(geoLoader->getReferenceResult());
        vector<float> errorVector = result.compareDiagonal(geoLoader->getReferenceResult());
        QString errorVectorString;
        for (unsigned i=0; i<errorVector.size(); ++i){
            errorVectorString += QString::number(errorVector[i]/percent) + " ";
        }

        QString msg = QString("# (") + QDateTime().currentDateTime().toString() + QString("): Caplet extraction is done!");
        logTime(msg);
        log(pathFileBaseName);
        log(resultText);
        log(QString("Diagonal Error(%) = ")+errorVectorString);
        log(QString("Error(%) = ")+QString::number(error/percent));
    }
}

void MainWindow::on_referenceButton_clicked()
{
    if (isLoaded == false){
        return;
    }
    geoLoader->clearResult();
    const float percent = 0.01;
    const string fastcapFlag = "-o4";

    float alpha   = ui->alphaLineEdit->text().toFloat();
    float epsilon = ui->epsilonLineEdit->text().toFloat() * percent;
    float pwcSize = ui->initPWCSizeLineEdit->text().toFloat() * unit;

    geoLoader->getPWCBasisFunction(unit, pwcSize);
    QString pathFileBaseName = canonicalPath+"/"+fileBaseName+"_"+QString::number(pwcSize/unit);
    ExtractionInfo *thisResult = &(geoLoader->runFastcap(pathFileBaseName.toUtf8().data(), fastcapFlag) );
    ExtractionInfo *prevResult;
    int thisNPanel = thisResult->nBasisFunction;
    int prevNPanel = 0;
    float err = 1;
    logTime("Start to compute the reference capacitance matrix...");
    int counter = 1;
    do{
        ++counter;
        prevResult  = thisResult;
        prevNPanel  = thisNPanel;
        pwcSize    /= alpha;
        geoLoader->getPWCBasisFunction(unit, pwcSize);
        QString pathFileBaseName = canonicalPath+"/"+fileBaseName+"_"+QString::number(pwcSize);
        thisResult  = &(geoLoader->runFastcap(pathFileBaseName.toUtf8().data(), fastcapFlag) );
        thisNPanel  = thisResult->nBasisFunction;
        try{
            err         = thisResult->compare(*prevResult);
        }
        catch(length_error &e){
            log(e.what());
            log("Dump two matrices...");
            log("Matrix this:");
            stringstream ssThis;
            thisResult->printMatrix(ssThis);
            log(QString(ssThis.str().c_str()));
            log("Matrix reference:");
            stringstream ssPrev;
            prevResult->printMatrix(ssPrev);
            log(QString(ssPrev.str().c_str()));
        }

        QString msg = QString::number(counter) + QString( ". (PWCSize, Error(%) ):")
                + QString::number(pwcSize/unit) + ": " + QString::number(err/percent);
        log(msg);
    }while(err>epsilon || thisNPanel==prevNPanel);

    geoLoader->storeLastAsReference();
    list<ExtractionInfo> resultList = geoLoader->compareAllAgainstReference();
    stringstream ssHeader;
    resultList.front().printBasicHeader(ssHeader);
    log(ssHeader.str().c_str());
    counter = 1;
    for ( list<ExtractionInfo>::const_iterator each = resultList.begin();
          each != --resultList.end(); ++each, ++counter){
        stringstream ssTableLine;
        each->printBasicLine(ssTableLine);
        log(QString::number(counter) + QString(", ") + ssTableLine.str().c_str());
    }

    logTime("Done comuting the reference capacitance matrix.");
}

void MainWindow::logTime(const QString &msg, bool flagBlue)
{
    if (flagBlue==true){
        ui->logText->appendHtml(QString("<font color=blue>")+msg+"</font>");
    }
    else{
        ui->logText->appendHtml(QString("<font color=black>")+msg+"</font>");
    }
}

void MainWindow::log(const QString &msg)
{
    ui->logText->appendPlainText(msg);
}


void MainWindow::on_actionPalette_triggered()
{
    colorpalette = new ColorPalette(panelRenderer->color, panelRenderer->nColor, this);
    colorpalette->setVisible(true);
}


void MainWindow::on_flatCheckBox_clicked(bool checked)
{
    if (checked==false) {
        ui->archCheckBox->setEnabled(false);
        ui->archLengthLineEdit->setEnabled(false);
    }
    else {
        ui->archCheckBox->setEnabled(true);
        ui->archLengthLineEdit->setEnabled(true);
    }
    on_actionCAPLET_triggered();
}

void MainWindow::on_archCheckBox_clicked(bool checked)
{
    if (checked==false) {
        ui->archLengthLineEdit->setEnabled(false);
    }
    else {
        ui->archLengthLineEdit->setEnabled(true);
    }
    on_actionCAPLET_triggered();
}

void MainWindow::on_archLengthLineEdit_returnPressed()
{
    on_actionCAPLET_triggered();
}

void MainWindow::on_actionAbout_triggered()
{
    QMessageBox msgBox(this);
    msgBox.setText("<br>"
                   "<b>CAPLET 1.0</b> <br>"
                   "2013.02 <br>"
                   "<br>"
                   "Yu-Chung Hsiao <br>"
                   "yuchsiao@mit.edu"
                   "<br>"
                   "GNU LGPL v3.0");
    QPixmap img("icons/icon2.png");
    img = img.scaledToWidth(120);
    msgBox.setIconPixmap(img);
    msgBox.setWindowTitle("About CAPLET");
    msgBox.exec();
}
