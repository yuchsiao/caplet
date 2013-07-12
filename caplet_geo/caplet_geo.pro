#-------------------------------------------------
#
# Project created by QtCreator 2012-01-02T20:09:36
#
#-------------------------------------------------

QT       += core gui opengl

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = caplet_geo
TEMPLATE = app


SOURCES += main.cpp\
        mainwindow.cpp \
    geoloader.cpp \
    panelrenderer.cpp \
    gdsgeometry.cpp \
    colorpalette.cpp

HEADERS  += mainwindow.h \
    geoloader.h \
    panelrenderer.h \
    gdsgeometry.h \
    debug.h \
    colorpalette.h

FORMS    += mainwindow.ui \
    colorpalette.ui

RESOURCES += \
    icons.qrc

