#-------------------------------------------------
#
# Project created by QtCreator 2013-06-17T07:13:00
#
#-------------------------------------------------

QT       -= gui

TARGET = Vtrlib
TEMPLATE = lib

DEFINES += VTRLIB_LIBRARY \
           VTRLIB_EXPORTS

DLLDESTDIR = ../../BinQt

INCLUDEPATH += ../../General

SOURCES += ../../Vtrlib/Vtr.cpp

HEADERS += ../../Vtrlib/Vtr.h

unix:!symbian {
    maemo5 {
        target.path = /opt/usr/lib
    } else {
        target.path = /usr/lib
    }
    INSTALLS += target
}
