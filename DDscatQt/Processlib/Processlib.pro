#-------------------------------------------------
#
# Project created by QtCreator 2013-06-17T07:00:43
#
#-------------------------------------------------

QT       -= gui

TARGET = Processlib
TEMPLATE = lib

DEFINES += PROCESSLIB_LIBRARY \
           PROCESSLIB_EXPORTS

DLLDESTDIR = ../../BinQt

INCLUDEPATH += ../../General \
               ../../Vtrlib

LIBS += ../../BinQt/General.dll \
        ../../BinQt/Vtrlib.dll

SOURCES += ../../Processlib/ProcessHelper.cpp \
           ../../Processlib/Linia.cpp \
           ../../Processlib/Processlib.cpp

HEADERS += ../../Processlib/ProcessHelper.h \
           ../../Processlib/Linia.h \
           ../../Processlib/Processlib.h

unix:!symbian {
    maemo5 {
        target.path = /opt/usr/lib
    } else {
        target.path = /usr/lib
    }
    INSTALLS += target
}
