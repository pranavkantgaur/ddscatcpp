#-------------------------------------------------
#
# Project created by QtCreator 2013-01-14T18:02:31
#
#-------------------------------------------------

QT       -= gui

TARGET = Fourierlib
TEMPLATE = lib

DEFINES += FOURIERLIB_LIBRARY \
           FOURIERLIB_EXPORTS

DLLDESTDIR = ../../BinQt

INCLUDEPATH += ../../General

LIBS += ../../BinQt/General.dll

SOURCES = ../../Fourierlib/AbstractFftEngine.cpp \
          ../../Fourierlib/FftEngineFftw.cpp \
          ../../Fourierlib/FftEngineGpfaft.cpp \
          ../../Fourierlib/FftEngineMkl.cpp \
          ../../Fourierlib/FftEngineCuda.cpp

HEADERS = ../../Fourierlib/Fourierlib.h \
          ../../Fourierlib/AbstractFftEngine.h \
          ../../Fourierlib/FftEngineFftw.h \
          ../../Fourierlib/FftEngineGpfaft.h \
          ../../Fourierlib/FftEngineMkl.h \
          ../../Fourierlib/FftEngineCuda.h

unix:!symbian {
    maemo5 {
        target.path = /opt/usr/lib
    } else {
        target.path = /usr/lib
    }
    INSTALLS += target
}
