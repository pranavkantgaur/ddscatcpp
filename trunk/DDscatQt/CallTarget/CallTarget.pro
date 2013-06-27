TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

DESTDIR = ../../BinQt

INCLUDEPATH += ../../General \
               ../../Targetlib

LIBS += ../../BinQt/General.dll \
        ../../BinQt/Targetlib.dll

SOURCES += ../../CallTarget/CallTarget.cpp \
           ../../CallTarget/CallTargetMain.cpp
