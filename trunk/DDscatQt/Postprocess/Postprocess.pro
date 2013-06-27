TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

DESTDIR = ../../BinQt

INCLUDEPATH += ../../General \
               ../../Processlib \
               ../../Vtrlib

LIBS += ../../BinQt/General.dll \
        ../../BinQt/Processlib.dll

SOURCES += ../../Postprocess/Postprocess.cpp


