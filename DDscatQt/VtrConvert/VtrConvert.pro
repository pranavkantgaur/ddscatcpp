TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

DESTDIR = ../../BinQt

INCLUDEPATH += ../../Vtrlib \
               ../../General

LIBS += ../../BinQt/Vtrlib.dll

SOURCES += ../../VtrConvert/VtrConvert.cpp

HEADERS += ../../VtrConvert/VtrConvert.h

