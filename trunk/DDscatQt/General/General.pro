#-------------------------------------------------
#
# Project created by QtCreator 2013-01-14T18:34:09
#
#-------------------------------------------------

QT       -= gui

TARGET = General
TEMPLATE = lib

DLLDESTDIR = ../../BinQt

DEFINES += GENERAL_LIBRARY \
           GENERAL_EXPORTS

INCLUDEPATH += ../../Xml/include

LIBS += ../../BinQt/libxml2.dll

SOURCES += ../../General/Definitions.cpp \
           ../../General/Enumerator.cpp \
           ../../General/Functions.cpp \
           ../../General/FourArray.cpp \
           ../../General/General.cpp \
           ../../General/Scaner.cpp \
           ../../General/XmlHelper.cpp

HEADERS += ../../General/Definitions.h \
           ../../General/Enumerator.h \
           ../../General/Functions.h \
           ../../General/FourArray.h \
           ../../General/General.h \
           ../../General/Scaner.h \
           ../../General/XmlHelper.h \
           ../../General/ArrayF.h \
           ../../General/Math3.h \
           ../../General/MathDouble.h \
           ../../General/MathFloat.h \
           ../../General/Vect3.h \
           ../../General/Vect6.h

unix:!symbian {
    maemo5 {
        target.path = /opt/usr/lib
    } else {
        target.path = /usr/lib
    }
    INSTALLS += target
}
