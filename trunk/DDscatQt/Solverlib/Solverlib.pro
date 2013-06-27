#-------------------------------------------------
#
# Project created by QtCreator 2013-01-14T18:01:06
#
#-------------------------------------------------

QT       -= gui

TARGET = Solverlib
TEMPLATE = lib

DEFINES += SOLVERLIB_LIBRARY \
           SOLVERLIB_EXPORTS

DLLDESTDIR = ../../BinQt

INCLUDEPATH += ../../General

LIBS += ../../BinQt/General.dll

SOURCES += ../../Solverlib/AbstractSolver.cpp \
           ../../Solverlib/SolverPim.cpp \
           ../../Solverlib/SolverQmrcg.cpp \
           ../../Solverlib/SolverTangcg.cpp \
           ../../Solverlib/SolverZbcg2.cpp \
           ../../Solverlib/SolverSbicgm.cpp \
           ../../Solverlib/Cgcommon2.cpp \
           ../../Solverlib/CgStruct.cpp \
           ../../Solverlib/Blas.cpp

HEADERS += ../../Solverlib/Solverlib.h \
           ../../Solverlib/AbstractSolver.h \
           ../../Solverlib/SolverPim.h \
           ../../Solverlib/SolverQmrcg.h \
           ../../Solverlib/SolverTangcg.h \
           ../../Solverlib/SolverZbcg2.h \
           ../../Solverlib/SolverSbicgm.h \
           ../../Solverlib/Cgcommon2.h \
           ../../Solverlib/CgStruct.h \
           ../../Solverlib/Blas.h

unix:!symbian {
    maemo5 {
        target.path = /opt/usr/lib
    } else {
        target.path = /usr/lib
    }
    INSTALLS += target
}
