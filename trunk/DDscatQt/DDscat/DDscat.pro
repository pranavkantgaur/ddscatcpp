TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

DESTDIR = ../../BinQt

INCLUDEPATH += ../../General \
               ../../Targetlib \
               ../../Fourierlib \
               ../../Solverlib \
               ../../Xml/include

LIBS += ../../BinQt/General.dll \
        ../../BinQt/Targetlib.dll \
        ../../BinQt/Fourierlib.dll \
        ../../BinQt/Solverlib.dll \
        ../../BinQt/libxml2.dll

SOURCES += ../../DDscat/Alphadiag.cpp \
           ../../DDscat/Cprod.cpp \
           ../../DDscat/DDscat.cpp \
           ../../DDscat/DDscatCommons.cpp \
           ../../DDscat/DDscatEngine.cpp \
           ../../DDscat/DDscatParameters.cpp \
           ../../DDscat/DielectricManager.cpp \
           ../../DDscat/DipoleData.cpp \
           ../../DDscat/FileNamer.cpp \
           ../../DDscat/Getfml.cpp \
           ../../DDscat/Getmueller.cpp \
           ../../DDscat/GreenFunctionManager.cpp \
           ../../DDscat/Matrix.cpp \
           ../../DDscat/Mkl_fake.cpp \
           ../../DDscat/Mpi_fake.cpp \
           ../../DDscat/Namid.cpp \
           ../../DDscat/Nearfield.cpp \
           ../../DDscat/Nuller.cpp \
           ../../DDscat/OutputManager.cpp \
           ../../DDscat/Pbcscavec.cpp \
           ../../DDscat/Rot2.cpp \
           ../../DDscat/Rotate.cpp \
           ../../DDscat/ScatterManager.cpp \
           ../../DDscat/Timeit.cpp \
           ../../DDscat/TimerManager.cpp \
           ../../DDscat/Version.cpp

HEADERS += ../../DDscat/Cprod.h \
           ../../DDscat/DDscatCommons.h \
           ../../DDscat/DDscatEngine.h \
           ../../DDscat/DDscatMain.h \
           ../../DDscat/DDscatParameters.h \
           ../../DDscat/DielectricManager.h \
           ../../DDscat/DipoleData.h \
           ../../DDscat/FileNamer.h \
           ../../DDscat/GreenFunctionManager.h \
           ../../DDscat/Matrix.h \
           ../../DDscat/Mkl_fake.h \
           ../../DDscat/Mpi_fake.h \
           ../../DDscat/OutputManager.h \
           ../../DDscat/ScatterManager.h \
           ../../DDscat/SumHolder.h \
           ../../DDscat/Timeit.h \
           ../../DDscat/TimerManager.h
