#-------------------------------------------------
#
# Project created by QtCreator 2013-01-14T18:01:38
#
#-------------------------------------------------

QT       -= gui

TARGET = Targetlib
TEMPLATE = lib

DEFINES += TARGETLIB_LIBRARY \
           TARGETLIB_EXPORTS

DLLDESTDIR = ../../BinQt

INCLUDEPATH += ../../General

LIBS += ../../BinQt/General.dll

SOURCES += ../../Targetlib/AbstractTarget.cpp \
           ../../Targetlib/AbstractDFData.cpp \
           ../../Targetlib/Item0.cpp \
           ../../Targetlib/Line.cpp \
           ../../Targetlib/TargetManager.cpp \
           ../../Targetlib/LoadableTarget.cpp \
           ../../Targetlib/Tar2el.cpp \
           ../../Targetlib/Tar3el.cpp \
           ../../Targetlib/Tar2sp.cpp \
           ../../Targetlib/Taranirec.cpp \
           ../../Targetlib/Tarblocks.cpp \
           ../../Targetlib/Tarcel.cpp \
           ../../Targetlib/Tarcyl.cpp \
           ../../Targetlib/Tarcylcap.cpp \
           ../../Targetlib/Tarell.cpp \
           ../../Targetlib/Targspher.cpp \
           ../../Targetlib/Tarhex.cpp \
           ../../Targetlib/Tarlyrslab.cpp \
           ../../Targetlib/Tarnas.cpp \
           ../../Targetlib/TarNel.cpp \
           ../../Targetlib/Tarnsp.cpp \
           ../../Targetlib/Taroct.cpp \
           ../../Targetlib/Tarpbxn.cpp \
           ../../Targetlib/TarPolyhedra.cpp \
           ../../Targetlib/Tarprsm.cpp \
           ../../Targetlib/Tarrctblk3.cpp \
           ../../Targetlib/Tarrctell.cpp \
           ../../Targetlib/Tarrec.cpp \
           ../../Targetlib/Tarrecrec.cpp \
           ../../Targetlib/Tarslabhole.cpp \
           ../../Targetlib/Tartet.cpp \
           ../../Targetlib/Target_From_File.cpp

HEADERS += ../../Targetlib/Targetlib.h \
           ../../Targetlib/TargetDefinitions.h \
           ../../Targetlib/AbstractTarget.h \
           ../../Targetlib/AbstractDFData.h \
           ../../Targetlib/Item0.h \
           ../../Targetlib/Line.h \
           ../../Targetlib/TargetManager.h \
           ../../Targetlib/LoadableTarget.h \
           ../../Targetlib/Tar2el.h \
           ../../Targetlib/Tar3el.h \
           ../../Targetlib/Tar2sp.h \
           ../../Targetlib/Taranirec.h \
           ../../Targetlib/Tarblocks.h \
           ../../Targetlib/Tarcel.h \
           ../../Targetlib/Tarcyl.h \
           ../../Targetlib/Tarcylcap.h \
           ../../Targetlib/Tarell.h \
           ../../Targetlib/Targspher.h \
           ../../Targetlib/Tarhex.h \
           ../../Targetlib/Tarlyrslab.h \
           ../../Targetlib/Tarnas.h \
           ../../Targetlib/TarNel.h \
           ../../Targetlib/Tarnsp.h \
           ../../Targetlib/Taroct.h \
           ../../Targetlib/Tarpbxn.h \
           ../../Targetlib/TarPolyhedra.h \
           ../../Targetlib/Tarprsm.h \
           ../../Targetlib/Tarrctblk3.h \
           ../../Targetlib/Tarrctell.h \
           ../../Targetlib/Tarrec.h \
           ../../Targetlib/Tarrecrec.h \
           ../../Targetlib/Tarslabhole.h \
           ../../Targetlib/Tartet.h \
           ../../Targetlib/Target_From_File.h

unix:!symbian {
    maemo5 {
        target.path = /opt/usr/lib
    } else {
        target.path = /usr/lib
    }
    INSTALLS += target
}
