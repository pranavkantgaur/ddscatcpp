#ifndef SOLVERLIB_GLOBAL_H
#define SOLVERLIB_GLOBAL_H

#include <QtCore/qglobal.h>

#if defined(SOLVERLIB_LIBRARY)
#  define SOLVERLIBSHARED_EXPORT Q_DECL_EXPORT
#else
#  define SOLVERLIBSHARED_EXPORT Q_DECL_IMPORT
#endif

#endif // SOLVERLIB_GLOBAL_H
