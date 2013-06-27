#ifndef VTRLIB_GLOBAL_H
#define VTRLIB_GLOBAL_H

#include <QtCore/qglobal.h>

#if defined(VTRLIB_LIBRARY)
#  define VTRLIBSHARED_EXPORT Q_DECL_EXPORT
#else
#  define VTRLIBSHARED_EXPORT Q_DECL_IMPORT
#endif

#endif // VTRLIB_GLOBAL_H
