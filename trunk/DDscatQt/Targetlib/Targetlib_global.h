#ifndef TARGETLIB_GLOBAL_H
#define TARGETLIB_GLOBAL_H

#include <QtCore/qglobal.h>

#if defined(TARGETLIB_LIBRARY)
#  define TARGETLIBSHARED_EXPORT Q_DECL_EXPORT
#else
#  define TARGETLIBSHARED_EXPORT Q_DECL_IMPORT
#endif

#endif // TARGETLIB_GLOBAL_H
