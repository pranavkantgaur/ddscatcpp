#ifndef GENERAL_GLOBAL_H
#define GENERAL_GLOBAL_H

#include <QtCore/qglobal.h>

#if defined(GENERAL_LIBRARY)
#  define GENERALSHARED_EXPORT Q_DECL_EXPORT
#else
#  define GENERALSHARED_EXPORT Q_DECL_IMPORT
#endif

#endif // GENERAL_GLOBAL_H
