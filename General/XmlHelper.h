#ifndef __XMLHELPER_H__
#define __XMLHELPER_H__

#include "Definitions.h"
#include "Scaner.h"
#include "libxml/xmlmemory.h"

GENERALLIB_API bool XmlHaveProperty(xmlNodePtr Node, const xmlChar *Name);

GENERALLIB_API int XmlGetIntProperty(xmlNodePtr Node, const xmlChar *Name);
GENERALLIB_API int XmlGetHexIntProperty(xmlNodePtr Node, const xmlChar *Name);
GENERALLIB_API real XmlGetDoubleProperty(xmlNodePtr Node, const xmlChar *Name);
GENERALLIB_API real XmlGetExponentProperty(xmlNodePtr Node, const xmlChar *Name);
GENERALLIB_API bool XmlGetBoolProperty(xmlNodePtr Node, const xmlChar *Name);

GENERALLIB_API bool XmlHaveChild(const xmlNodePtr Father, const xmlChar *Name);
GENERALLIB_API xmlNodePtr XmlFindChild(const xmlNodePtr Father, const xmlChar *Name); 
GENERALLIB_API xmlNodePtr XmlFindChildWithProp(const xmlNodePtr Father, const xmlChar *Name, const xmlChar *Prop, const xmlChar *PropValue);

#endif // __XMLHELPER_H__

