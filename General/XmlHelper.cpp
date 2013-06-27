#include "StdAfx.h"

#include "XmlHelper.h"

bool XmlHaveProperty(xmlNodePtr Node, const xmlChar *Name)
{
	bool bResult = false;
	xmlChar *uri = xmlGetProp(Node, Name);
	if (uri) 
		bResult = true;
	xmlFree(uri);
	return bResult;
}

int XmlGetIntProperty(xmlNodePtr Node, const xmlChar *Name)
{
	int Result = 0;
	xmlChar *uri = xmlGetProp(Node, Name);
	if (uri)
	{
		int len = strlen((char *)uri);
		if (len) 
			Result = Scaner::ScanInt((const char *)uri, len);  
		xmlFree(uri);
	}
	return Result;
}

int XmlGetHexIntProperty(xmlNodePtr Node, const xmlChar *Name)
{
	int Result = 0;
	xmlChar *uri = xmlGetProp(Node, Name);
	if (uri)
	{
		int len = strlen((char *)uri);
		Result = Scaner::ScanIntHex((const char *)uri, len);  
		xmlFree(uri);
	}
	return Result;
}

real XmlGetDoubleProperty(xmlNodePtr Node, const xmlChar *Name)
{
	real Result = (real)0.;
	xmlChar *uri = xmlGetProp(Node, Name);
	if (uri)
	{
		int len = (int)strlen((char *)uri);
		Result = Scaner::ScanDouble((char *)uri, len);  
		xmlFree(uri);
	}
	return Result;
}

real XmlGetExponentProperty(xmlNodePtr Node, const xmlChar *Name)
{
	real Result = (real)0.;
	xmlChar *uri = xmlGetProp(Node, Name);  
	if (uri)
	{
		int len = strlen((char *)uri);
		Result = Scaner::ScanExponent((char *)uri, len);  
		xmlFree(uri);
	}
	return Result;
}

bool XmlGetBoolProperty(xmlNodePtr Node, const xmlChar *Name)
{
	bool Result = false;
	xmlChar *uri = xmlGetProp(Node, Name);
	if (uri)
	{
		char tmp = (char)tolower(uri[0]);
		if ((tmp == 'f') || (tmp == 'n') || (tmp == '0')) Result = false;
		if ((tmp == 't') || (tmp == 'y') || (tmp == '1')) Result = true;
		xmlFree(uri);
	}
	return Result;
}

bool XmlHaveChild(const xmlNodePtr Father, const xmlChar *Name)
{
	xmlNodePtr Cur = Father->xmlChildrenNode;
	while(Cur)
	{
		if (Cur->type == XML_ELEMENT_NODE) 
			if (!strcmp((const char *)Cur->name, (const char *)Name)) 
				return true;
		Cur = Cur->next;
	}
	return false;
}

xmlNodePtr XmlFindChild(const xmlNodePtr Father, const xmlChar *Name)
{
	xmlNodePtr Cur = Father->xmlChildrenNode;
	while(Cur)
	{
		if (Cur->type == XML_ELEMENT_NODE) 
			if (!strcmp((const char *)Cur->name, (const char *)Name)) 
				return Cur;
		Cur = Cur->next;
	}
	return NULL;
}

xmlNodePtr XmlFindChildWithProp(const xmlNodePtr Father, const xmlChar *Name, const xmlChar *Prop, const xmlChar *PropValue)
{
	xmlNodePtr Cur = Father->xmlChildrenNode;
	while(Cur)
	{
		if (Cur->type == XML_ELEMENT_NODE) 
		{
			if (!strcmp((const char *)Cur->name, (const char *)Name))
			{
				xmlChar *uri = xmlGetProp(Cur, Prop);
				if (!uri) return NULL;
				if (!strcmp((const char *)uri, (const char *)PropValue))
				{
					xmlFree(uri);
					return Cur;
				}
				xmlFree(uri);
			}
		}
		Cur = Cur->next;
	}
	return NULL;
}
