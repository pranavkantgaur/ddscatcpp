#include "StdAfx.h"

#include "Enumerator.h"

FftMethod FftEnumerator(const char *title)
{
	if (!strcmp(title, "GPFAFT"))
		return FftMethod_GPFAFT;
	if (!strcmp(title, "FFTW21"))
		return FftMethod_FFTW21;
	if (!strcmp(title, "FFTMKL"))
		return FftMethod_FFTMKL;
	if (!strcmp(title, "CUDAFX"))
		return FftMethod_CUDAFX;
	return FftMethod_End;
}

const char *FftEnumerator(FftMethod method)
{
	switch(method)
	{
	case FftMethod_GPFAFT:
		return "GPFAFT";
	case FftMethod_FFTW21:
		return "FFTW21";
	case FftMethod_FFTMKL:
		return "FFTMKL";
	case FftMethod_CUDAFX:
		return "CUDAFX";
	default:
        return NULL;
	}
}

SolMethod SolEnumerator(const char *title)
{
	if (!strcmp(title, "PETRKP"))
		return SolMethod_PETRKP;
	if (!strcmp(title, "PBCGST"))
		return SolMethod_PBCGST;
	if (!strcmp(title, "PBCGS2"))
		return SolMethod_PBCGS2;
	if (!strcmp(title, "GPBICG"))
		return SolMethod_GPBICG;
	if (!strcmp(title, "QMRCCG"))
		return SolMethod_QMRCCG;
	if (!strcmp(title, "SBICGM"))
		return SolMethod_SBICGM;  

	return SolMethod_End;
}

const char *SolEnumerator(SolMethod method)
{
	switch(method)
	{
	case SolMethod_PETRKP:
		return "PETRKP"; 
	case SolMethod_PBCGST:
		return "PBCGST";
	case SolMethod_PBCGS2:
		return "PBCGS2";
	case SolMethod_GPBICG:
		return "GPBICG";
	case SolMethod_QMRCCG:
		return "QMRCCG";
	case SolMethod_SBICGM:
		return "SBICGM";
	default:
		return NULL;
	}
}

TorqMethod TorqEnumerator(const char *title)
{
	if (!strcmp(title, "NOTORQ"))
		return TorqMethod_NOTORQ;
	if (!strcmp(title, "DOTORQ"))
		return TorqMethod_DOTORQ;
	return TorqMethod_End;
}

const char *TorqEnumerator(TorqMethod method)
{
	switch(method)
	{
	case TorqMethod_NOTORQ:
		return "NOTORQ";
	case TorqMethod_DOTORQ:
		return "DOTORQ";
	default:
		return NULL;
	}
}

AlphaMethod AlphaEnumerator(const char *title)
{
	if (!strcmp(title, "GKDLDR"))
		return AlphaMethod_GKDLDR;
	if (!strcmp(title, "LATTDR"))
		return AlphaMethod_LATTDR;
	if (!strcmp(title, "FLTRCD"))
		return AlphaMethod_FLTRCD;
	return AlphaMethod_End;
}

const char *AlphaEnumerator(AlphaMethod method)
{
	switch(method)
	{
	case AlphaMethod_GKDLDR:
		return "GKDLDR";
	case AlphaMethod_LATTDR:
		return "LATTDR";
	case AlphaMethod_FLTRCD:
		return "FLTRCD";
	default:
		return NULL;
	}
}

BinflagMethod BinflagEnumerator(const char *title)
{
	if (!strcmp(title, "NOTBIN"))
		return BinflagMethod_NOTBIN;
	if (!strcmp(title, "ORIBIN"))
		return BinflagMethod_ORIBIN;
	if (!strcmp(title, "ALLBIN"))
		return BinflagMethod_ALLBIN;
	return BinflagMethod_End;
}

const char *BinflagEnumerator(BinflagMethod method)
{
	switch(method)
	{
	case BinflagMethod_NOTBIN:
		return "NOTBIN";
	case BinflagMethod_ORIBIN:
		return "ORIBIN";
	case BinflagMethod_ALLBIN:
		return "ALLBIN";
	default:
		return NULL;
	}
}

FrameCode FrameEnumerator(const char *title)
{
	if (!strcmp(title, "LFRAME"))
		return FrameCode_LFRAME;
	if (!strcmp(title, "TFRAME"))
		return FrameCode_TFRAME;
	return FrameCode_End;
}

const char *FrameEnumerator(FrameCode method)
{
	switch(method)
	{
	case FrameCode_LFRAME:
		return "LFRAME";
	case FrameCode_TFRAME:
		return "TFRAME";
	default:
		return NULL;
	}
}

const char *FrameEnumeratorVerbose(FrameCode method)
{
	switch(method)
	{
	case FrameCode_LFRAME:
		return "Laboratory Frame";
	case FrameCode_TFRAME:
		return "Target Frame";
	default:
		return "Unknown Frame";
	}
}

CdividMethod CdividEnumerator(const char *title)
{
	if (!strcmp(title, "LIN"))
		return CdividLin;
	if (!strcmp(title, "INV"))
		return CdividInv;
	if (!strcmp(title, "LOG"))
		return CdividLog;
	if (!strcmp(title, "TAB"))
		return CdividTab;
	return CdividEnd;
}

const char *CdividEnumerator(CdividMethod method)
{
	switch(method)
	{
	case CdividLin:
		return "LIN";
	case CdividInv:
		return "INV";
	case CdividLog:
		return "LOG";
	case CdividTab:
		return "TAB";
	default:
		return NULL;
	}
}
