#ifndef __ENUMERATOR_H__
#define __ENUMERATOR_H__

#include "General.h"

typedef enum _NearfieldMethod
{
	NearfieldMethodNo, NearfieldMethodPrepare, NearfieldMethodDo, NearfieldMethodEnd
} NearfieldMethod;

typedef enum _NearfieldBMethod
{
	NearfieldBMethodElectric, NearfieldBMethodBoth, NearfieldBMethodEnd
} NearfieldBMethod;

typedef enum _FftMethod
{
	FftMethod_GPFAFT, FftMethod_FFTW21, FftMethod_FFTMKL, FftMethod_CUDAFX, FftMethod_End
} FftMethod;

GENERALLIB_API FftMethod FftEnumerator(const char *title);
GENERALLIB_API const char *FftEnumerator(FftMethod method);

typedef enum _SolMethod
{
	SolMethod_PBCGS2, SolMethod_PBCGST, SolMethod_GPBICG, SolMethod_PETRKP, SolMethod_QMRCCG, SolMethod_SBICGM, SolMethod_End
} SolMethod;

GENERALLIB_API SolMethod SolEnumerator(const char *title);
GENERALLIB_API const char *SolEnumerator(SolMethod method);

typedef enum _TorqMethod
{
	TorqMethod_NOTORQ, TorqMethod_DOTORQ, TorqMethod_End
} TorqMethod;

GENERALLIB_API TorqMethod TorqEnumerator(const char *title);
GENERALLIB_API const char *TorqEnumerator(TorqMethod method);

typedef enum _AlphaMethod
{
	AlphaMethod_GKDLDR, AlphaMethod_LATTDR, AlphaMethod_FLTRCD, AlphaMethod_End
} AlphaMethod;

GENERALLIB_API AlphaMethod AlphaEnumerator(const char *title);
GENERALLIB_API const char *AlphaEnumerator(AlphaMethod method);

typedef enum _BinflagMethod
{
	BinflagMethod_NOTBIN, BinflagMethod_ORIBIN, BinflagMethod_ALLBIN, BinflagMethod_End
} BinflagMethod;

GENERALLIB_API BinflagMethod BinflagEnumerator(const char *title);
GENERALLIB_API const char *BinflagEnumerator(BinflagMethod method);

typedef enum _FrameCode
{
	FrameCode_LFRAME, FrameCode_TFRAME, FrameCode_End
} FrameCode;

GENERALLIB_API FrameCode FrameEnumerator(const char *title);
GENERALLIB_API const char *FrameEnumerator(FrameCode method);
GENERALLIB_API const char *FrameEnumeratorVerbose(FrameCode method);

typedef enum _CdividMethod
{
	CdividLin, CdividInv, CdividLog, CdividTab, CdividEnd
} CdividMethod;

GENERALLIB_API CdividMethod CdividEnumerator(const char *title);
GENERALLIB_API const char *CdividEnumerator(CdividMethod method);

typedef enum _PeriodicBoundaryFlag
{
	PeriodicNo, PeriodicY, PeriodicZ, PeriodicBoth, PeriodicEnd
} PeriodicBoundaryFlag;

typedef enum _IsotropicFlag
{
	TargetIsIsotropic, TargetIsAnisotropic, TargetIsAnisotropicDisoriented, TargetIsEnd
} IsotropicFlag;

#endif
