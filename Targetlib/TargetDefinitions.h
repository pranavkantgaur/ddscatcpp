#ifndef __TARGET_DEFINITIONS_H__
#define __TARGET_DEFINITIONS_H__

typedef enum _TargetType
{
	TargetType_Aniellips,	TargetType_AniEll2,		TargetType_AniEll3,		TargetType_AniEllN,
	TargetType_Ellipsoid,	TargetType_Ellipso2,	TargetType_Ellipso3,	TargetType_EllipsoN,
	TargetType_Anirctngl,	TargetType_Bislinpbc,	TargetType_Conellips,	TargetType_Cylinder1, 
	TargetType_Cylndrcap,	TargetType_Cylndrpbc,	TargetType_Dskblypbc,	TargetType_Dskrctngl,	
	TargetType_Dskrctpbc,	TargetType_Dw1996tar,	TargetType_Elinrct,		TargetType_HexPrism,
	TargetType_Hexgonpbc,	TargetType_Layrdslab,	TargetType_Lyrslbpbc,	TargetType_Mltblocks,	
	TargetType_Rctglpbc,	TargetType_Rctglprsm,	TargetType_Rctglblk3,	TargetType_Recrecpbc,	
	TargetType_Slabhole,	TargetType_Slbholpbc,	TargetType_SpheresN,	TargetType_Sphroid2,
	TargetType_Sphrn_pbc,	TargetType_SphAniN,		TargetType_Tetrahdrn,	TargetType_Trilyrpbc,	
	TargetType_Trnglprsm,	TargetType_Gausssph,	TargetType_Anifilpbc,	TargetType_Anifrmfil,	
	TargetType_From_file,	TargetType_Frmfilpbc,	TargetType_Rctgrctg,	TargetType_Uniaxicyl,
	TargetType_OctPrism,	TargetType_Octahedron,	TargetType_Icosahedron,	TargetType_Dodecahedron,
	TargetType_End
} TargetType;

#endif