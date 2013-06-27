	/* **int DDscatParameters::LoadXml(const char *fileName, int &nrfld)
{
	printf("' ========= LaodXml ==================='\n");
	const char *ReaparLabel = "Reapar";


		
	xmlNodePtr nod = xmlDocGetRootElement(doc);//верхний нод (корень)
		
	printf("' ========= Parameter file for v%s", xmlGetProp(nod, (const xmlChar *)"ver") );//аттрибут DDScatParameterFile (ver)
	printf(" ==================='\n\n");
	
	nod = nod->children;
	while(xmlStrcmp(nod->name,(const xmlChar *) "Preliminaries")) nod = nod->next;
	
	printf("'**** Preliminaries ****'\n\n");

	strncpy(cmdtrq, (char *) xmlGetProp(nod, (const xmlChar *)"CMTORQ"), 6);
	cmdtrq[6] = '\0';
	if(*cmdtrq == *"DOTORQ")
	{
		printf("Reapar DOTORQ - compute torques \n");
	}
	if (*cmdtrq == *"NOTORQ")
	{
		printf("Reapar NOTORQ - do not compute torques \n");
	}	
//
// Define method used for iterative solution of Complex linear equations
//	CMDSOL*6	= 'PETRKP' -  Petravic & Kuo-Petravic method
//				= 'PBCGST' -  PIM BiConjugate Gradient with Stabilization
//				= 'PBCGS2' -  M.A Botcheve implementation of BiCGstab enhanced to improve convergence properties with finite precision arithmetic
//				= 'GPBICG' -  Tang et al conjugate gradient solver implemented by P.C. Chaumet & A. Rahmani
//				= 'QMRCCG' -  PIM interface to QMR solver implemented by P.C. Chaumet and A. Rahmani
	strncpy(cmdsol, (char *) xmlGetProp(nod, (const xmlChar *)"CMDSOL"), 6);
	cmdsol[6] = '\0';
	printf("Reapar %s - CCG Method  \n", cmdsol);
//
// Define FFT method:
//	CMDFFT*6	= 'GPFAFT' -  GPFA code of Temperton
//				= 'FFTW21' -  FFTW 2.1.x code of Frigo & Johnson
//				= 'FFTMKL' -  Use DFTI from Intel MKL
	strncpy(cmdfft, (char *) xmlGetProp(nod, (const xmlChar *)"CMDFFT"), 6);
	cmdfft[6] = '\0';
	if (*cmdfft == *"FFTW21")
	{
		printf("Reapar FFTW21 - using FFTW 2.1.x package from Frigo & Johnson\n");
	}
	if (*cmdfft == *"GPFAFT")
	{
		printf("Reapar GPFAFT - using GPFA package from Clive Temperton\n");
	}
	if (*cmdfft == *"FFTMKL")
	{
		printf("Reapar FFTMKL - using DFTI from Intel Math Kernel Library (MKL)\n");
	}
//
// Define prescription for computing polarizabilities:
//	CALPHA*6	= 'LATTDR' - Lattice Dispersion Relation of Draine & Goodman (1993)
//				= 'GKDLDR' - Lattice Dispersion Relation of Gutkowicz-Krusin & Draine (2004)
	strncpy(calpha, (char *) xmlGetProp(nod, (const xmlChar *)"CALPHA"), 6);
	calpha[6] = '\0';
	if (*calpha == *"LATTDR")
	{
		printf("Reapar LATTDR - Draine & Goodman (1993) LDR for alpha\n");
	}
	if (*calpha == *"GKDLDR")
	{
		printf("Reapar GKDLDR - Gutkowicz-Krusin & Draine (2004) LDR for alpha\n");
	}
//
// Binary file flag
	strncpy(cbinflag, (char *) xmlGetProp(nod, (const xmlChar *)"CBINFLAG"), 6);
	cbinflag[6] = '\0';
	printf("Reapar %s - Unformatted binary dump option\n", cbinflag);
//
// Read upper bound on target extent

	while(xmlStrcmp(nod->name,(const xmlChar *) "InitialMemoryAllocation")) nod = nod->next;
	
	mxnx = atoi( (char *) xmlGetProp(nod, (const xmlChar *)"X") );
	mxny = atoi( (char *) xmlGetProp(nod, (const xmlChar *)"Y") );
	mxnz = atoi( (char *) xmlGetProp(nod, (const xmlChar *)"Z") );

	printf("Reapar '**** Initial Memory Allocation ****'\n \nReapar allow MXNX, MXNY, MXNZ = %d %d %d for target generation", mxnx, mxny, mxnz);

// 
//	CSHAPE*9	= 'FROM_FILE' shape and composition data will later be read from file CFLSHP for dielectric tensors that are diagonal in Target Frame
//				= 'ANIFRMFIL' read shape and composition data from file, for general anisotropic dielectric tensors with orientation angles 
//							  THETADF,PHIDF,BETAD relative to Target Frame.
//				= 'FRMFILPBC' shape and composition data for TUC will later be read from file CFLSHP for dielectric tensors that are diagonal in Target Frame 
//							  PYD and PZD are input via ddscat.par
//				= 'ANIFILPBC' shape and composition data for TUC will later be read from file CFLSHP for general anisotropic dielectricc tensors with orientation 
//							  ngles THETADF,PHIDF,BETAD relative to Target Frame. PYD and PZD are input via ddscat.par
//				= 'ANIELLIPS' ellipsoid of anisotropic material
//				= 'ANIRCTNGL' homogeneous anisotropic rectangular target
//				= 'ANI_ELL_2' two touching anisotropic ellipsoids of materials 1-6
//				= 'ANI_ELL_3' three touching anisotropic ellipsoids of materials 1-9
//				= 'BISLINPBC' bilayer slab with periodic grid of lines parallel to z on top, with y-period/d=PYD [or, if PYD=0, a single line]
//				= 'CONELLIPS' two concentric ellipsoids of materials 1,2
//				= 'CYLINDER1' homogeneous finite cylinder
//				= 'CYLNDRCAP' homogeneous cylinder with hemispherical endcaps
//				= 'CYLNDRPBC' 1-d or 2-d array of finite cylinders
//				= 'DSKBLYPBC' 1-d or 2-d array of disk on bilayer rect. slab
//				= 'DSKRCTNGL' single disk on rectangular slab
//				= 'DSKRCTPBC' 1-d or 2-d array of disk on rectangular slab
//				= 'DW1996TAR' 13-cube target used by Draine & Weingartner 1996
//				= 'EL_IN_RCT' ellipsoid embedded in rectangular block
//				= 'ELLIPSOID' ellipsoid (homogeneous and isotropic)
//				= 'ELLIPSPBC' 1-d or 2-d array of ellipsoids
//				= 'ELLIPSO_2' two touching isotropic ellipsoids of materials 1 and 2
//				= 'ELLIPSO_3'  three touching isotropic ellipsoids of materials 1,2,3
//				= 'GAUSS_SPH'  gaussian sphere target
//				= 'HEXGONPBC'  1-d or 2-d array of hexagonal prisms
//				= 'HEX_PRISM'  homogeneous hexagonal prism
//				= 'LYRD_SLAB'  layered slab target, with up to 4 separate material layers
//				= 'LYRSLBPBC'  1-d or 2-d array of layered rect. slab targets, with up to 4 material layers
//				= 'MLTBLOCKS'  collection of cubic blocks defined by data in file 'blocks.par'
//				= 'RCTGLBLK3'  isolated target: 3 rectangular blocks with centers on x-axis
//				= 'RCTGLPRSM'  homogeneous rectangular prism
//				= 'RCTGL_PBC'  1-d or 2-d array of rectangular targets
//				= 'SLAB_HOLE'  rectangular block with cylindrical hole
//				= 'SLBHOLPBC'  1-d or 2-d array of rectangular blocks with cylindrical hole
//				= 'SPHERES_N'  multisphere target = union of N spheres
//				= 'SPHRN_PBC'  1-d or 2-d array of multisphere target
//				= 'SPHROID_2'  two touching spheroids with symmetry axes at specified angle!
//				= 'SPH_ANI_N'  multisphere target, with spheres that can have different, anisotropic, composition
//				= 'TETRAHDRN'  regular tetrahedron
//				= 'TRILYRPBC'  periodic target: 3 layer rectangular structure
//				= 'TRNGLPRSM'  triangular prism (homogeneous and isotropic)
//				= 'UNIAXICYL'  cylinder of unixaxial material
	while(xmlStrcmp(nod->name,(const xmlChar *) "TargetGeometryAndComposition")) nod = nod->next;

	printf("\nReapar '**** Target Geometry and Composition ****'\n\n");
	strncpy(cshape, (char *) xmlGetProp(nod, (const xmlChar *)"CSHAPE"), 9);
	cshape[9] = '\0';
	strlwr(cshape);

	ShapeNumber shape = Enumerator(cshape);

	if (shape != Shape_End)
	{
		strupr(cshape);	
		printf("Reapar %s - Shape definition \n", cshape);
		strlwr(cshape);
	}
	else
	{
		printf("Reapar CSHAPE = %S  Unrecognized shape directive\n", cshape);
	}
//
// Read shape parameters
	nod = nod->children;

	char *ch;

	//shape = Shape_uniaxicyl;
	switch(shape)											// "anifilpbc", "ellipspbc", "sphrn_pbc"			<- not represented ?
	{
	case Shape_anifrmfil:			// ANIFRMFIL
	case Shape_from_file:			// FROM_FILE
	case Shape_mltblocks:			// MLTBLOCKS
		numShpar = 0;						// targets needing 0 shape parameters (skip one line)
		break;

	case Shape_dw1996tar:			// DW1996TAR
	case Shape_tetrahdrn:			// TETRAHDRN
		numShpar = 1;						// targets needing 1 shape parameter
		while(xmlStrcmp(nod->name,(const xmlChar *) "SP")) nod = nod->next;
		ch = (char*) xmlGetProp(nod, nod->properties->name);
		shpar[0] = atof(ch);
		printf("Reapar %f = sphar[0]\n",shpar[0]);		
		break;

	case Shape_cylndrcap:	// CYLNDRCAP
	case Shape_uniaxicyl:	// UNIAXICYL						
		numShpar = 2;						// targets needing 2 shape parameters
		TakeSP(numShpar, nod);	
		break;

	case Shape_frmfilpbc:	// FRMFILPBC
	case Shape_spheres_n:	// SPHERES_N
	case Shape_sph_ani_n:	// SPH_ANI_N
		{
			numShpar = 2;					// targets needing 2 shape parameters and file name
			TakeSP(numShpar, nod);
			free (cflshp);

			nod = nod->next;
			while(xmlStrcmp(nod->name,(const xmlChar *) "SP")) nod = nod->next;
			
			cflshp = (char *)malloc(strlen((char*) xmlGetProp(nod, nod->properties->name)) + 1);
			strcpy(cflshp, (char*) xmlGetProp(nod, nod->properties->name));

			printf("Reapar shape file= %s\n", cflshp);
		}
		break;

	case Shape_aniellips:	// ANIELLIPS
	case Shape_ani_ell_2:	// ANI_ELL_2
	case Shape_ani_ell_3:	// ANI_ELL_3
	case Shape_anirctngl:	// ANIRCTNGL
	case Shape_cylinder1:	// CYLINDER1
	case Shape_ellipsoid:	// ELLIPSOID
	case Shape_ellipso_2:	// ELLIPSO_2
	case Shape_ellipso_3:	// ELLIPSO_3
	case Shape_hex_prism:	// HEX_PRISM
	case Shape_rctglprsm:	// RCTGLPRSM					// targets needing 3 shape parameters
		numShpar = 3;
		
		TakeSP(numShpar, nod);
		break;

	case Shape_sphrn_pbc:		// SPHRN_PBC
		{
			
			numShpar = 3;					// targets needing 3 shape parameters and file name
			TakeSP(numShpar, nod);

			nod = nod->next;
			while(xmlStrcmp(nod->name,(const xmlChar *) "SP")) nod = nod->next;
			

			free (cflshp);
			cflshp = (char *)malloc(strlen((char*) xmlGetProp(nod, nod->properties->name)) + 1);
			strcpy(cflshp, (char*) xmlGetProp(nod, nod->properties->name));

			printf("Reapar shape file= %s\n", cflshp);
		}
		break;

	case Shape_slab_hole:	// SLAB_HOLE
	case Shape_trnglprsm:	// TRNGLPRSM
		numShpar = 4;					// targets needing 4 shape parameters 
		TakeSP(numShpar, nod);
		break;
	case Shape_cylndrpbc:	// CYLNDRPBC
	case Shape_dskrctngl:	// DSKRCTNGL
	case Shape_hexgonpbc:	// HEXGONPBC
	case Shape_rctgl_pbc:	// RCTGL_PBC
		numShpar = 5;							// targets needing 5 shape parameters
		TakeSP(numShpar, nod);
		break;
	case Shape_bislinpbc:	// BISLINPBC
	case Shape_conellips:	// CONELLIPS
	case Shape_el_in_rct:	// EL_IN_RCT
	case Shape_gauss_sph:	// GAUSS_SPH
	case Shape_rctg_rctg:	// RCTG_RCTG
	case Shape_slbholpbc:	// SLBHOLPBC					
		numShpar = 6;											// targets needing 6 shape parameters
		TakeSP(numShpar, nod);
		break;
	case Shape_dskrctpbc:	// DSKRCTPBC
	case Shape_layrdslab:	// LAYRDSLAB
		numShpar = 7;						// targets needing 7 shape parameters
		TakeSP(numShpar, nod);
		break;
	case Shape_dskblypbc:	// DSKBLYPBC
	case Shape_recrecpbc:	// RECRECPBC
								// targets needing 8 shape parameters
		numShpar = 8;
		TakeSP(numShpar, nod);
		break;

	case Shape_lyrslbpbc:	// LYRSLBPBC
	case Shape_rctglblk3:	// RCTGLBLK3
								// targets needing 9 shape parameters
		numShpar = 9;
		TakeSP(numShpar, nod);
		break;

// targets needing 10 shape parameters: none

	case Shape_trilyrpbc:	// TRILYRPBC				// targets needing 11 shape parameters
		numShpar = 11;
		TakeSP(numShpar, nod);
		break;

	default:
		Errmsg("Fatal", ReaparLabel, " No instructions for reading shape parameters for this shape");
		break;
	}

//	
// Need to set JPBC for PBC targets
// JPBC	= 0 if PBC are not used
//		= 1 if PBC in y direction only
//		= 2 if PBC in z direction only
//		= 3 if PBC in y and z directions
	jpbc = 0;
	switch(shape)
	{
	case Shape_bislinpbc:	// BISLINPBC
		if (shpar[5] > zero_) jpbc = 1;
		jpbc += 2;
	  break;

	case Shape_anifilpbc:	// ANIFILPBC
		if (shpar[0] > zero_) jpbc = 1;
		if (shpar[1] > zero_) jpbc += 2;
	  break;

	case Shape_cylndrpbc:	// CYLNDRPBC
	case Shape_hexgonpbc:	// HEXGONPBC
	case Shape_rctgl_pbc:	// RCTGL_PBC
		if (shpar[3] > zero_) jpbc = 1;
		if (shpar[4] > zero_) jpbc += 2;
	  break;

	case Shape_dskblypbc:	// DSKBLYPBC
	case Shape_recrecpbc:	// RECRECPBC
		if (shpar[6] > zero_) jpbc = 1;
		if (shpar[7] > zero_) jpbc += 2;
	  break;

	case Shape_dskrctpbc:	// DSKRCTPBC
		if (shpar[5] > zero_) jpbc = 1;
		if (shpar[6] > zero_) jpbc += 2;
	  break;

	case Shape_frmfilpbc:	// FRMFILPBC
		if (shpar[0] > zero_) jpbc = 1;
		if (shpar[1] > zero_) jpbc += 2;
	  break;

	case Shape_lyrslbpbc:	// LYRSLBPBC
		if (shpar[7] > zero_) jpbc = 1;
		if (shpar[8] > zero_) jpbc += 2;
	  break;

	case Shape_slbholpbc:	// SLBHOLPBC
		if (shpar[4] > zero_) jpbc = 1;
		if (shpar[5] > zero_) jpbc += 2;
	  break;

	case Shape_sphrn_pbc:	// SPHRN_PBC
		if (shpar[1] > zero_) jpbc = 1;
		if (shpar[2] > zero_) jpbc += 2;
	  break;

	case Shape_trilyrpbc:	// TRILYRPBC
		if (shpar[9]  > zero_) jpbc = 1;
		if (shpar[10] > zero_) jpbc += 2;
	  break;

	default:
		break;
	}
//
// following code disabled 03.01.29
// (retain for use in future noncubic version)
// 
//   Obtain lattice anisotropy parameters DX(1-3)
//   For cubic lattice, DX(1)=DX(2)=DX(3)
//   Note that we do not require here that DX(1)*DX(2)*DX(3)=1 upon
//   input; DX is renormalized here before being returned to DDSCAT
// 
//      READ(IOPAR,FMT=*,ERR=99)DX(1),DX(2),DX(3)
//      WRITE(CMSGNM,FMT='(F8.3,F8.3,F8.3,A)')DX,' = relative lattice spacings dx,dy,dz'
//      CALL WRIMSG('REAPAR',CMSGNM)
// 
//      DELTA=(DX(1)*DX(2)*DX(3))**(1./3.)
//      DX(1)=DX(1)/DELTA
//      DX(2)=DX(2)/DELTA
//      DX(3)=DX(3)/DELTA
//      WRITE(CMSGNM,FMT='(F8.3,F8.3,F8.3,A)')DX,' = normalized lattice spacings dx,dy,dz'
//      CALL WRIMSG('REAPAR',CMSGNM)
// and replaced by following:
	dx.Set(onex_, onex_, onex_);
//
// Obtain names of file(s) containing dielectric function(s)
//	NCOMP = number of different dielectric functions
//	CFLEPS(1-NCOMP) = names of files containing dielectric functions
	nod = nod->parent;
	ch = (char*) xmlGetProp(nod, (const xmlChar *)"NCOMP");
	ncomp = atof(ch);
	printf("Reapar NCOM P= %d\n", ncomp);
	
//
// *** Check that NCOMP=2 if CSHAPE=UNIAXICYL *******************************
//                      3           ANIELLIPS
//                      2           ELLIPSO_2
//                      3           ELLIPSO_3
//                      6           ANI_ELL_2
//                      9           ANI_ELL_3
//                      2           CONELLIPS
//                      2           SPHROID_2
//                      2-4         LYRD_SLAB
//                      1-4         LYRSLBPBC

	int j;
	switch(shape)
	{
	case Shape_uniaxicyl:	// UNIAXICYL
		if (ncomp != 2)
		{
			printf("Reapar NCOMP must be 2 for option UNIAXICYL");
			return lineCounter;
		}
		break;

	case Shape_aniellips:		// ANIELLIPS
		if (ncomp != 3)
		{
			printf("Reapar NCOMP must be 3 for option ANIELLIPS");
			return lineCounter;
		}
		break;

	case Shape_ellipso_2:	// ELLIPSO_2
		if (ncomp != 2)
		{
			printf("Reapar NCOMP must be 2 for option ELLIPSO_2");
			return lineCounter;
		}
		break;

	case Shape_ellipso_3:	// ELLIPSO_3
		if (ncomp != 3)
		{
			printf("Reapar NCOMP must be 3 for option ELLIPSO_3");
			return lineCounter;
		}
		break;

	case Shape_ani_ell_2:		// ANI_ELL_2
		if (ncomp != 6)
		{
			printf("Reapar NCOMP must be 6 for option ANI_ELL_2");
			return lineCounter;
		}
		break;

	case Shape_ani_ell_3:	// ANI_ELL_3
		if (ncomp != 9)
		{
			printf("Reapar NCOMP must be 3 for option ANI_ELL_3");
			return lineCounter;
		}
		break;

	case Shape_conellips:		// CONELLIPS
		if (ncomp != 2)
		{
			printf("Reapar NCOMP must be 2 for option CONELLIPS");
			return lineCounter;
		}
		break;

	case Shape_sphroid_2:	// SPHROID_2
		if (ncomp != 2)
		{
			printf("Reapar NCOMP must be 2 for option SPHROID_2");
			return lineCounter;
		}
		break;

	case Shape_lyrslbpbc:	// LYRSLBPBC
		j = 4;
		if (shpar[6] <= zero_) 
		{
			j = 3;
			if (shpar[5] <= zero_) 
			{
				j = 2;
				if (shpar[4] <= zero_) 
					j = 1;
			}
		}
		
		if (ncomp < j) 
		{
			printf("Reapar LYRSLBPBC: need ncomp=%d", j);
			return lineCounter;
		}
	  break;

	case Shape_layrdslab:	// LYRD_SLAB
		j = 4;
		if (shpar[6] <= zero_)
		{
			j = 3;
			if (shpar[5] <= zero_)
			{
				j = 2;
				if (shpar[5] <= zero_)
					j = 1;
			}
		}
		
		if(ncomp < j)
		{
			printf("Reapar LYRD_SLAB: need ncomp=%d", j);
			return lineCounter;
		}
	  break;

	case Shape_spheres_n:	// SPHERES_N
		if (ncomp > 1)
		{
			printf("Reapar SPHERES_N is for single isotropic  material -- try option SPH_ANI_N ?");
			return lineCounter;
		}
		break;

	default:
		break;
	}
//
//	
	nod = nod->children;
	for(j=0; j<ncomp; ++j)
	{
		nod = nod->next;
		while(xmlStrcmp(nod->name,(const xmlChar *) "FWRI")) nod = nod->next;
		ch = (char*) xmlGetProp(nod, (const xmlChar *)"FILE");
		cfleps.push_back(string(ch));
		printf("Reapar 0 %s\n", ch);
	}
	nod = nod->parent;
	
	////////////////////////////////////////////////////////////////////
//
// Specify whether NEARFIELD calculation is to be done
// and specify fractional expansion of computational volume
//
	while(xmlStrcmp(nod->name,(const xmlChar *) "NearfieldCalculation")) nod = nod->next;

	ch = (char*) xmlGetProp(nod, (const xmlChar *)"NRFLD");

	j = atoi(ch);
	switch(j)
	{
	case 0:
		printf("Reapar 0 = nrfld : nearfield calculation not desired\n");
		break;

	case 1:
		printf("Reapar 1 = nrfld : calculate nearfield in specified volume\n");
		break;

	default:
		printf("Reapar %6d = nrfld : calculate nearfield in specified volume\n", j);
		return lineCounter;
		break;
	}
	nrfld += j;
//
// Upon return to DDSCAT:
// NRFLD	= 0: no interest in nearfield calculation
//			= 1: prepare to do nearfield calculation on next pass
//			= 2: perform nearfield calculation
	
	if(nrfld > 0)
	{
		//extendxyz.Load(Buffer, "%lf"); 
		//char Buffer[255]="";
		nod = nod->children;

		while(xmlStrcmp(nod->name,(const xmlChar *) "Xm")) nod = nod->next;
		nod = nod->children;
		ch = (char*)nod->content;
		extendxyz.Data(0) = atof(ch);
		//strcat(Buffer, ch);
		//strcat(Buffer, " ");
		nod = nod->parent;

		while(xmlStrcmp(nod->name,(const xmlChar *) "Xp")) nod = nod->next;
		nod = nod->children;
		ch = (char*)nod->content;
		extendxyz.Data(1) = atof(ch);
		nod = nod->parent;

		while(xmlStrcmp(nod->name,(const xmlChar *) "Ym")) nod = nod->next;
		nod = nod->children;
		ch = (char*)nod->content;
		extendxyz.Data(2) = atof(ch);
		nod = nod->parent;

		while(xmlStrcmp(nod->name,(const xmlChar *) "Yp")) nod = nod->next;
		nod = nod->children;
		ch = (char*)nod->content;
		extendxyz.Data(3) = atof(ch);
		nod = nod->parent;

		while(xmlStrcmp(nod->name,(const xmlChar *) "Zm")) nod = nod->next;
		nod = nod->children;
		ch = (char*)nod->content;
		extendxyz.Data(4) = atof(ch);
		nod = nod->parent;

		while(xmlStrcmp(nod->name,(const xmlChar *) "Zp")) nod = nod->next;
		nod = nod->children;
		ch = (char*)nod->content;
		extendxyz.Data(5) = atof(ch);
		nod = nod->parent;
		
		//extendxyz.Load(Buffer, "%lf");
	}
	else
	{
		extendxyz.Clear();
	}

	printf("Reapar %6.3lf%6.3lf%6.3lf%6.3lf%6.3lf%6.3lf\n", extendxyz.Data(0), extendxyz.Data(1), extendxyz.Data(2), extendxyz.Data(3), extendxyz.Data(4), extendxyz.Data(5));
	printf("Reapar fractional extension in -x,+x,-y,+y,-z,+z directions\n", j);
//
//     Define INITIALIZATION:
//     INIT = 0 to start with |X0>=0
//            1 to obtain |X0> from 4th-order expansion in polarizability
//            2 to read |X0> from file solvp.in
//  disabled 08.03.12 v7.0.5 since we always use INIT=0
//  skip line:
//       READ(IOPAR,FMT=*,ERR=99)CLINE
//       CALL WRIMSG(' ',CLINE)
//       CWHERE='error reading INIT in ddscat.par'
//       READ(IOPAR,FMT=*,ERR=99)INIT
//       IF(INIT==0)THEN
//          CALL WRIMSG('REAPAR','INIT=0 to start with |X0>=0 (CCG method)')
//       ELSEIF(INIT==1)THEN
//          CALL WRIMSG('REAPAR','INIT=1 to obtain |X0> from 4th-order expansion (CCG)')
//       ELSEIF(INIT==2)THEN
//          CALL WRIMSG('REAPAR','INIT=2 to read |X0> from file solvp.in')
//       ELSE
//          CALL ERRMSG('FATAL','REAPAR',' Wrong value of INIT')
//       ENDIF
// !***********************************************************************

//     Define error tolerance:
//     TOL= maximum acceptable value of |Ax-E|/|E|
	printf("Reapar '**** Error Tolerance ****'\n\n");
	const real tolLowerLimit = (real)1.e-5;
	while(xmlStrcmp(nod->name,(const xmlChar *) "TOL")) nod = nod->next;
	nod = nod->children;
	ch = (char*) nod->content;
	nod = nod->parent;
	tol = atof(ch);
	if (tol < tolLowerLimit)
	{
		tol = tolLowerLimit;
		printf("Reapar TOL corrected to be on lower limit.\n");
	}
	printf("Reapar %10.3e = TOL = max. acceptable normalized residual |Ax-E|/|E|\n", tol);
//
	printf("\n'Reapar **** maximum number of iterations allowed ****'\n\n");
	while(xmlStrcmp(nod->name,(const xmlChar *) "MXITER")) nod = nod->next;
	nod = nod->children;
	ch = (char*) nod->content;
	nod = nod->parent;
	mxiter = atoi(ch);
	printf("Reapar %10d = MXITER\n", mxiter);
//
// 2009.08.27 (BTD) add code to catch problem with
// we expect TOL to be between 1e-10 and <1.
// if outside this range, there has probably been an error reading TOL
	if(tol < 1.e-10 || tol > onex_)
	{
// Note: if the number of diel.fn. files is less than NCOMP or greater than NCOMP, this will lead to subsequent errors
//       in reading ddscat.par, which will probably affect reading of TOL
		printf("Reapar %10.3e = TOL", tol);
		printf("Reapar Appears that there has been an error reading TOL");
		printf("Reapar Check whether there are NCOMP diel.fn. files");
		printf("Reapar error reading ddscat.par file");
        return lineCounter;
	}
//	
	printf("Reapar %10.3le  = TOL = max. acceptable normalized residue |Ax-E|/|E|", tol);
	if (tol < tolLowerLimit || tol > 1.e-1)
	{
		printf("Reapar strange value of tol ");
		return lineCounter;
	}
//
// Define summation limit parameter GAMMA for PBC calculations summations will be carried out to distance r=2/(k*alpha)
// with suppression factor exp[-(alpha*kr)^4]=exp[-16] GAMMA is used only when JPBC=1,2, or 3
					// skip line
	while(xmlStrcmp(nod->name,(const xmlChar *) "GAMMA")) nod = nod->next;
	nod = nod->children;
	ch = (char*) nod->content;
	nod = nod->parent;
	gamma = atof(ch);	
//
	printf("\n'**** Interaction cutoff parameter for PBC calculations ****'\n\n");
	if (jpbc > 0)
	{
		printf("Reapar %10.3lf = GAMMA = replica dipole summation limiter for PBC", gamma);
		if (gamma > 1.e-1 || gamma < 1.e-4)
		{
			printf("Reapar strange value of gamma ");
			return lineCounter;
		}
	}
	else
	{
		printf("Reapar [GAMMA is not used in present  (non-PBC) calculation]");
	}
//
// Define angular resolution used for calculating radiation force,
//     <cos>, <cos^2>, and radiation torque (if CMTORQ='DOTORQ')
//      ETASCA = parameter controlling number of scattering angles
//             = 1 is a good choice, but can reduce this if higher accuracy is required for angular averages
//      number of scattering angles will be approximately 6*[(3+x)/ETASCA]^2
	printf("\n'**** Angular resolution for calculation of <cos>, etc. ****'\n\n");

	while(xmlStrcmp(nod->name,(const xmlChar *) "ETASCA")) nod = nod->next;
	nod = nod->children;
	ch = (char*) nod->content;
	nod = nod->parent;
	etasca = atof(ch);

	if (etasca > 1.e-3)
	{
		printf("Reapar %10.3lf = ETASCA (parameter controlling number of scattering angles\n", etasca);
	}
	else
	{
		printf("Reapar %10.3lf  is not an appropriate value for ETASCA: check ddscat.par\n", etasca);
	}

//
// Define WAVELENGTH :
//	WAVINI	= first wavelength (physical units)
//	WAVEND	= last wavelength (physical units)
//	NWAV	= number of wavelengths
//	CDIVID	= 'LIN', 'LOG', or 'INV' for equal increments in wavelength, log(wavelength), or frequency
//			= 'TAB' to read wavelengths from file 'wave.tab'
	real wavini, wavend;
	printf("'**** Vacuum wavelengths (micron) ****'\n\n");

	while(xmlStrcmp(nod->name,(const xmlChar *) "VacuumWavelengths")) nod = nod->next;
	wavini = atof((char *) xmlGetProp(nod, (const xmlChar *)"First"));
	wavend = atof((char *) xmlGetProp(nod, (const xmlChar *)"Last"));
	nwav  = atof((char *) xmlGetProp(nod, (const xmlChar *)"HowMany"));
	char* cdivid = (char *) xmlGetProp(nod, (const xmlChar *)"How");

	if (strcmp(cdivid, "LIN") && strcmp(cdivid, "LOG") && strcmp(cdivid, "INV") && strcmp(cdivid, "TAB"))
	{
		printf("Reapar CDIVID for wavelengths must be LIN,LOG,INV, or TAB\n");
		return lineCounter;
	}
	if (strcmp(cdivid, "TAB"))
	{
		wavea = (real *)malloc(nwav * sizeof(real));
		printf("Reapar %3d wavelengths from  %7.4lf to %7.4lf\n", nwav, wavini, wavend);
		Divide(cdivid, wavini, wavend, nwav, wavea);
	}
	else
	{
		char Buffer[255];
		FILE *file28 = fopen("wave.tab", "r");
		dgets(Buffer, 255, file28);

		int curSize = 16;
		int deltaSize = 16;
		wavea = (real *)malloc(curSize * sizeof(real));

		nwav = 0;
		while(1)
		{
			char *ia = dgets(Buffer, 255, file28);
			if (!ia)
				break;
			sscanf(Buffer, "%lf", &wavini);
			wavea[nwav++] = wavini;
			if (nwav >= curSize)
			{
				wavea = (real *)realloc(wavea, (curSize + deltaSize) * sizeof(real));
				curSize += deltaSize;
			}
		}
		fclose(file28);
		if (curSize != nwav)
			wavea = (real *)realloc(wavea, nwav * sizeof(real));

		printf("Reapar %d  wavelengths from %lf  to %lf", nwav, wavea[0], wavea[nwav]);
	}
//
// Define NAMBIENT = refractive index of ambient medium
	real aefini, aefend; 
	printf("'**** Refractive index of ambient medium'\n\n");
//
	while(xmlStrcmp(nod->name,(const xmlChar *) "NAMBIENT")) nod = nod->next;
	nod = nod->children;
	ch = (char*) nod->content;
	nod = nod->parent;
	nambient = atof(ch);
	
	printf("Reapar %7.4f = NAMBIENT = refractive index of ambient medium\n", nambient);
//
// Define	AEFF	= effective radius a_eff 
//					= radius of equal-solid-volume sphere (physical units)
//			AEFINI	= first a_eff (physical units)
//			AEFEND	= last a_eff (physical units)
//			NRAD	= number of a_eff
//			CDIVID	= 'LIN', 'LOG', or 'INV' for equal increments in a_eff, log(a_eff), or 1/a_eff
//					= 'TAB' to read aeff from file 'aeff.tab'
	printf("'**** Effective Radii (micron) **** '\n\n");

	while(xmlStrcmp(nod->name, (const xmlChar *) "Aeff")) nod = nod->next;
	aefini = atof((char *) xmlGetProp(nod, (const xmlChar *)"First"));
	aefend = atof((char *) xmlGetProp(nod, (const xmlChar *)"Last"));
	nrad  = atof((char *) xmlGetProp(nod, (const xmlChar *)"HowMany"));
	aeffa = (real *)malloc(nrad * sizeof(real));
	cdivid = (char *) xmlGetProp(nod, (const xmlChar *)"How");
	cdivid[3] = '\0';

	if (strcmp(cdivid, "LIN") && strcmp(cdivid, "LOG") && strcmp(cdivid, "INV") && strcmp(cdivid, "TAB"))
	{
		printf("Reapar CDIVID for wavelengths must be LIN,LOG,INV, or TAB\n");
		return lineCounter;
	}
	if (strcmp(cdivid, "TAB"))
	{
		printf("Reapar %d eff. radii from %lf to %lf\n", nrad, aefini, aefend);
		Divide(cdivid, aefini, aefend, nrad, aeffa);
	}
	else
	{
		char Buffer[255];
		FILE *file28 = fopen("aeff.tab", "r");
		dgets(Buffer, 255, file28);					// skip one header line

		nrad = 0;
		while(1)
		{
			char *ia = dgets(Buffer, 255, file28);
			if (!ia)
				break;
			nrad++;
			sscanf(Buffer, realFormat, &aefini);
			aeffa[nwav] = aefini;
		}
		fclose(file28);
		printf("Reapar %3d' eff. radii from %7.4lf to %7.4lf", nrad, aeffa[0], aeffa[nrad-1]);
	}
//
// Define incident polarizations (in Lab Frame)
// It is assumed that incident radiation is along x axis in Lab Frame
	printf("'**** Define Incident Polarizations ****'\n\n");
//
// Read Complex polarization vector CXE01_LF=e01 (normalize if necessary).
	while(xmlStrcmp(nod->name,(const xmlChar *) "IncidentPolarization")) nod = nod->next;
	nod = nod->children;
	while(xmlStrcmp(nod->name,(const xmlChar *) "PolarizationState")) nod = nod->next;
	nod = nod->children;

	while(xmlStrcmp(nod->name,(const xmlChar *) "X")) nod = nod->next;
	cxe01_lf.Data(0).Re() = atof((char *) xmlGetProp(nod, (const xmlChar *)"Re"));
	cxe01_lf.Data(0).Im() = atof((char *) xmlGetProp(nod, (const xmlChar *)"Im"));

	while(xmlStrcmp(nod->name,(const xmlChar *) "Y")) nod = nod->next;
	cxe01_lf.Data(1).Re() = atof((char *) xmlGetProp(nod, (const xmlChar *)"Re"));
	cxe01_lf.Data(1).Im() = atof((char *) xmlGetProp(nod, (const xmlChar *)"Im"));

	while(xmlStrcmp(nod->name,(const xmlChar *) "Z")) nod = nod->next;
	cxe01_lf.Data(2).Re() = atof((char *) xmlGetProp(nod, (const xmlChar *)"Re"));
	cxe01_lf.Data(2).Im() = atof((char *) xmlGetProp(nod, (const xmlChar *)"Im"));

	real e1 = cxe01_lf.Data(0).mod();
	if (e1 != zero_)
	{
		printf("Reapar cxe01_lf(1) must be zero!\n");
		return lineCounter;
	}
// Normalize:
	e1 = sqrt(cxe01_lf.Data(1).modSquared() + cxe01_lf.Data(2).modSquared());
	cxe01_lf.Data(1) /= e1;
	cxe01_lf.Data(2) /= e1;
// Construct orthogonal normalized polarization vector CXE02=e02 using xhat cross e01* = e02
	cxe02_lf.Data(0).clear();
	cxe02_lf.Data(1) = -(cxe01_lf.Data(2).conjg());
	cxe02_lf.Data(2) =   cxe01_lf.Data(1).conjg();
// IORTH	= 1 to calculate only for single polarization
//			= 2 to also calculate for orthogonal polarization
	nod = nod->parent;
	nod = nod->parent;
	
	iorth = atoi((char *) xmlGetProp(nod, (const xmlChar *)"IORTH"));
	printf("Reapar IORTH=%2d\n", iorth);
//
// Specify whether or not to write ".sca" files  skip line:
	printf("'**** Specify which output files to write ****'\n\n");
//
// IWRKSC	= 0 to NOT write ".sca" file for each target orientation
//			= 1 to write ".sca" file for each target orientation
	while(xmlStrcmp(nod->name,(const xmlChar *) "IWRKSC")) nod = nod->next;
	nod = nod->children;
	ch = (char*) nod->content;
	nod = nod->parent;
	iwrksc = atoi(ch);
	printf("Reapar IWRKSC=%2d\n", iwrksc);
//
// Specify whether or not to write ".pol" files
// 12.01.23 (BTD) * eliminate this -- no longer needed since we now include support for fast nearfield calculations, and program call readE to read the nearfield output files.
//                * now initialize IWRPOL=0, although this can then be over-ridden to set IWRPOL=1 when NRFLD=1
//
// IWRPOL	= 0 to NOT write ".pol1" and ".pol2" file for each target orientation (beta,theta)
//			= 1 to write ".pol1" and ".pol2" file for each target orientation (beta,theta)
//
//	strcpy(cwhere, "error reading IWRPOL in ddscat.par");
//	fgets(Buffer, 255, iopar);
//	sscanf(Buffer, "%d", &iwrpol);
//	sprintf(cmsgnm, "IWRPOL=%2d", iwrpol);
//	Wrimsg(ReaparLabel, cmsgnm);
	iwrpol = 0;
//
// In the event that user set NRFLD=1 but IWRPOL=0, set IWRPOL to 1
	if(nrfld == 1 && iwrpol <= 0)
	{
		iwrpol = 1;
		printf("Reapar set iwrpol=1 because nrfld=1\n");
	}
//
// Read information determining target rotations
	printf("'**** Prescribe Target Rotations ****'\n\n");					
//
	while(xmlStrcmp(nod->name,(const xmlChar *) "PrescribeTargetRotations")) nod = nod->next;
	nod = nod->children;

	while(xmlStrcmp(nod->name,(const xmlChar *) "BETA")) nod = nod->next;
	betami = atof((char *) xmlGetProp(nod, (const xmlChar *)"BETAMI"));
	betamx = atof((char *) xmlGetProp(nod, (const xmlChar *)"BETAMX"));
	nbeta = atoi((char *) xmlGetProp(nod, (const xmlChar *)"NBETA"));

	while(xmlStrcmp(nod->name,(const xmlChar *) "THETA")) nod = nod->next;
	thetmi = atof((char *) xmlGetProp(nod, (const xmlChar *)"THETMI"));
	thetmx = atof((char *) xmlGetProp(nod, (const xmlChar *)"THETMX"));
	ntheta = atoi((char *) xmlGetProp(nod, (const xmlChar *)"NTHETA"));

	while(xmlStrcmp(nod->name,(const xmlChar *) "PHI")) nod = nod->next;
	phimin = atof((char *) xmlGetProp(nod, (const xmlChar *)"PHIMIN"));
	phimax = atof((char *) xmlGetProp(nod, (const xmlChar *)"PHIMAX"));
	nphi = atoi((char *) xmlGetProp(nod, (const xmlChar *)"NPHI"));
//
// check that user has not requested more than 1000 orientations and IWRKSC=1
	if (iwrksc > 0 && (nbeta*ntheta*nphi) > 1000)
	{
		printf("Reapar error: if iwrksc=1, nbeta*ntheta*nphi must be .le. 1000\n");
		return lineCounter;
	}
	printf("'**** Specify first IWAV, IRAD, IORI (normally 0 0 0) ****'\n\n");
	
	nod = nod->parent;
	while(xmlStrcmp(nod->name,(const xmlChar *) "SpecifyFirst")) nod = nod->next;

	iwav0 = atoi((char *) xmlGetProp(nod, (const xmlChar *)"IWAV"));
	irad0 = atoi((char *) xmlGetProp(nod, (const xmlChar *)"IRAD"));
	iori0 = atoi((char *) xmlGetProp(nod, (const xmlChar *)"IORI"));
//
// Check that IWAV0,IRAD0,IORI0 are OK:
	if (iwav0+1 > nwav)
	{
		printf("Reapar iwav0+1 > nwav");
		return lineCounter;
	}
	if (irad0+1 > nrad) 
	{
		printf("Reapar irad0+1 > nrad");
		return lineCounter;
	}
	if (iori0+1 > nbeta*ntheta*nphi) 
	{
		printf("Reapar iori0+1 > nbeta*ntheta*nphi");
		return lineCounter;
	}
//
// If NPHI>1, then set IORTH=2 regardless of value input.
	if ((iorth == 1) && (nphi > 1))
	{
		iorth = 2;
		printf("Reapar set iorth=2 since nphi>1 ");
	}	
	switch(iorth)
	{
	case 1:
		printf("Reapar Calculate only for single polarization \n");
		break;

	case 2:
		printf("Reapar Do orthogonal polarization for each case \n");
		break;

	default:
		printf("Reapar WRONG VALUE OF IORTH \n");
		return lineCounter;
		break;
	}
	printf("Reapar %7.2lf%7.2lf  Range of BETA values ; NBETA =%d\n", betami, betamx, nbeta);
	printf("Reapar %7.2lf%7.2lf  Range of THETA values; NTHETA=%d\n", thetmi, thetmx, ntheta);
	printf("Reapar %7.2lf%7.2lf  Range of PHI values ;   NPHI =%d\n", phimin, phimax, nphi);
//
// Convert from degrees to radians
	betami /= Degrad;
	betamx /= Degrad;
	thetmi /= Degrad;
	thetmx /= Degrad;
	phimin /= Degrad;
	phimax /= Degrad;
//
// Specify elements of scattering matrix to be printed out
	printf("'**** Select Elements of S_ij Matrix to Print ****'\n\n");
	while(xmlStrcmp(nod->name,(const xmlChar *) "S_ijMatrix")) nod = nod->next;
	nsmelts = atoi((char *) xmlGetProp(nod, (const xmlChar *)"NSMELTS"));

	if (nsmelts > 9)
	{
		printf("Reapar Error: %d > 9 is not permitted by output routines\n", nsmelts);
		printf("Reapar NSMELTS > 9 is not allowed in ddscat.par\n");
		return lineCounter;
	}
	if (nsmelts <= 0)
	{
		smind1[0] = 11;
		smind1[1] = 21;
		smind1[2] = 31;
		smind1[3] = 41;
		smind1[4] = 12;
		smind1[5] = 13;
		nsmelts = 6;
		//dgets(Buffer, 255, iopar);
	}
	else
	{
		nod = nod->children;
		for(j=0; j<nsmelts; ++j)
		{
			nod = nod->next;
			while(xmlStrcmp(nod->name,(const xmlChar *) "ij")) nod = nod->next;
			nod = nod->children;
			ch = (char*) nod->content;
			nod = nod->parent;
			smind1[j] = atoi(ch);
		}
	}
	for(j=0; j<nsmelts; ++j)
	{
		smind2[j] = smind1[j] % 10;
		smind1[j] = smind1[j] / 10;
	}
	nod = nod->parent;
//
// Specify scattering directions to be calculated
// Two options:
// CMDFRM	= 'LFRAME': specify scattering directions n in Lab Frame (where incident beam is in x-direction)
//			  THETAN,PHIN = Direction of n from n0 in calculation of the scattering matrix; cos(THETAN) is n0 \dot n
//			  and PHIN is azimuthal angle of n from Lab xy pla 
// Note: THETA1, THETA2, and PHI for each scatterin plane are entered in degrees, and immediately converted to radians.
// Arbitrary number of scattering planes may be considered.
//
// CMDFRM	= 'TFRAME': specify scattering directions in Target Frame, defined by target axes a1,a2,a3
// If JPBC = 0:
//			THETAN,PHIN = Direction of n relative to a1,a2,a3: 
//			THETAN = angle between n and a1
//			PHIN   = angle between a1-n plane and a1-a2 plane
//    JPBC = 1:
//			THETAN = Diffraction order along y_TF axis
//			PHIN   = azimuthal angle around y_TF
//	  JPBC = 2:
//			THETAN = Diffraction order along z_TF axis
//			PHIN   = azimuthal angle around z_TF
//	  JPBC = 3:
//			THETAN = Diffraction order along y_TF
//			PHIN   = Diffraction order along z_TF
//where we first run through transmitted waves 1 -> NSCAT/2 
//      and then run through reflected waves NSCAT/2+1 -> NSCAT
// 
// Three cases:
//	JPBC = 0: single isolated target.
//				specify
//				phi for scattering plane
//				thetamin, thetamax, dtheta for scattering plane
//	JPBC = 1,2: periodic in one dimension.
//				specify
//				diffraction order in direction of target periodicity
//				phimin, phimax, dphi for scattering cone
//	JPBC=3:		periodic in two dimensions
//				specify
//				difraction order in y direction and order in z direction
	real phi1, theta1, theta2, dtheta, delta;
	int nplanes, jplane, nsca0;
	nsca = 0;
// 
	while(xmlStrcmp(nod->name,(const xmlChar *) "ScatteredDirections")) nod = nod->next;
	strncpy(cmdfrm, (char *) xmlGetProp(nod, (const xmlChar *)"CMDFRM"), 6);
	cmdfrm[6] = '\0';
	strupr(cmdfrm);
	//jpbc = 0;
	if (strcmp(cmdfrm, "LFRAME") && strcmp(cmdfrm, "TFRAME"))
	{
		printf("Reapar Error reading ddscat.par file\n");
		return lineCounter;
	}
	if (jpbc > 0)
	{
		if (!strcmp(cmdfrm, "LFRAME"))
		{
			printf("Reapar cannot use LFRAME when JPBC=>0\n");
			return lineCounter;
		}
	}

	if (jpbc == 0)
	{
		if (!strcmp(cmdfrm, "LFRAME"))
		{
			printf("Reapar cmdfrm=%s : scattering directions given in lab frame\n", cmdfrm);
		}
		else
		{
			printf("Reapar cmdfrm=%s : scattering directions given in target frame\n", cmdfrm);
		}
	}

	nplanes = atoi((char *) xmlGetProp(nod, (const xmlChar *)"NPLANES"));

	switch(jpbc)
	{
	case 0:
		printf("Reapar %4d = number of scattering planes\n", nplanes);
		break;

	case 1:
	case 2:
		printf("Reapar %4d = number of scattering cones\n", nplanes);
		break;

	case 3:
		printf("Reapar %4d = number of diffraction orders for transmission\n", nplanes);
		break;

	default:
		break;
	}

	if (nplanes > 0)
	{
		nod = nod->children;
		if (jpbc < 3)
		{
			for(jplane=0; jplane<nplanes; ++jplane)
			{
				nod = nod->next;
				while(xmlStrcmp(nod->name,(const xmlChar *) "Plane")) nod = nod->next;
				phi1 = atof((char *) xmlGetProp(nod, (const xmlChar *)"phi"));
				theta1 = atof((char *) xmlGetProp(nod, (const xmlChar *)"thetan_min"));
				theta2 = atof((char *) xmlGetProp(nod, (const xmlChar *)"thetan_max"));
				dtheta = atof((char *) xmlGetProp(nod, (const xmlChar *)"dtheta"));

				if (jpbc < 1)
					printf("Reapar %7.1lf%7.1lf%7.1f = phi, theta_min, theta_max for scatt. plane%4d\n", phi1, theta1, theta2, jplane);
				else
					printf("Reapar %7.1lf%7.1lf%7.1f = order, zeta_min, zeta_max for scattering cone\n", phi1, theta1, theta2);
				if (jpbc == 0)
				{
					if(Fabs(theta1-theta2) > zero_ && dtheta == zero_)
					{
						printf("Reapar DTHETA=0 in ddscat.par!\n");
						return lineCounter;
					}
				}
// ! Convert to radians
				phi1 /= Degrad;
				theta1 /= Degrad;
				theta2 /= Degrad;
				dtheta /= Degrad;
// ! Allow for possibility that user did not enter correct sign for
// ! DTHETA (if I did it, others will too...)
				if (theta2 < theta1)
					dtheta = -Fabs(dtheta);
// ! compute theta values for this scattering plane/cone
				nsca0 = 1;
				delta = theta2 - theta1;
				if (delta)
					nsca0 = 1 + nint_(delta/dtheta);
// ChB: ! realloc to have sufficient space
				thetan = (real *)realloc(thetan, (nsca+nsca0)*sizeof(real));
				phin = (real *)realloc(phin, (nsca+nsca0)*sizeof(real));
				thetan[nsca] = theta1;
				if(nsca0 > 2)
				{
					for(j=1; j<nsca0-1; ++j)
					{
						thetan[nsca+j] = theta1 + j*dtheta;
					}
				}
				thetan[nsca+nsca0-1] = theta2;
				for(j=0; j<nsca0; ++j)
				{
					phin[nsca+j] = phi1;
				}
				nsca += nsca0;
				printf("Reapar %4d = number of scattering angles in this scattering %s\n", nsca0, (jpbc < 1 ? "plane" : "cone"));
			}
		}
		else
		{
			sprintf(curFormat, "%s%s", realFormat, realFormat);
			for(jplane=0; jplane<nplanes; ++jplane)
			{
				fgets(Buffer, 255, iopar);
				sscanf(Buffer, curFormat, phin[nsca], thetan[nsca]);
				++nsca;
			}
			for(j=0; j<nsca; ++j)
			{
				phin[j+nsca] = phin[j];
				thetan[j+nsca] = thetan[j];
			}
			nsca = 2*nsca;
		}
	}


	//if(xmlStrcmp(nod->name,(const xmlChar *)"request")) return -1;
	//nod = nod->xmlChildrenNode;
	//xmlErrorPtr err =  xmlGetLastError();
	//ShowPar();
	// cxe01_lf; cxe02_lf;	cfleps;
	xmlUnlinkNode(nod);
	xmlFreeNode(nod);
	xmlFreeDoc(doc);
	system("pause");
	return 0;
}

FOR DDscatParameters.h

void TakeSP(int numSP, xmlNodePtr &nod1)
{
	char*ch;
	for (int i = 0; i < numSP ; i++)
	{
		nod1 = nod1->next;
		while(xmlStrcmp(nod1->name,(const xmlChar *) "SP")) nod1 = nod1->next;
		ch = (char*) xmlGetProp(nod1, nod1->properties->name);
		shpar[i] = atof(ch);
		printf("Reapar %f = sphar[%d]\n",shpar[i] , i);
	}
}
** */
