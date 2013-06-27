#include "StdAfx.h"

#include "ProcessHelper.h"

ProcessHelper::ProcessHelper(void)
{
	cxbsca = cxeps = cxpol = cxesca = cxadia = cxeinc = cxbinc = NULL;
	s = NULL;
}

ProcessHelper::~ProcessHelper(void)
{
	CleanDelete2(cxbsca);
	CleanDelete2(cxeps);
	CleanDelete2(cxpol);
	CleanDelete2(cxesca);
	CleanDelete2(cxadia);
	CleanDelete2(cxeinc);
	CleanDelete2(cxbinc);
	CleanDelete2(s);
}

bool ProcessHelper::Load(const char *cflename, int nrword)
{
#ifdef _WIN32
	int mode = O_BINARY|O_RDONLY;
#else
	int mode = O_RDONLY;
#endif	
	int binaryFile = open(cflename, mode);
	if (binaryFile == -1)
		return false;
//
	char Buffer[256];
	int i = -1;
	do
	{
		++i;
		read(binaryFile, Buffer+i, sizeof(char));
	} while(Buffer[i]);
	stringVersion = string(Buffer);
	read(binaryFile, &intVersion, sizeof(int));
//
	read(binaryFile, &nrword_nf, sizeof(int));
	if (nrword_nf != nrword)
	{
        fprintf(stderr, "word length=%d in file %s\n", nrword_nf, cflename);
		fprintf(stderr, "word length=%d in program readnf\n", nrword);
		close(binaryFile);
		return false;
	}
//
// Start loading the field data
	read(binaryFile, &i, sizeof(int));
	nrfldb = (NearfieldBMethod)i;
	read(binaryFile, &nxyz, sizeof(int));
	read(binaryFile, &nat0, sizeof(int));
	read(binaryFile, &nat3, sizeof(int));
	read(binaryFile, &ncomp, sizeof(int));
	read(binaryFile, &nx, sizeof(int));
	read(binaryFile, &ny, sizeof(int));
	read(binaryFile, &nz, sizeof(int));
	x0.Read(binaryFile);
	read(binaryFile, &aeff, sizeof(real));
	read(binaryFile, &nambient, sizeof(real));
	read(binaryFile, &wave, sizeof(real));
	ak_tf.Read(binaryFile);
	cxe0_tf.Read(binaryFile);
	cxb0_tf.Read(binaryFile);
//
	cxeps = new Complex[ncomp];
	read(binaryFile, cxeps, ncomp * sizeof(Complex));
	icomp.Dimension(nxyz, 3);
	for(i=0; i<nxyz; ++i)
		icomp.ReadLine(binaryFile, i);
	cxpol = new Complex[nat3];
	read(binaryFile, cxpol, nat3 * sizeof(Complex));
	cxesca = new Complex[nat3];
	read(binaryFile, cxesca, nat3 * sizeof(Complex));
	cxadia = new Complex[nat3];
	read(binaryFile, cxadia, nat3 * sizeof(Complex));
	cxeinc = new Complex[nat3];
//
	if(nrfldb == NearfieldBMethodBoth)
	{
		cxbsca = new Complex[nat3];
		cxbinc = new Complex[nat3];
		read(binaryFile, cxbsca, nat3 * sizeof(Complex));
	}
	close(binaryFile);
//
// compute phase for JX=JY=JZ=0
	real phi0 = ak_tf.Scalar(x0);
	for(int jx=0; jx<nx; ++jx)
	{
		real phix = (jx + 1) * ak_tf.data[0] + phi0;
		for(int jy=0; jy<ny; ++jy)
		{
			real phiy = phix + (jy + 1) * ak_tf.data[1];
			for(int jz=0; jz<nz; ++jz)
			{
				real arg = phiy + (jz + 1) * ak_tf.data[2];
				Complex cxphas = Complex(Cos(arg), Sin(arg));
				int index = (jx * ny + jy) * nz + jz;
// compute E_macro contributed by incident wave
				for(int k=0; k<3; ++k)
				{
					short ic = icomp.Value(index, k);
					Complex cxfac((real)1., (real)0.);
					if ((ic > 0) && (ic <= ncomp))
						cxfac = Complex((real)3., (real)0.) / (cxeps[ic-1] + (real)2.);
					cxeinc[3*index + k] = cxfac * cxe0_tf.data[k] * cxphas;
				}
				if(nrfldb == NearfieldBMethodBoth)
				{
					for(int k=0; k<3; ++k)
						cxbinc[3*index + k] = cxb0_tf.data[k] * cxphas;
				}
			}
		}
	}

	return true;
}

void ProcessHelper::PrepareNormsAndPointing()					// VTK supplementary arrays
{
	s_inc.Set((cxe0_tf.data[1] * cxe0_tf.data[2].conjg() - cxe0_tf.data[2] * cxe0_tf.data[1].conjg()).re, 
			  (cxe0_tf.data[2] * cxe0_tf.data[0].conjg() - cxe0_tf.data[0] * cxe0_tf.data[2].conjg()).re,
			  (cxe0_tf.data[0] * cxe0_tf.data[1].conjg() - cxe0_tf.data[1] * cxe0_tf.data[0].conjg()).re);
	snorm = s_inc.ModSquared();
//
// determine E_inc, E_sca, and P at points along defined track
	const real dd = (real)0.5001;
	dphys = aeff * Pow(((real)4. * Pi / nat0 / (real)3.), (real)1./(real)3.);
	xmin = (real)((x0.data[0] + (real)1. - dd) * dphys); 
	xmax = (real)((x0.data[0] +       nx + dd) * dphys);
	ymin = (real)((x0.data[1] + (real)1. - dd) * dphys);
	ymax = (real)((x0.data[1] +       ny + dd) * dphys);
	zmin = (real)((x0.data[2] + (real)1. - dd) * dphys);
	zmax = (real)((x0.data[2] +       nz + dd) * dphys);
//
// calculate normalized Poynting vector at all points
	if (nrfldb != NearfieldBMethodBoth)
		return;
	s = new real[nat3];
	for(int jx=0; jx<nx; ++jx)
	{
		for(int jy=0; jy<ny; ++jy)
		{
			for(int jz=0; jz<nz; ++jz)
			{
				int index = 3 * ((jx * ny + jy) * nz + jz);
				Complex cxex = cxeinc[index    ] + cxesca[index    ];
				Complex cxey = cxeinc[index + 1] + cxesca[index + 1];
				Complex cxez = cxeinc[index + 2] + cxesca[index + 2];
				Complex cxbx = cxbinc[index    ] + cxbsca[index    ];
				Complex cxby = cxbinc[index + 1] + cxbsca[index + 1];
				Complex cxbz = cxbinc[index + 2] + cxbsca[index + 2];
				s[index    ] = (cxey * cxbz.conjg() - cxez * cxby.conjg()).re / snorm;
				s[index + 1] = (cxez * cxbx.conjg() - cxex * cxbz.conjg()).re / snorm;
				s[index + 2] = (cxex * cxby.conjg() - cxey * cxbx.conjg()).re / snorm;
			}
		}
	}
}

int ProcessHelper::SanityCheck()
{
//
// Check solution
// this check neglects the off-diagonal elements of A
// this will not be a valid assumption for anisotropic materials
// that do not have optical axes aligned with TF xyz axes

    int j;
	real einc2 = (real)0.;
    for(j=0; j<3; ++j)
	{
		einc2 += ((cxe0_tf.data[j] * cxe0_tf.data[j].conjg())).re;
	}
//
	sumerr2 = (real)0.;
	int jj1 = 0;
	for(j=0; j<nxyz; ++j)
	{
		if(icomp.Value(j, 0) > 0)
		{
			for(int jx=0; jx<3; ++jx)
			{
				int index = 3*j + jx;
				int ii = icomp.Value(j, jx) - 1;
				if ((ii < ncomp) && (ii >= 0))
				{
					Complex cxerr = cxpol[index] * cxadia[index] - (cxeinc[index] + cxesca[index]) * (cxeps[ii] + (real)2.) / (real)3.;
					sumerr2 += cxerr.modSquared();
				}
			}
			++jj1;
		}
	}
	sumerr2 /= (real(nat0) * einc2);

	return jj1;
}

bool ProcessHelper::PrepareVTKFile(const char *cflvtr, int ivtr)
{
// write VTK file for graphics
// define mesh for graphics
	Vtr vtr;
	vtr.OpenFile(cflvtr);
//
	int j;
	real *vtrx = new real[nx];
	for(j=0; j<nx; ++j)
		vtrx[j] = (j + 1 + x0.data[0]) * dphys;
	real *vtry = new real[ny];
	for(j=0; j<ny; ++j)
		vtry[j] = (j + 1 + x0.data[1]) * dphys;
	real *vtrz = new real[nz];
	for(j=0; j<nz; ++j)
		vtrz[j] = (j + 1 + x0.data[2]) * dphys;
	vtr.WriteMesh3d(vtrx, vtry, vtrz, nx, ny, nz);
	delete [] vtrx;
	delete [] vtry;
	delete [] vtrz;

// mesh x,y,z assuming that we are on rectangular grid)
// variables are dimensioned nx,ny,nz
	Array3F<real> vtr8, vtrvx, vtrvy, vtrvz;
	vtr8.Dimension(nx, ny, nz);

	if ((ivtr == 3) && (nrfldb == NearfieldBMethodBoth))
	{
		vtrvx.Dimension(nx, ny, nz);
		vtrvy.Dimension(nx, ny, nz);
		vtrvz.Dimension(nx, ny, nz);
	}	
//
// electric field intensity or energy density
// fprintf(stdout, "Ivtr = %d\n", ivtr);
	int jx, jy, jz;
	j = 0;
	for(jx=0; jx<nx; ++jx)
	{
		for(jy=0; jy<ny; ++jy)
		{
			for(jz=0; jz<nz; ++jz)
			{
				Complex cxee[3];
				cxee[0] = cxeinc[j    ] + cxesca[j    ];
				cxee[1] = cxeinc[j + 1] + cxesca[j + 1];
				cxee[2] = cxeinc[j + 2] + cxesca[j + 2];
// intrinsic function dot_product(cx,cy)=sum(conjg(cx),cy)
//         vtr8(ix,iy,iz)=sqrt(sum(abs(cxee)**2))
				if (ivtr < 3)
				{
					real module = (cxee[0]*cxee[0].conjg() + cxee[1]*cxee[1].conjg() + cxee[2]*cxee[2].conjg()).re;
					vtr8.Value(jx, jy, jz) = (ivtr == 1) ? Sqrt(module) : module;
// if (icomp.Value(j/3, 0) != 0)
//   fprintf(stdout, "%3d%3d%3d%12.6f\n", jx, jy, jz, vtr8.Value(jx, jy, jz));
//fprintf(stderr, "%3d%3d%3d%3d%12.6f%12.6f%12.6f%12.6f%12.6f%12.6f%12.6f%12.6f%12.6f%12.6f%12.6f%12.6f\n", jx, jy, jz, icomp.Value(j/3, 0), 
//		cxeinc[j].re, cxeinc[j].im, cxeinc[j+1].re, cxeinc[j+1].im, cxeinc[j+2].re, cxeinc[j+2].im,
//		cxesca[j].re, cxesca[j].im, cxesca[j+1].re, cxesca[j+1].im, cxesca[j+2].re, cxesca[j+2].im);
				}
				else
				{
//					Complex cxbb[3];
//					cxbb[0] = cxbinc[j    ] + cxbsca[j    ];
//					cxbb[1] = cxbinc[j + 1] + cxbsca[j + 1];
//					cxbb[2] = cxbinc[j + 2] + cxbsca[j + 2];
					vtrvx.Value(jx, jy, jz) = s[j    ];
					vtrvy.Value(jx, jy, jz) = s[j + 1];
					vtrvz.Value(jx, jy, jz) = s[j + 2];
				}
				j += 3;
			}
		}
	}
// 
	switch(ivtr)
	{
	case 1:
        vtr.WriteScalar3d("Intensity", vtr8);
		break;

	case 2:
		vtr.WriteScalar3d("Esquared", vtr8);
		break;

	case 3:
		vtr.WriteVector3d("Pointing", vtrvx, vtrvy, vtrvz, nx, ny, nz);
		break;

	default:
		break;
	}
// Composition 
// This can be used to display nicely inhomogeneous objects ICOMP=0 outside the target
	j = 0;
	for(jx=0; jx<nx; ++jx)
	{
		for(jy=0; jy<ny; ++jy)
		{
			for(jz=0; jz<nz; ++jz)
			{
// let us output only one component of icomp
				vtr8.Value(jx, jy, jz) = (real)icomp.Value(j, 0);
				++j;
			}
		}
	}
	vtr.WriteScalar3d("Composition", vtr8);
	vtr.CloseFile();
	vtr.CollectFile();				// Produces "output.vtr" and "output.pvd" files
//
// We are done with 3d objects (large files)
// Now VTR output for smaller (3)-vectors defining geometry
// we probably should also output CXE0 (but it is complex, what to do?)
// copy to kind=8 precision for graphics and output

	Array3F<real> u, v, w;
	u.Dimension(1, 1, 1);
	v.Dimension(1, 1, 1);
	w.Dimension(1, 1, 1);
	real xvect = (real)0.;
	real yvect = (real)0.;
	real zvect = (real)0.;
	vtr.OpenFile("akr");
	vtr.WriteMesh3d(&xvect, &yvect, &zvect, 1, 1, 1);
	u.Value(0, 0, 0) = ak_tf.data[0];
	v.Value(0, 0, 0) = ak_tf.data[1];
	w.Value(0, 0, 0) = ak_tf.data[2];
	vtr.WriteVector3d("akr", u, v, w, 1, 1, 1);
	u.Value(0, 0, 0) = cxe0_tf.data[0].re;
	v.Value(0, 0, 0) = cxe0_tf.data[1].re;
	w.Value(0, 0, 0) = cxe0_tf.data[2].re;
	vtr.WriteVector3d("e0", u, v, w, 1, 1, 1);
	vtr.CloseFile();
	vtr.CollectFile();
//
	return true;
}

bool ProcessHelper::IsLiniaOk(const Linia &op)
{
	return ((op.GetData(0) >= xmin) && (op.GetData(3) <= xmax) && (op.GetData(1) >= ymin) && (op.GetData(4) <= ymax) && (op.GetData(2) >= zmin) && (op.GetData(5) <= zmax));
}

void ProcessHelper::OutMinMax(void)
{
	fprintf(stderr, "xmin,ymin,zmin = %13.2e%13.2e%13.2e\n", xmin, ymin, zmin);
	fprintf(stderr, "xmax,ymax,zmax = %13.2e%13.2e%13.2e\n", xmax, ymax, zmax);
}

bool ProcessHelper::WriteLiniaData(const char *filename, const Linia &linia)
{
	FILE *outputFile = fopen(filename, "w");
	if (!outputFile)
		return false;
//
	fprintf(outputFile, "%10.4e = a_eff (radius of equal-volume sphere)\n", aeff);
	fprintf(outputFile, "%10.4e = d = dipole spacing\n", dphys);
	fprintf(outputFile, "%10d = N = number of dipoles in target or TUC\n", nat0);
	fprintf(outputFile, "%8d%8d%8d = NX,NY,NZ = extent of computational volume\n", nx, ny, nz);
	fprintf(outputFile, "%10.5e = wavelength in vacuo\n", wave);
	fprintf(outputFile, "%10.5lf = refractive index of ambient medium\n", nambient);
	ak_tf.Fprintf(outputFile, "%10.5lf", NULL, " = k_{inc,TF} * d\n");
	fprintf(outputFile, "%10.5lf %10.5lf = Re(E_inc,x) Im(E_inc,x) at x_TF=0,y_TF=0,z_TF=0\n", cxe0_tf.data[0].re, cxe0_tf.data[0].im);
	fprintf(outputFile, "%10.5lf %10.5lf = Re(E_inc,y) Im(E_inc,y) \"\n", cxe0_tf.data[1].re,  cxe0_tf.data[1].im);
	fprintf(outputFile, "%10.5lf %10.5lf = Re(E_inc,z) Im(E_inc,z) \"\n", cxe0_tf.data[2].re,  cxe0_tf.data[2].im);
	if (nrfldb != NearfieldBMethodElectric)
	{
		fprintf(outputFile, "%10.5lf %10.5lf = Re(B_inc,x) Im(B_inc,x) at x_TF=0,y_TF=0,z_TF=0\n", cxb0_tf.data[0].re, cxb0_tf.data[0].im);
		fprintf(outputFile, "%10.5lf %10.5lf = Re(B_inc,y) Im(B_inc,y) \"\n", cxb0_tf.data[1].re,  cxb0_tf.data[1].im);
		fprintf(outputFile, "%10.5lf %10.5lf = Re(B_inc,z) Im(B_inc,z) \"\n", cxb0_tf.data[2].re,  cxb0_tf.data[2].im);
	}
	fprintf(outputFile, "%10.5lf %10.5lf %10.5lf  = 2*(4pi/c)*<S_inc> where <S_inc>=time-averaged incident Poynting vector                       [Poynting vector S = (Sx,Sy,Sz) = [Poynting vector S = (Sx,Sy,Sz) =  (c/4pi)*Re(E)xRe(B) ]\n", s_inc.data[0], s_inc.data[1], s_inc.data[2]);
	fprintf(outputFile, "%10.5lf = 2*(4pi/c)*|<S_inc>|\n", snorm);
	if (nrfldb == NearfieldBMethodElectric)
	{
		fprintf(outputFile, "                                 ----------------------- E field ---------------------------\n");
		fprintf(outputFile, "   x_TF       y_TF       z_TF      Re(E_x)   Im(E_x)   Re(E_y)   Im(E_y)   Re(E_z)   Im(E_z)\n");
	}
	else
	{
		fprintf(outputFile, "                                 ----------------------- E field ---------------------------  ----------------------- B field --------------------------  -normalized Poynting vector-\n");
		fprintf(outputFile, "   x_TF       y_TF       z_TF      Re(E_x)   Im(E_x)   Re(E_y)   Im(E_y)   Re(E_z)   Im(E_z)   Re(B_x)   Im(B_x)   Re(B_y)   Im(B_y)   Re(B_z)   Im(B_z)  -- <(Sx,Sy,Sz)>/|<S_inc>| --\n");
	}
//
	const real onex_ = (real)1.;
	Complex cxe_inc[3], cxe_sca[3], cxp[3], cxe[3], cxb_inc[3], cxb_sca[3], cxb[3];
	real xtf[3], w[8], svec[3];
	int jj[8], j;
	int nzy = nz * ny;
	for(int ja=0; ja<linia.GetNab(); ++ja)
	{
		real zeta = (real)ja / (real)linia.GetNab();
		linia.ParametricPoint(zeta, xtf);
//
// XTF is assumed to be in physical units
// XTF/DPHYS is in dipole units
// I + X0 = XTF/DPHYS
// I = XTF/DPHYS - X0

		int ix1 = (int)(xtf[0] / dphys - x0.data[0]);
		int iy1 = (int)(xtf[1] / dphys - x0.data[1]);
		int iz1 = (int)(xtf[2] / dphys - x0.data[2]);
		if (ix1 < 1) ix1 = 1;
		if (ix1 >= nx) ix1 = nx - 1;
		if (iy1 < 1) iy1 = 1;
		if (iy1 >= ny) iy1 = ny - 1;
		if (iz1 < 1) iz1 = 1;
		if (iz1 >= nz) iz1 = nz - 1;
//
// -------- determine indices ----------------- zyx
		jj[0] = (iz1 - 1) + nz * (iy1 - 1) + nzy * (ix1 - 1);		// 000
		jj[1] = jj[0] + nzy;										// 001
		jj[2] = jj[0] + nz;											// 010
		jj[3] = jj[2] + nzy;										// 011
		jj[4] = jj[0] + 1;											// 100
		jj[5] = jj[4] + nzy;										// 101
		jj[6] = jj[4] + nz;											// 110
		jj[7] = jj[6] + nzy;										// 111
		real wx = xtf[0] / dphys - x0.data[0] - ix1;
		real wy = xtf[1] / dphys - x0.data[1] - iy1;
		real wz = xtf[2] / dphys - x0.data[2] - iz1;
//
// ------------ determine weights ---------- zyx
		w[0] = (onex_ - wx) * (onex_ - wy) * (onex_ - wz);			// 000
		w[1] = wx * (onex_ - wy) * (onex_ - wz);					// 001
		w[2] = (onex_ - wx) * wy * (onex_ - wz);					// 010
		w[3] = wx * wy * (onex_ - wz);								// 011
		w[4] = (onex_ - wx) * (onex_ - wy) * wz;					// 100
		w[5] = wx * (onex_ - wy) * wz;								// 101
		w[6] = (onex_ - wx) * wy * wz;								// 110
		w[7] = wx * wy * wz;										// 111
//
		switch(nrfldb)
		{
			case 0:
// -------- evaluate weighted averages -----------------------------------
				for(j=0; j<3; ++j)
				{
					cxe_inc[j].clear();
					cxe_sca[j].clear();
					cxp[j].clear();
					for(int ji=0; ji<8; ++ji)
					{
						int index = 3 * jj[ji];
						cxe_inc[j] += (cxeinc[index + j] * w[ji]);
						cxe_sca[j] += (cxesca[index + j] * w[ji]);
						cxp[j]     += (cxpol[index + j] * w[ji]);
					}
					cxe[j] = cxe_inc[j] + cxe_sca[j];
				}
//
// write out total E field at point on track
				fprintf(outputFile, "%11.3e%11.3e%11.3e%10.5lf%10.5lf%10.5lf%10.5lf%10.5lf%10.5lf\n", 
					xtf[0], xtf[1], xtf[2], cxe[0].re, cxe[0].im, cxe[1].re, cxe[1].im, cxe[2].re, cxe[2].im);
				break;

			case 1:
// -------- evaluate weighted averages -----------------------------------
				for(j=0; j<3; ++j)
				{
					cxb_inc[j].clear();
					cxb_sca[j].clear();
					for(int ji=0; ji<8; ++ji)
					{
						int index = 3 * jj[ji];
						cxb_inc[j] += (cxbinc[index + j] * w[ji]);
						cxb_sca[j] += (cxbsca[index + j] * w[ji]);
						svec[j] += (s[index + j] * w[ji]);
					}
					cxb[j] = cxb_inc[j] + cxb_sca[j];
				}
//
// write out total E and B fields and Pointing vector at point on track
				fprintf(outputFile, "%11.3e%11.3e%11.3e%10.5lf%10.5lf%10.5lf%10.5lf%10.5lf%10.5lf%10.5lf%10.5lf%10.5lf%10.5lf%10.5lf%10.5lf%10.5lf%10.5lf%10.5lf\n", 
					xtf[0], xtf[1], xtf[2], cxe[0].re, cxe[0].im, cxe[1].re, cxe[1].im, cxe[2].re, cxe[2].im, 
					cxb[0].re, cxb[0].im, cxb[1].re, cxb[1].im, cxb[2].re, cxb[2].im,
					svec[0], svec[1], svec[2]);
				break;

				default:
					break;
		}

	}

	return true;
}
