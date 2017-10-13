#include "cpml.h"
#include "logger.h"
#include "forwardmodeling.h"
#define SIGN(a,b) (b >= 0 ? fabs(a) : -fabs(a))
CPML::CPML() {
	init = false;
	count = 0;
}

void CPML::GetXBoundaryMPos(int xPos, int zPos, int *xMPos, int *zMPos, const int nx, const int nz, const ForwardModeling &fm) const 
{
	int d = fm.getFDLEN();
	int bx0 = fm.getbx0();
	int bxn = fm.getbxn();
	int nb = bx0 - d;
	int cnx = nx - bx0 - bxn;
	if(xPos - d < nb)
		*xMPos = xPos - d;
	else
		*xMPos = xPos - d - cnx;
	*zMPos = zPos - d;
}

void CPML::GetZBoundaryMPos(int xPos, int zPos, int *xMPos, int *zMPos, const int nx, const int nz, const ForwardModeling &fm) const
{
	int d = fm.getFDLEN();
	int bx0 = fm.getbx0();
	int bz0 = fm.getbz0();
	int bzn = fm.getbzn();
	int nb = bx0 - d;
	int cnz = nz - bz0 - bzn;
	*xMPos = xPos - d;
	if(zPos - d < nb)
		*zMPos = zPos - d;
	else
		*zMPos = zPos - d - cnz;
}

void CPML::initCPML(const int nx, const int nz, const ForwardModeling &fm) {
	int d = fm.getFDLEN();
	int bx0 = fm.getbx0();
	int nb = bx0 - d;

	psiX.resize(2 * nb * (nz - d * 2));
	psiXLa.resize(2 * nb * (nz - d * 2));
	phiX.resize(2 * nb * (nz - d * 2));
	phiXLa.resize(2 * nb * (nz - d * 2));
	EtaX.resize(2 * nb * (nz - d * 2));
	EtaXLa.resize(2 * nb * (nz - d * 2));
	psi2X.resize(2 * nb * (nz - d * 2));
	psi2XLa.resize(2 * nb * (nz - d * 2));
	phi2X.resize(2 * nb * (nz - d * 2));
	phi2XLa.resize(2 * nb * (nz - d * 2));
	u020BXLa.resize(2 * nb * (nz - d * 2));

	psiX.assign(2 * nb * (nz - d * 2),0);
	psiXLa.assign(2 * nb * (nz - d * 2),0);
	phiX.assign(2 * nb * (nz - d * 2),0);
	phiXLa.assign(2 * nb * (nz - d * 2),0);
	EtaX.assign(2 * nb * (nz - d * 2),0);
	EtaXLa.assign(2 * nb * (nz - d * 2),0);
	psi2X.assign(2 * nb * (nz - d * 2),0);
	psi2XLa.assign(2 * nb * (nz - d * 2),0);
	phi2X.assign(2 * nb * (nz - d * 2),0);
	phi2XLa.assign(2 * nb * (nz - d * 2),0);
	u020BXLa.assign(2 * nb * (nz - d * 2),0);

	psiZ.resize(2 * nb * (nx - d * 2));
	psiZLa.resize(2 * nb * (nx - d * 2));
	phiZ.resize(2 * nb * (nx - d * 2));
	phiZLa.resize(2 * nb * (nx - d * 2));
	EtaZ.resize(2 * nb * (nx - d * 2));
	EtaZLa.resize(2 * nb * (nx - d * 2));
	psi2Z.resize(2 * nb * (nx - d * 2));
	psi2ZLa.resize(2 * nb * (nx - d * 2));
	phi2Z.resize(2 * nb * (nx - d * 2));
	phi2ZLa.resize(2 * nb * (nx - d * 2));
	u002BZLa.resize(2 * nb * (nx - d * 2));

	psiZ.assign(2 * nb * (nx - d * 2),0);
	psiZLa.assign(2 * nb * (nx - d * 2),0);
	phiZ.assign(2 * nb * (nx - d * 2),0);
	phiZLa.assign(2 * nb * (nx - d * 2),0);
	EtaZ.assign(2 * nb * (nx - d * 2),0);
	EtaZLa.assign(2 * nb * (nx - d * 2),0);
	psi2Z.assign(2 * nb * (nx - d * 2),0);
	psi2ZLa.assign(2 * nb * (nx - d * 2),0);
	phi2Z.assign(2 * nb * (nx - d * 2),0);
	phi2ZLa.assign(2 * nb * (nx - d * 2),0);
	u002BZLa.assign(2 * nb * (nx - d * 2),0);

	ux.resize(nx * nz);
	uxLa.resize(nx * nz);
	uz.resize(nx * nz);
	uzLa.resize(nx * nz);

	ux.assign(nx * nz,0);
	uxLa.assign(nx * nz,0);
	uz.assign(nx * nz,0);
	uzLa.assign(nx * nz,0);

	psixlen = nz - d * 2;
	psizlen = 2 * nb;

	count = 0;
	init = true;
}

void CPML::applyCPML(float *uLa, float *u, float *uNe, const float *vel, const int nx, const int nz, const ForwardModeling &fm) {
	if(count == 0) {
		initCPML(nx, nz, fm);
		INFO() << "Warning: the same CPML variables can not be used in two different forward modeling at the same time!!!";
		if(!init) {
			INFO() << "CPML no initialization!!";
			exit(1);
		}
		else {
			init = false;
		}
	}
	count = (count + 1) % fm.getnt();
	int d = fm.getFDLEN();
	int bx0 = fm.getbx0();
	int bxn = fm.getbxn();
	int bz0 = fm.getbz0();
	int bzn = fm.getbzn();
	float dx = fm.getdx();
	float dt = fm.getdt();
	int nBMPosX, nBMPosZ;
	float lB, dDlB, dD2lB, alphaDlB, alphaD2lB;
	float DlB0; //Larger ValuedxB leads to stronger attenation
	float alphaDlB0; //Larger ValueAlphaxB leads to faster phase variation
	int cnx = nx - bx0 - bxn;
	int cnz = nz - bz0 - bzn;
	float blengthx = cnx * dx;
	float blengthz = cnz * dx;
	float u020 = 0.0f, u002 = 0.0f, aB, bB;
#ifdef USE_OPENMP
	#pragma omp parallel for private(nBMPosX, nBMPosZ, lB, dDlB, dD2lB, alphaDlB, alphaD2lB, DlB0, alphaDlB0, u020, u002, aB, bB)
#endif
	for(int ix = d ; ix < nx - d ; ix ++) { 
		for(int iz = d ; iz < nz - d ; iz ++) {
			ux[ix * nz + iz] = (2. / 3. * (u[(ix + 1) * nz + iz] - u[(ix - 1) * nz + iz]) - 1. / 12. * (u[(ix + 2) * nz + iz] - u[(ix - 2) * nz + iz])) / dx;
			uz[ix * nz + iz] = (2. / 3. * (u[ix * nz + iz + 1] - u[ix * nz + iz - 1]) - 1. / 12. * (u[ix * nz + iz + 2] - u[ix * nz + iz - 2])) / dx;

			if (ix < bx0 || ix >= nx - bxn || iz < bz0 || iz >= nz - bzn)
			{
				u020 = 1.0 / 12.0 * (-30 * u[ix * nz + iz] + 16 *(u[(ix - 1) * nz + iz]+u[(ix + 1) * nz + iz]) - (u[(ix - 2) * nz + iz]+u[(ix + 2) * nz + iz])) / dx / dx;	
				u002 = 1.0 / 12.0 * (-30 * u[ix * nz + iz] + 16 *(u[ix * nz + iz - 1]+u[ix * nz + iz + 1]) - (u[ix * nz + iz - 2]+u[ix * nz + iz + 2])) / dx / dx;

				//X boundaries
				if(ix < bx0 || ix >= nx - bxn)
				{
					DlB0 = 80000; //Larger ValuedxB leads to stronger attenation
					alphaDlB0 = 0; //Larger ValueAlphaxB leads to faster phase variation
					//Generate boundary coordinates
					GetXBoundaryMPos(ix, iz, &nBMPosX, &nBMPosZ, nx, nz, fm);
					//Compute the convolutions
					if(ix < bx0)
						lB = (ix - bx0) * dx;
					else
						lB = (ix - (nx - bxn - 1)) * dx;
					dDlB = DlB0 * (lB / blengthx) * (lB / blengthx);
					dD2lB = DlB0 * 2 * lB / (blengthx * blengthx);
					//alphaDlB = alphaDlB0 * fabs(lB) / blengthx;
					//alphaD2lB = SIGN(alphaDlB0 / blengthx, lB);
					alphaDlB = alphaDlB0 * (1 - fabs(lB) / blengthx);
					alphaD2lB = -SIGN(alphaDlB0 / blengthx, lB);
					
					bB = exp(-(dDlB + alphaDlB) * dt);
					aB = (1 - bB) / (dDlB + alphaDlB);
					psi2X[nBMPosX * psixlen + nBMPosZ] = bB * psi2XLa[nBMPosX * psixlen + nBMPosZ] \
						+ aB * 0.5 * (u020 + u020BXLa[nBMPosX * psixlen + nBMPosZ]);
					phi2X[nBMPosX * psixlen + nBMPosZ] = bB * phi2XLa[nBMPosX * psixlen + nBMPosZ] \
						+ aB * 0.5 * (psi2X[nBMPosX * psixlen + nBMPosZ] + psi2XLa[nBMPosX * psixlen + nBMPosZ]);
					psiX[nBMPosX * psixlen + nBMPosZ] = bB * psiXLa[nBMPosX * psixlen + nBMPosZ] \
						+ aB * 0.5 * (ux[ix * nz + iz] + uxLa[ix * nz + iz]);
					phiX[nBMPosX * psixlen + nBMPosZ] = bB * phiXLa[nBMPosX * psixlen + nBMPosZ] \
						+ aB * 0.5 * (psiX[nBMPosX * psixlen + nBMPosZ] + psiXLa[nBMPosX * psixlen + nBMPosZ]);
					EtaX[nBMPosX * psixlen + nBMPosZ] = bB * EtaXLa[nBMPosX * psixlen + nBMPosZ] \
						+ aB * 0.5 * (phiX[nBMPosX * psixlen + nBMPosZ] + phiXLa[nBMPosX * psixlen + nBMPosZ]);

					//Update the former status
					u020BXLa[nBMPosX * psixlen + nBMPosZ] = u020;
					psi2XLa[nBMPosX * psixlen + nBMPosZ] = psi2X[nBMPosX * psixlen + nBMPosZ];
					phi2XLa[nBMPosX * psixlen + nBMPosZ] = phi2X[nBMPosX * psixlen + nBMPosZ];
					psiXLa[nBMPosX * psixlen + nBMPosZ] = psiX[nBMPosX * psixlen + nBMPosZ];
					phiXLa[nBMPosX * psixlen + nBMPosZ] = phiX[nBMPosX * psixlen + nBMPosZ];
					EtaXLa[nBMPosX * psixlen + nBMPosZ] = EtaX[nBMPosX * psixlen + nBMPosZ];

					u020 = u020\
								 - 2 * dDlB * psi2X[nBMPosX * psixlen + nBMPosZ] + (dDlB * dDlB) * phi2X[nBMPosX * psixlen + nBMPosZ]\
								 - dD2lB * psiX[nBMPosX * psixlen + nBMPosZ]\
								 + dDlB * (2 * dD2lB + alphaD2lB) * phiX[nBMPosX * psixlen + nBMPosZ]\
								 - (dDlB * dDlB) * (dD2lB+alphaD2lB) * EtaX[nBMPosX * psixlen + nBMPosZ];
				}

				//Z boundaries
				if(iz < bz0 || iz >= nz - bzn)
				{
					if(iz < bz0) {
						DlB0 = 2000; //Larger ValuedxB leads to stronger attenation
						alphaDlB0 = 0; //Larger ValueAlphaxB leads to faster phase variation
					}
					else {
						DlB0 = 4000;
						alphaDlB0 = 0; //Larger ValueAlphaxB leads to faster phase variation
					}
					//Generate boundary coordinates
					GetZBoundaryMPos(ix, iz, &nBMPosX, &nBMPosZ, nx , nz, fm);
					//Compute the convolutions
					if(iz < bz0)
						lB = (iz - bz0) * dx;
					else
						lB = (iz - (nz - bzn - 1)) * dx;
					dDlB = DlB0 * (lB / blengthz) * (lB / blengthz);
					dD2lB = DlB0 * 2 * lB / (blengthz * blengthz);
					alphaDlB = alphaDlB0 * fabs(lB) / blengthz;
					alphaD2lB = SIGN(alphaDlB0 / blengthz, lB);
					
					bB = exp(-(dDlB + alphaDlB) * dt);
					aB = (1 - bB) / (dDlB + alphaDlB);
					psi2Z[nBMPosX * psizlen + nBMPosZ] = bB * psi2ZLa[nBMPosX * psizlen + nBMPosZ] \
						+ aB * 0.5 * (u002 + u002BZLa[nBMPosX * psizlen + nBMPosZ]);
					phi2Z[nBMPosX * psizlen + nBMPosZ] = bB * phi2ZLa[nBMPosX * psizlen + nBMPosZ] \
						+ aB * 0.5 * (psi2Z[nBMPosX * psizlen + nBMPosZ] + psi2ZLa[nBMPosX * psizlen + nBMPosZ]);
					psiZ[nBMPosX * psizlen + nBMPosZ] = bB * psiZLa[nBMPosX * psizlen + nBMPosZ] \
						+ aB * 0.5 * (uz[ix * nz + iz] + uzLa[ix * nz + iz]);
					phiZ[nBMPosX * psizlen + nBMPosZ] = bB * phiZLa[nBMPosX * psizlen + nBMPosZ] \
						+ aB * 0.5 * (psiZ[nBMPosX * psizlen + nBMPosZ] + psiZLa[nBMPosX * psizlen + nBMPosZ]);
					EtaZ[nBMPosX * psizlen + nBMPosZ] = bB * EtaZLa[nBMPosX * psizlen + nBMPosZ] \
						+ aB * 0.5 * (phiZ[nBMPosX * psizlen + nBMPosZ] + phiZLa[nBMPosX * psizlen + nBMPosZ]);

					//Update the former status
					u002BZLa[nBMPosX * psizlen + nBMPosZ] = u002;
					psi2ZLa[nBMPosX * psizlen + nBMPosZ] = psi2Z[nBMPosX * psizlen + nBMPosZ];
					phi2ZLa[nBMPosX * psizlen + nBMPosZ] = phi2Z[nBMPosX * psizlen + nBMPosZ];
					psiZLa[nBMPosX * psizlen + nBMPosZ] = psiZ[nBMPosX * psizlen + nBMPosZ];
					phiZLa[nBMPosX * psizlen + nBMPosZ] = phiZ[nBMPosX * psizlen + nBMPosZ];
					EtaZLa[nBMPosX * psizlen + nBMPosZ] = EtaZ[nBMPosX * psizlen + nBMPosZ];
					//Update u002
					u002 = u002 \
						- 2 * dDlB * psi2Z[nBMPosX * psizlen + nBMPosZ] + (dDlB * dDlB) * phi2Z[nBMPosX * psizlen + nBMPosZ] \
						- dD2lB * psiZ[nBMPosX * psizlen + nBMPosZ] \
						+ dDlB * (2 * dD2lB + alphaD2lB) * phiZ[nBMPosX * psizlen + nBMPosZ] \
						- (dDlB * dDlB) * (dD2lB + alphaD2lB) * EtaZ[nBMPosX * psizlen + nBMPosZ];
				}
				uNe[ix * nz + iz] = 2 * u[ix * nz + iz] - uLa[ix * nz + iz] + (1.0f / vel[ix * nz + iz]) * (u002 + u020) * dx * dx;
			}
		}
	}
	std::swap(ux, uxLa);
	std::swap(uz, uzLa);
}

