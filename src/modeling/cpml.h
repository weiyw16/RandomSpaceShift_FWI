#ifndef CPML_H
#define CPML_H
#include <vector>
class ForwardModeling;
class CPML {
	public:
		CPML();
		void initCPML(const int nx, const int nz, const ForwardModeling &fm);
		void applyCPML(float *uLa, float *u, float *uNe, const float *vel, const int nx, const int nz, const ForwardModeling &fm);

	private:
		void GetXBoundaryMPos(int xPos, int zPos, int *xMPos, int *zMPos, const int nx, const int nz, const ForwardModeling &fm) const;
		void GetZBoundaryMPos(int xPos, int zPos, int *xMPos, int *zMPos, const int nx, const int nz, const ForwardModeling &fm) const;

	private:
		std::vector<float> psiX, psiXLa, phiX, phiXLa, EtaX, EtaXLa, psi2X, psi2XLa, phi2X, phi2XLa, u020BXLa;
		std::vector<float> psiZ, psiZLa, phiZ, phiZLa, EtaZ, EtaZLa, psi2Z, psi2ZLa, phi2Z, phi2ZLa, u002BZLa;
		std::vector<float> ux, uz, uxLa, uzLa;
		int psixlen, psizlen;
		bool init;
		int count;
};
#endif
