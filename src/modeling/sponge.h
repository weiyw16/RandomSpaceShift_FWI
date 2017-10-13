#ifndef SPONGE_H
#define SPONGE_H
#include <vector>
#include <cmath>
class Sponge {
	public:
		void initbndr(int nb);
		void applySponge(float* p, const float *vel, int nx, int nz, int nb, float dt, float dx, int freeSurface);
	private:
		std::vector<float> bndr;
};
#endif
