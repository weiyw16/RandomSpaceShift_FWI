/*
 * forwardmodeling.cpp
 *
 *  Created on: Feb 29, 2016
 *      Author: rice
 */

#include <cmath>
#include <functional>
#include "forwardmodeling.h"
#include "logger.h"
#include "sum.h"
#include "sfutil.h"
#include "common.h"

extern "C" {
#include <rsf.h>
#include "fdutil.h"
#include "fd4t10s-damp-zjh.h"
#include "fd4t10s-zjh-born.h"
#include "fd4t10s-zjh.h"
#include "fd4t10s-nobndry.h"
}
#include <sys/time.h>

#ifdef USE_SW
extern "C" {
#include "fd4t10s-damp-zjh-cg.h"
}
#endif

CPML* ForwardModeling::getCPML(int cpmlId) const{
	return cpml[cpmlId];
}

void ForwardModeling::initFdUtil(sf_file &vinit, Velocity *v, int nb, float dx, float dt) {
	const int ompchunk = 8;
	fd = fdutil_init(false, false, sf_iaxa(vinit, 1), sf_iaxa(vinit, 2), bx0, ompchunk); 
	sp = sponge_make(bx0);
	float *v2 = (float *)malloc(sizeof(float) * v->nx * v->nz);
	float **vv;
	for(int i = 0 ; i < v->nx * v->nz ; i ++)
		v2[i] = 1.0 / v->dat[i] * dx * dx / dt / dt;
	vv = f1dto2d(v2, v->nx, v->nz);
	abc = abcone2d_make(EXFDBNDRYLEN, dt, vv, false, fd);
}

static void fillForStencil(Velocity &exvel, int halo) {
  int nxpad = exvel.nx;
  int nzpad = exvel.nz;
  int nx = nxpad - 2 * halo;
  int nz = nzpad - 2 * halo;

  std::vector<float> &vel_e = exvel.dat;

  //expand z direction first
  for (int ix = halo; ix < nx + halo; ix++) {
    for (int iz = 0; iz < halo; iz++) {
      vel_e[ix * nzpad + iz] = vel_e[ix * nzpad + halo];               // top
    }
    for (int iz = nz + halo; iz < nzpad; iz ++) {
      vel_e[ix * nzpad + iz] = vel_e[ix * nzpad + (nz + halo - 1)]; // bottom
    }
  }

  //Then x direction
  for (int iz = 0; iz < nzpad; iz++) {
    for (int ix = 0; ix < halo; ix++) {
      vel_e[ix * nzpad + iz] = vel_e[halo * nzpad + iz];               // left
    }
    for (int ix = nx + halo; ix < nxpad; ix++) {
      vel_e[ix * nzpad + iz] = vel_e[(nx + halo - 1) * nzpad + iz]; // right
    }
  }
}

static void expandForStencil(Velocity &exvel, const Velocity &v0, int halo) {
  int nx = v0.nx;
  int nz = v0.nz;
  int nzpad = nz + 2 * halo;

  const std::vector<float> &vel = v0.dat;
  std::vector<float> &vel_e = exvel.dat;

  //copy the vel into vel_e
  for (int ix = halo; ix < nx + halo; ix++) {
    for (int iz = halo; iz < nz + halo; iz++) {
      vel_e[ix * nzpad + iz] = vel[(ix - halo) * nz +  (iz - halo)];
    }
  }

  fillForStencil(exvel, halo);

}

static void expandBndry(Velocity &exvel, const Velocity &v0, int nb, int freeSurface) {
  int nx = v0.nx;
  int nz = v0.nz;
  int nxpad = nx + 2 * nb;
	int nzpad;
	if(freeSurface)
		nzpad = nz + nb;
	else
		nzpad = nz + 2 * nb;
  const std::vector<float> &a = v0.dat;
  std::vector<float> &b = exvel.dat;

	if(freeSurface) {
		/// internal
		for (int ix = 0; ix < nx; ix++) {
			for (int iz = 0; iz < nz; iz++) {
				b[(nb + ix) * nzpad + iz] = a[ix * nz + iz];
			}
		}

		/// boundary, free surface
		for (int ix = 0; ix < nb; ix++) {
			for (int iz = 0; iz < nz; iz++) {
				b[ix * nzpad + iz] = a[iz];                              /* left */
				b[(nb + nx + ix) * nzpad + iz] = a[(nx - 1) * nz + iz];  /* right */
			}
		}

		for (int ix = 0; ix < nxpad; ix++) {
			for (int iz = 0; iz < nb; iz++) {
				b[ix * nzpad + (nz + iz)] = b[ix * nzpad + (nz - 1)];  /* bottom*/
			}
		}
	}
	else {
		/// internal
		for (int ix = 0; ix < nx; ix++) {
			for (int iz = 0; iz < nz; iz++) {
				b[(nb + ix) * nzpad + (nb + iz)] = a[ix * nz + iz];
			}
		}

		/// boundary
		for (int ix = 0; ix < nb; ix++) {
			for (int iz = 0; iz < nz; iz++) {
				b[ix * nzpad + nb + iz] = a[iz];                              /* left */
				b[(nb + nx + ix) * nzpad + nb + iz] = a[(nx - 1) * nz + iz];  /* right */
			}
		}

		for (int ix = 0; ix < nxpad; ix++) {
			for (int iz = 0; iz < nb; iz++) {
				b[ix * nzpad + iz] = b[ix * nzpad + nb];												 /* top*/
				b[ix * nzpad + (nb + nz + iz)] = b[ix * nzpad + (nb + nz - 1)];  /* bottom*/
			}
		}
	}
}

static void transvel(std::vector<float> &vel, float dx, float dt) {
  for (size_t i = 0; i < vel.size(); i ++) {
    vel[i] = (dx * dx) / (dt * dt * vel[i] * vel[i]);
  }
}

static void recoverVel(std::vector<float> &vel, float dx, float dt) {
  for (size_t i = 0; i < vel.size(); i ++) {
    vel[i] = std::sqrt(dx*dx / (dt*dt*vel[i]));
  }
}

Velocity ForwardModeling::expandDomain(const Velocity& _vel) {
  // expand for boundary, free surface
  int nb = bx0 - EXFDBNDRYLEN;

	/*
	sf_file sf_v1 = sf_output("v0_before.rsf");
	sf_putint(sf_v1, "n1", _vel.nz);
	sf_putint(sf_v1, "n2", _vel.nx);
	sf_floatwrite(const_cast<float*>(&_vel.dat[0]), _vel.nz * _vel.nx, sf_v1);
	*/

  Velocity exvelForBndry;
	if(freeSurface) {
		exvelForBndry.resize(_vel.nx + 2 * nb, _vel.nz + nb);
		expandBndry(exvelForBndry, _vel, nb, freeSurface);
	}
	else {
		exvelForBndry.resize(_vel.nx + 2 * nb, _vel.nz + 2 * nb);
		expandBndry(exvelForBndry, _vel, nb, freeSurface);
	}

	/*
	sf_file sf_v2 = sf_output("v0_after.rsf");
	sf_putint(sf_v2, "n1", _vel.nz + 2 * nb);
	sf_putint(sf_v2, "n2", _vel.nx + 2 * nb);
	sf_floatwrite(const_cast<float*>(&exvelForBndry.dat[0]), (_vel.nz + 2 * nb) * (_vel.nx + 2 * nb), sf_v2);
	exit(1);
	*/

	/*
	sf_file sf_v1 = sf_output("v0_bndry.rsf");
	sf_putint(sf_v1, "n1", _vel.nz + nb);
	sf_putint(sf_v1, "n2", _vel.nx + 2 * nb);
	sf_floatwrite(const_cast<float*>(&exvelForBndry.dat[0]), (_vel.nz + nb) * (_vel.nx + 2 * nb), sf_v1);
	exit(1);
	*/

  transvel(exvelForBndry.dat, dx, dt);

	/*
	sf_file sf_v1 = sf_output("v0_before.rsf");
	sf_putint(sf_v1, "n1", _vel.nz + nb);
	sf_putint(sf_v1, "n2", _vel.nx + 2 * nb);
	sf_floatwrite(const_cast<float*>(&exvelForBndry.dat[0]), (_vel.nz + nb) * (_vel.nx + 2 * nb), sf_v1);
	exit(1);
	*/

  // expand for stencil
  Velocity ret(exvelForBndry.nx+2*EXFDBNDRYLEN, exvelForBndry.nz+2*EXFDBNDRYLEN);
  expandForStencil(ret, exvelForBndry, EXFDBNDRYLEN);

	/*
	sf_file sf_v1 = sf_output("v0_after.rsf");
	sf_putint(sf_v1, "n1", ret.nz);
	sf_putint(sf_v1, "n2", ret.nx);
	sf_floatwrite(const_cast<float*>(&ret.dat[0]), ret.nx * ret.nz, sf_v1);
	exit(1);
	*/

  return ret;
}

void ForwardModeling::addBornwv(float *fullwv_t0, float *fullwv_t1, float *fullwv_t2, const float *exvel_m, float dt, int it, float *rp1) const {
	int nx = vel->nx;
	int nz = vel->nz;

	if(it == 0) {
		for(int i = bx0 ; i < nx - bxn ; i ++)
			for(int j = bz0 ; j < nz - bzn ; j ++)
				rp1[i * nz + j] += 2 * (fullwv_t2[i * nz + j] - fullwv_t1[i * nz + j]) / vel->dat[i * nz + j] * exvel_m[i * nz + j] / dt;
	}
	else if(it == nt - 1) {
		for(int i = bx0 ; i < nx - bxn ; i ++)
			for(int j = bz0 ; j < nz - bzn ; j ++)
				rp1[i * nz + j] += 2 * (fullwv_t1[i * nz + j] - fullwv_t0[i * nz + j]) / vel->dat[i * nz + j] * exvel_m[i * nz + j] / dt;
	}
	else {
		for(int i = bx0 ; i < nx - bxn ; i ++)
			for(int j = bz0 ; j < nz - bzn ; j ++)
				rp1[i * nz + j] -= 2 * (fullwv_t2[i * nz + j] - 2 * fullwv_t1[i * nz + j] + fullwv_t0[i * nz + j]) / vel->dat[i * nz + j] * exvel_m[i * nz + j];
	}
}

void ForwardModeling::swStepForward(std::vector<float> &p0, std::vector<float> &p1) const {
  static std::vector<float> u2(vel->nx * vel->nz, 0);

	
  struct timeval t1, t2;	
	float gflop = 0;
	//damp
  gettimeofday(&t1, NULL);
  fd4t10s_damp_zjh_2d_vtrans_test(&p0[0], &p1[0], &vel->dat[0], &u2[0], vel->nx, vel->nz, bx0, freeSurface);
  gettimeofday(&t2, NULL);
	printf("Finish stepForward on CPU. Time: %0.9lfs FLOPS:%.5f GFLOPS\n\n\n", TIME(t1, t2), gflop/(TIME(t1, t2))); 
	//sponge
  //fd4t10s_nobndry_2d_vtrans(&p0[0], &p1[0], &vel->dat[0], &u2[0], vel->nx, vel->nz, bx0, freeSurface);
	//spng->applySponge(&p0[0], &vel->dat[0], vel->nx, vel->nz, bx0, dt, dx, freeSurface);
	//spng->applySponge(&p1[0], &vel->dat[0], vel->nx, vel->nz, bx0, dt, dx, freeSurface);

	//fdUtil
	//float **pp0 = f1dto2d(p0, vel->nx, vel->nz);
	//float **pp1 = f1dto2d(p1, vel->nx, vel->nz);
	//abcone2d_apply(pp0, pp1, EXFDBNDRYLEN, abc, fd);
	//sponge2d_apply(pp0, sp, fd);
	//sponge2d_apply(pp1, sp, fd);
	

}

void ForwardModeling::stepForward(std::vector<float> &p0, std::vector<float> &p1) const {

  static std::vector<float> u2(vel->nx * vel->nz, 0);
//#ifdef USE_SW
  //struct timeval t1, t2;	
	//float gflop = 0;
	//damp
  //gettimeofday(&t1, NULL);
  //static std::vector<float> p2(vel->nx * vel->nz, 0);
  //fd4t10s_nobndry_zjh_2d_vtrans_cg(&p0[0], &p1[0], &p2[0], &vel->dat[0], &u2[0], vel->nx, vel->nz, bx0, nt, freeSurface);
//	std::swap(p0, p2);
//#else
	//spng->applySponge(&p0[0], &vel->dat[0], vel->nx, vel->nz, bx0, dt, dx, freeSurface);
	//spng->applySponge(&p1[0], &vel->dat[0], vel->nx, vel->nz, bx0, dt, dx, freeSurface);
  //gettimeofday(&t2, NULL);
	//printf("Finish stepForward on accelerator. Time: %0.9lfs FLOPS:%.5f GFLOPS\n\n\n", TIME(t1, t2), gflop/(TIME(t1, t2))); 
	
	//sponge
  fd4t10s_nobndry_2d_vtrans(&p0[0], &p1[0], &vel->dat[0], &u2[0], vel->nx, vel->nz, bx0, freeSurface);
	spng->applySponge(&p0[0], &vel->dat[0], vel->nx, vel->nz, bx0, dt, dx, freeSurface);
	spng->applySponge(&p1[0], &vel->dat[0], vel->nx, vel->nz, bx0, dt, dx, freeSurface);
//#endif 

	//fdUtil
	//float **pp0 = f1dto2d(p0, vel->nx, vel->nz);
	//float **pp1 = f1dto2d(p1, vel->nx, vel->nz);
	//abcone2d_apply(pp0, pp1, EXFDBNDRYLEN, abc, fd);
	//sponge2d_apply(pp0, sp, fd);
	//sponge2d_apply(pp1, sp, fd);
}
void ForwardModeling::stepForward(std::vector<float> &p0, std::vector<float> &p1, bool vtrans) const {
	static std::vector<float> u2(vel->nx * vel->nz, 0);
	if(vtrans){
	fd4t10s_nobndry_2d_vtrans(&p0[0], &p1[0], &vel->dat[0], &u2[0], vel->nx, vel->nz, bx0, freeSurface);
	spng->applySponge(&p0[0], &vel->dat[0], vel->nx, vel->nz, bx0, dt, dx, freeSurface);
	spng->applySponge(&p1[0], &vel->dat[0], vel->nx, vel->nz, bx0, dt, dx, freeSurface);
	}
	else{
	int ix,iz;
	int nx=vel->nx;
	int nz=vel->nz;
	//float *vel_trans[nx*nz] = 0.0;
	std::vector<float> vel_trans(nx * nz, 0);
	//float *vel0 = &vel->dat[0];
	for(ix = 0; ix < nx; ix++){
	for(iz = 0; iz < nz; iz++){
		vel_trans[ix*nz + iz] = dx * dx / (dt * dt * vel->dat[ix*nz + iz] * vel->dat[ix*nz + iz]);
	}
	}
	fd4t10s_nobndry_2d_vtrans(&p0[0], &p1[0], &vel_trans[0], &u2[0], vel->nx, vel->nz, bx0, freeSurface);
	spng->applySponge(&p0[0], &vel->dat[0], vel->nx, vel->nz, bx0, dt, dx, freeSurface);
	spng->applySponge(&p1[0], &vel->dat[0], vel->nx, vel->nz, bx0, dt, dx, freeSurface);
	}	
}

	std::vector<float> ForwardModeling::getBornCoff(const Velocity &localvel, const Velocity &localvel_real , float dx, float dt) 
{
	int i;
	//float *velm=NULL;
	//float *vel0 = &vel->dat[0];
	//float *bco[nx*nz] = 0.0;
 	//std::vector<float> vel_m(vel->nx * vel->nz, 0.0f);
	std::vector<float> vel_trans(vel->nx * vel->nz, 0.0f);
	std::vector<float> vel_real_trans(vel->nx * vel->nz, 0.0f);
	//std::vector<float> vel_0(nx * nz, vel->dat);
	std::vector<float> bff(vel->nx * vel->nz, 0.0f);
	//vectorMinus(vel_real.dat, vel.dat, vel_m);
	vel_trans = localvel.dat;
	vel_real_trans = localvel_real.dat;
	//recoverVel(vel_trans, dx, dt);
	//recoverVel(vel_real_trans, dx, dt);
	//vectorMinus(vel_real_trans, vel_trans, vel_m);

	//velm = &vel_m[0];

	for ( i = 0; i < vel->nx * vel->nz; i++ ){
		bff[ i ] = 2.0 * (  vel_real_trans[ i] - vel_trans[i] ) / ( vel_trans[i]  *dx * dx);
	}
	//bcoff.assign(&bco[0],&bco[0]+nx*nz);
	return bff;
}

void ForwardModeling::stepbornForward(std::vector<float> &p0, std::vector<float> &p1) const {
	fd4t10s_zjh_born(&p0[0], &p1[0], &bcoff[0], vel->nx, vel->nz);
	spng->applySponge(&p0[0], &vel->dat[0], vel->nx, vel->nz, bx0, dt, dx, freeSurface);
	spng->applySponge(&p1[0], &vel->dat[0], vel->nx, vel->nz, bx0, dt, dx, freeSurface);
}

ForwardModeling::~ForwardModeling() {
#ifdef USE_SW
	fd4t10s_4cg_exit();
#endif
}

void ForwardModeling::stepForward(std::vector<float> &p0, std::vector<float> &p1, int cpmlId) const {
  static std::vector<float> u2(vel->nx * vel->nz, 0);
  static std::vector<float> p2(vel->nx * vel->nz, 0);

  fd4t10s_nobndry_2d_vtrans_3vars(&p0[0], &p1[0], &p2[0], &vel->dat[0], &u2[0], vel->nx, vel->nz, bx0, freeSurface);
	cpml[cpmlId]->applyCPML(&p0[0], &p1[0], &p2[0], &vel->dat[0], vel->nx, vel->nz, *this);
	std::swap(p0, p2);
}

void ForwardModeling::bindVelocity(const Velocity& _vel) {
  this->vel = &_vel;
}

void ForwardModeling::bindRealVelocity(const Velocity& _vel) {
  this->vel_real = &_vel;
}

void ForwardModeling::bindBornCoff( std::vector<float> &b) {
  this->bcoff = b;
	sf_file fbcoff = sf_output("bcoff.rsf");
	sf_floatwrite(&bcoff[0], b.size(), fbcoff);
}

void ForwardModeling::recordSeis(float* seis_it, const float* p,
    const ShotPosition& geoPos) const {

  int ng = geoPos.ns;
  int nzpad = vel->nz;

//  DEBUG() << format("ng %d") % ng;
//  float sum = 0;
  for (int ig = 0; ig < ng; ig++) {
    int gx = geoPos.getx(ig) + bx0;
    int gz = geoPos.getz(ig) + bz0;	
    int idx = gx * nzpad + gz;
    seis_it[ig] = p[idx];
//    DEBUG() << format("ig %d, idx %d, v %.20f") % ig % idx % seis_it[ig];
  }
//  DEBUG() << format("sum %.20f") % sum;

}



const Velocity& ForwardModeling::getVelocity() const {
  return *vel;
}

/*
void ForwardModeling::stepBackward(std::vector<float> &p0, std::vector<float> &p1) const {
  static std::vector<float> u2(vel->nx * vel->nz, 0);
  fd4t10s_zjh_2d_vtrans(&p0[0], &p1[0], &vel->dat[0], &u2[0], vel->nx, vel->nz);
}
*/

void ForwardModeling::stepBackward(std::vector<float> &p0, std::vector<float> &p1) const {
  static std::vector<float> u2(vel->nx * vel->nz, 0);
#ifdef USE_SW
  static std::vector<float> p2(vel->nx * vel->nz, 0);
  fd4t10s_nobndry_zjh_2d_vtrans_cg(&p0[0], &p1[0], &p2[0], &vel->dat[0], &u2[0], vel->nx, vel->nz, bx0, nt, freeSurface);
	std::swap(p0, p2);
#else
  fd4t10s_zjh_2d_vtrans(&p0[0], &p1[0], &vel->dat[0], &u2[0], vel->nx, vel->nz);
#endif
}


void ForwardModeling::addSource(float* p, const float* source,
    const ShotPosition& pos) const
{
  manipSource(p, source, pos, std::plus<float>());
}

void ForwardModeling::subSource(float* p, const float* source,
    const ShotPosition& pos) const {
  manipSource(p, source, pos, std::minus<float>());
}

void ForwardModeling::manipSource(float* p, const float* source,
    const ShotPosition& pos, boost::function2<float, float, float> op) const {
  int nzpad = vel->nz;

  for (int is = 0; is < pos.ns; is++) {
    int sx = pos.getx(is) + bx0;
    int sz = pos.getz(is) + bz0; 
    p[sx * nzpad + sz] = op(p[sx * nzpad + sz], source[is]);
//    DEBUG() << format("sx %d, sz %d, source[%d] %f") % sx % sz % is % source[is];
  }
}

void ForwardModeling::bornMaskGradient(float* grad, int H) const {
  int nxpad = vel->nx;
  int nzpad = vel->nz;
#ifdef USE_OPENMP
	#pragma omp parallel for
#endif
	for (int h = -H ; h <= H ; h ++) {
		int ind = h + H;
		for (int ix = 0; ix < nxpad; ix++) {
			for (int iz = 0; iz < nzpad; iz++) {
				if (ix < bx0 || iz < bz0 || ix >= nxpad - bxn || iz >= nzpad - bzn) {
					grad[ind * nxpad * nzpad + ix  * nzpad + iz] = 0.f;
				}
			}
		}
	}
}

void ForwardModeling::maskGradient(float* grad) const {
  int nxpad = vel->nx;
  int nzpad = vel->nz;
  for (int ix = 0; ix < nxpad; ix++) {
    for (int iz = 0; iz < nzpad; iz++) {
      if (ix < bx0 || iz < bz0 || ix >= nxpad - bxn || iz >= nzpad - bzn) {
        grad[ix  * nzpad + iz] = 0.f;
      }
    }
  }
}

void ForwardModeling::refillBoundary(float* gradient) const {
  int nzpad = vel->nz;
  int nxpad = vel->nx;

  for (int ix = 0; ix < nxpad; ix++) {
    for (int iz = 0; iz < bz0; iz++) {
      gradient[ix * nzpad + iz] = gradient[ix * nzpad + bz0];           // top
    }
    for (int iz = nzpad - bzn; iz < nzpad; iz++) {
      gradient[ix * nzpad + iz] = gradient[ix * nzpad + nzpad - bzn - 1];  // bottom
    }
  }

  for (int ix = 0; ix < bx0; ix++) {                              // left
    for (int iz = 0; iz < nzpad; iz++) {
      gradient[ix * nzpad + iz] = gradient[bx0 * nzpad + iz];
    }
  }

  for (int ix = nxpad - bxn; ix < nxpad; ix++) {                        // right
    for (int iz = 0; iz < nzpad; iz++) {
      gradient[ix * nzpad + iz] = gradient[(nxpad - bxn - 1) * nzpad + iz];
    }
  }
}

void ForwardModeling::sfWriteVel(const std::vector<float> &exvel, sf_file file) const {
  //assert(exvel.size() == vel->dat.size());
  int nzpad = vel->nz;
  int nxpad = vel->nx;
  int nz = nzpad - bz0 - bzn;

  std::vector<float> vv = exvel;
  recoverVel(vv, dx, dt);
  for (int ix = bx0; ix < nxpad - bxn; ix++) {
    sf_floatwrite(&vv[ix * nzpad + EXFDBNDRYLEN], nz, file);
  }
}

void ForwardModeling::refillVelStencilBndry() {
  Velocity &exvel = getVelocity();
  fillForStencil(exvel, EXFDBNDRYLEN);
}

void ForwardModeling::FwiForwardModeling(const std::vector<float>& encSrc,
    std::vector<float>& dcal, int shot_id) const {
  int nx = getnx();
  int nz = getnz();
  int ns = getns();
  int ng = getng();

  std::vector<float> p0(nz * nx, 0);
  std::vector<float> p1(nz * nx, 0);
  ShotPosition curSrcPos = allSrcPos->clipRange(shot_id, shot_id);

  /*
	sf_file sf_v0 = sf_input("/home/cbw/fwijob/fwi_test2/v0.rsf");
	sf_floatread(const_cast<float*>(&vel->dat[0]), nz * nx, sf_v0);
  */

  /*
	sf_file sf_v1 = sf_output("v1.rsf");
	sf_putint(sf_v1, "n1", nz);
	sf_putint(sf_v1, "n2", nx);
	sf_floatwrite(const_cast<float*>(&vel->dat[0]), nz * nx, sf_v1);
  exit(1);
  */

  /*
	sf_file sf_p0 = sf_input("/home/rice/cbw/pfwi/job/p0.rsf");
	sf_floatread(const_cast<float*>(&p0[0]), nz * nx, sf_p0);
  
	sf_file sf_p1 = sf_input("/home/rice/cbw/pfwi/job/p1.rsf");
	sf_floatread(const_cast<float*>(&p1[0]), nz * nx, sf_p1);
  */

  for(int it=0; it<nt; it++) {
    addSource(&p1[0], &encSrc[it], curSrcPos);

    /*
    sf_file sf_p0 = sf_output("pp0.rsf");
    sf_putint(sf_p0, "n1", nz);
    sf_putint(sf_p0, "n2", nx);
    sf_floatwrite(const_cast<float*>(&p0[0]), nz * nx, sf_p0);

    sf_file sf_p1 = sf_output("pp1.rsf");
    sf_putint(sf_p1, "n1", nz);
    sf_putint(sf_p1, "n2", nx);
    sf_floatwrite(const_cast<float*>(&p1[0]), nz * nx, sf_p1);
    exit(1);
    */

    stepForward(p0,p1);
    std::swap(p1, p0);
    recordSeis(&dcal[it*ng], &p0[0]);
	}
}



void ForwardModeling::BornForwardModeling(const std::vector<float> &exvel_m, const std::vector<float>& encSrc,
    std::vector<float>& dcal, int shot_id) const {
  int nx = getnx();
  int nz = getnz();
  int ns = getns();
  int ng = getng();

  std::vector<float> fullwv(3 * nz * nx, 0);
  std::vector<float> p0(nz * nx, 0);
  std::vector<float> p1(nz * nx, 0);
  std::vector<float> rp0(nz * nx, 0);
  std::vector<float> rp1(nz * nx, 0);
	float *fullwv_t0, *fullwv_t1, *fullwv_t2, *fullwv_t;	
	fullwv_t0 = &fullwv[0];
	fullwv_t1 = &fullwv[nz * nx];
	fullwv_t2 = &fullwv[2 * nz * nx];

  ShotPosition curSrcPos = allSrcPos->clipRange(shot_id, shot_id);
	int it = 0;
	for(int it0 = 0 ; it0 < nt + 1 ; it0 ++) {
		addSource(&p1[0], &encSrc[it0], curSrcPos);
		stepForward(p0,p1);
		std::swap(p1, p0);
		swap3(fullwv_t0, fullwv_t1, fullwv_t2);
		std::copy(p0.begin(), p0.end(), fullwv_t2);

		it = it0 - 1;
		if(it < 0) 
			continue;
		addBornwv(fullwv_t0, fullwv_t1, fullwv_t2, &exvel_m[0], dt, it, &rp1[0]);
		//fmMethod.addSource(&p1[0], &wlt[it], curSrcPos);
		stepForward(rp0,rp1);
		std::swap(rp1, rp0);
		recordSeis(&dcal[it*ng], &rp0[0]);
	}
}

void ForwardModeling::EssForwardModeling(const std::vector<float>& encSrc,
    std::vector<float>& dcal) const {
  int nx = getnx();
  int nz = getnz();
  int ns = getns();
  int ng = getng();

  std::vector<float> p0(nz * nx, 0);
  std::vector<float> p1(nz * nx, 0);

  for(int it=0; it<nt; it++) {
    addEncodedSource(&p1[0], &encSrc[it * ns]);
    stepForward(p0,p1);
    std::swap(p1, p0);
    recordSeis(&dcal[it*ng], &p0[0]);
  }
}

void ForwardModeling::bornScaleGradient(float* grad, int H) const {
  int nxpad = vel->nx;
  int nzpad = vel->nz;
  int nz = nzpad - bz0 - bzn;

#ifdef USE_OPENMP
	#pragma omp parallel for
#endif
	for (int h = -H ; h < H ; h ++) {
		int ind = h + H;
		for (int ix = bx0; ix < nxpad - bxn; ix++) {
			for (int iz = 1; iz < nz; iz++) {
				grad[ind * nxpad * nzpad + ix*nzpad + iz+bz0] *= std::sqrt(static_cast<float>(iz));
			}
		}
	}
}

void ForwardModeling::scaleGradient(float* grad) const {
  int nxpad = vel->nx;
  int nzpad = vel->nz;
  int nz = nzpad - bz0 - bzn;

  for (int ix = bx0; ix < nxpad - bxn; ix++) {
    for (int iz = 1; iz < nz; iz++) {
      grad[ix*nzpad + iz+bz0] *= std::sqrt(static_cast<float>(iz));
    }
  }
}

void ForwardModeling::removeDirectArrival(const ShotPosition &allSrcPos, const ShotPosition &allGeoPos, float* data, int nt, float t_width) const {
  int half_len = t_width / dt;
  int sx = allSrcPos.getx(0) + bx0;
  int sz = allSrcPos.getz(0) + bz0;
  int gz = allGeoPos.getz(0) + bz0; // better to assume all receivers are located at the same depth

  float vel_average = 0.0;
  int gmin = (sz < gz) ? sz : gz;
  int gmax = (sz > gz) ? sz : gz;

//  printf("dt %f, half_len %d, sx %d, selav %d, gelav %d\n", dt, half_len, sx, sz, gz);
//  printf("gmin %d, gmax %d\n", gmin, gmax);

  const std::vector<float> &vv = this->vel->dat;
  int nx = this->vel->nx;
  int nz = this->vel->nz;
  for (int i = 0; i < nx; i ++) {
    for (int k = gmin; k <= gmax; k ++) {
      vel_average += vv[i * nz + k];
    }
  }
  vel_average /= nx * (gmax - gmin + 1);

  //printf("vel_average: %.20f\n", vel_average);
  //exit(1);

  int ng = allGeoPos.ns;

  for (int itr = 0; itr < ng; itr ++) {
    int gx = allGeoPos.getx(itr) + bx0;
    int gz = allGeoPos.getz(itr) + bz0;

    float dist = (gx-sx)*(gx-sx) + (gz-sz)*(gz-sz);
    int t = (int)sqrt(dist * vel_average);
    int start = t;
    int end = ((t + 2 * half_len) > nt) ? nt : (t + 2 * half_len);

    for (int j = start; j < end; j ++) {
      data[itr * nt + j] = 0.f;
    }
  }

}

ForwardModeling::ForwardModeling(const ShotPosition& _allSrcPos, const ShotPosition& _allGeoPos,
    float _dt, float _dx, float _fm, int _nb, int _nt, int _freeSurface) :
      vel(NULL),vel_real(NULL), bcoff(NULL), allSrcPos(&_allSrcPos), allGeoPos(&_allGeoPos),
      dt(_dt), dx(_dx), fm(_fm),  nt(_nt), freeSurface(_freeSurface)
{
	if(freeSurface)
		bz0 = EXFDBNDRYLEN;
	else
		bz0 = _nb + EXFDBNDRYLEN;
  bx0 = bxn = bzn = _nb + EXFDBNDRYLEN;
  bndr.resize(bx0);
	spng = new Sponge();
  spng->initbndr(bndr.size());

	const int MAXCPML = 2;
	cpml = new CPML*[MAXCPML];
	for(int i = 0 ; i < MAXCPML ; i ++)
		cpml[i] = new CPML();
}

void ForwardModeling::addEncodedSource(float* p, const float* encsrc) const {
  this->addSource(p, encsrc, *this->allSrcPos);
}

void ForwardModeling::subEncodedSource(float* p, const float* source) const {
  this->subSource(p, source, *this->allSrcPos);
}

void ForwardModeling::recordSeis(float* seis_it, const float* p) const {
  this->recordSeis(seis_it, p, *this->allGeoPos);
}

void ForwardModeling::fwiRemoveDirectArrival(float* data, int shot_id) const {
  float t_width = 1.5 / fm;
  ShotPosition curSrcPos = allSrcPos->clipRange(shot_id, shot_id);
  this->removeDirectArrival(curSrcPos, *this->allGeoPos, data, nt, t_width);
}

void ForwardModeling::removeDirectArrival(float* data) const {
  float t_width = 1.5 / fm;
  this->removeDirectArrival(*this->allSrcPos, *this->allGeoPos, data, nt, t_width);
}

void ForwardModeling::addSource(float* p, const float* source, int is) const {

}

int ForwardModeling::getns() const {
  return allSrcPos->ns;
}

int ForwardModeling::getng() const {
  return allGeoPos->ns;
}

float ForwardModeling::getdt() const {
  return dt;
}

float ForwardModeling::getdx() const {
  return dx;
}

int ForwardModeling::getnt() const {
  return nt;
}

const ShotPosition& ForwardModeling::getAllSrcPos() const {
  return *allSrcPos;
}

const ShotPosition& ForwardModeling::getAllGeoPos() const {
  return *allGeoPos;
}

int ForwardModeling::getnx() const {
  assert(vel != NULL);
  return vel->nx;
}

int ForwardModeling::getnz() const {
  assert(vel != NULL);
  return vel->nz;
}

int ForwardModeling::getbx0() const {
  assert(vel != NULL);
  return bx0;
}

int ForwardModeling::getbxn() const {
  assert(vel != NULL);
  return bxn;
}

int ForwardModeling::getbz0() const {
  assert(vel != NULL);
  return bz0;
}

int ForwardModeling::getbzn() const {
  assert(vel != NULL);
  return bzn;
}

int ForwardModeling::getFDLEN() const {
	return EXFDBNDRYLEN;
}

Velocity& ForwardModeling::getVelocity() {
  return *const_cast<Velocity *>(vel);
}

const std::vector<float> ForwardModeling::getVelocityDiff() const
{
	std::vector<float> vel_m(vel->nx * vel->nz, 0.0f);
	vectorMinus(vel_real->dat, vel->dat, vel_m);
	return vel_m;
}

std::vector<float> ForwardModeling::initBndryVector(int nt) const {
  if (vel == NULL) {
    ERROR() << __PRETTY_FUNCTION__ << ": you should bind velocity first";
    exit(1);
  }
  int nxpad = vel->nx;
  int nzpad = vel->nz;
//  int nx = nxpad - 2*nb;
//  int nz = nzpad - nb;

//  int nx = nxpad - bx0 - bxn;
//  int nz = nzpad - bz0 - bzn;
//
//  bndrLen = FDLEN + 1;
//  bndrSize = bndrLen * (
//             nx +   /* bottom */
//             2 * nz /* left + right */
//             );

  bndrWidth = 6;
  int nx = nxpad - (bx0 - bndrWidth + bxn - bndrWidth);
  int nz = nzpad - bz0 - bzn;

  bndrSize = bndrWidth * (
      nx +   /* bottom */
      2 * nz /* left + right */
  );

  return std::vector<float>(nt*bndrSize, 0);
}

void ForwardModeling::writeBndry(float* _bndr, const float* p, int it) const {
  /**
     * say the FDLEN = 2, then the boundary we should save is mark by (*)
     * we omit the upper layer
     *
     *    **+-------------+**
     *    **-             -**
     *    **-             -**
     *    **-             -**
     *    **-             -**
     *    **-             -**
     *    **+-------------+**
     *    *******************
     *    *******************
     *
     */
    int nxpad = vel->nx;
    int nzpad = vel->nz;

    int nx = nxpad - (bx0 - bndrWidth + bxn - bndrWidth);
    int nz = nzpad - bz0 - bzn;

    float *bndr = &_bndr[it * bndrSize];

    for (int ix = 0; ix < nx; ix++) {
      for(int iz = 0; iz < bndrWidth; iz++) {
        bndr[iz + bndrWidth*ix] = p[(ix+bx0-bndrWidth)*nzpad + (nzpad - bzn + iz)]; // bottom
      }
    }

    for (int iz = 0; iz < nz; iz++) {
      for(int ix=0; ix < bndrWidth; ix++) {
        bndr[bndrWidth*nx+iz+nz*ix]         = p[(bx0-bndrWidth + ix)*nzpad + (bz0 + iz)];   // left
        bndr[bndrWidth*nx+iz+nz*(ix+bndrWidth)] = p[(nxpad - bxn + ix)*nzpad + (bz0 + iz)];  // right
      }
    }
}

void ForwardModeling::readBndry(const float* _bndr, float* p, int it) const {
  int nxpad = vel->nx;
  int nzpad = vel->nz;

  int nx = nxpad - (bx0 - bndrWidth + bxn - bndrWidth);
  int nz = nzpad - bz0 - bzn;
  const float *bndr = &_bndr[it * bndrSize];

  for (int ix = 0; ix < nx; ix++) {
    for(int iz = 0; iz < bndrWidth; iz++) {
      p[(ix+bx0-bndrWidth)*nzpad + (nzpad - bzn + iz)] = bndr[iz + bndrWidth*ix]; // bottom
    }
  }

  for (int iz = 0; iz < nz; iz++) {
    for(int ix=0; ix < bndrWidth; ix++) {
      p[(bx0-bndrWidth + ix)*nzpad + (bz0 + iz)] = bndr[bndrWidth*nx+iz+nz*ix];   // left
      p[(nxpad - bxn + ix)*nzpad + (bz0 + iz)] = bndr[bndrWidth*nx+iz+nz*(ix+bndrWidth)];  // right
    }
  }
}
