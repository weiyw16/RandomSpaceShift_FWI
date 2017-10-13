/*
 * forwardmodeling.h
 *
 *  Created on: Feb 29, 2016
 *      Author: rice
 */

#ifndef SRC_FM2D_DAMP4T10D_H_
#define SRC_FM2D_DAMP4T10D_H_

extern "C" {
#include <rsf.h>
#include "fdutil.h"
}
#include <boost/function.hpp>
#include "velocity.h"
#include "shot-position.h"
#include "sponge.h"
#include "cpml.h"

class ForwardModeling {
public:
  ForwardModeling(const ShotPosition &allSrcPos, const ShotPosition &allGeoPos, float dt, float dx, float fm, int nb, int nt, int freeSurface);
	~ForwardModeling();

  Velocity expandDomain(const Velocity &vel);
  Velocity expandDomain_notrans(const Velocity &vel);


	void addBornwv(float *fullwv_t0, float *fullwv_t1, float *fullwv_t2, const float *exvel_m, float dt, int it, float *rp1) const;
  void stepForward(std::vector<float> &p0, std::vector<float> &p1) const;
  void swStepForward(std::vector<float> &p0, std::vector<float> &p1) const;
  void stepbornForward(std::vector<float> &p0, std::vector<float> &p1) const;
  void stepForward(std::vector<float> &p0, std::vector<float> &p1, bool vtrans) const;
  void stepForward(std::vector<float> &p0, std::vector<float> &p1, int cpmlId) const;
  void stepBackward(std::vector<float> &p0, std::vector<float> &p1) const;
  void bindVelocity(const Velocity &_vel);
  void bindRealVelocity(const Velocity &_vel);
  void bindBornCoff(std::vector<float> &b);
  void addSource(float *p, const float *source, const ShotPosition &pos) const;
  void addSource(float *p, const float *source, int is) const;
  void subSource(float *p, const float *source, const ShotPosition &pos) const;
  void addEncodedSource(float *p, const float *encsrc) const;
  void recordSeis(float *seis_it, const float *p) const;
  void bornMaskGradient(float *grad, int H) const;
  void bornScaleGradient(float *grad, int H) const;
  void maskGradient(float *grad) const;
  void scaleGradient(float *grad) const;
  void refillBoundary(float *vel) const;
  void sfWriteVel(const std::vector<float> &exvel, sf_file file) const;

  void fwiRemoveDirectArrival(float* data, int shot_id) const;
  void removeDirectArrival(float* data) const;
  void subEncodedSource(float *p, const float *source) const;
  void refillVelStencilBndry();

  std::vector<float> initBndryVector(int nt) const;
  std::vector<float> getBornCoff(const Velocity &localvel, const Velocity &localvel_real, float dx, float dt);
  void writeBndry(float* _bndr, const float* p, int it) const;
  void readBndry(const float* _bndr, float* p, int it) const;

  void FwiForwardModeling(const std::vector<float> &encsrc, std::vector<float> &dcal, int shot_id) const;
  void EssForwardModeling(const std::vector<float> &encsrc, std::vector<float> &dcal) const;
	void BornForwardModeling(const std::vector<float>& exvel, const std::vector<float>& encSrc, std::vector<float>& dcal, int shot_id) const;

public:
  const Velocity &getVelocity() const;
  Velocity &getVelocity();
	const std::vector<float> getVelocityDiff() const;
  const ShotPosition &getAllSrcPos() const;
  const ShotPosition &getAllGeoPos() const;
  int getns() const;
  int getng() const;
  float getdt() const;
  float getdx() const;
  int getnt() const;
  int getnx() const;
  int getnz() const;
  int getbx0() const;
  int getbxn() const;
  int getbz0() const;
  int getbzn() const;
	int getFDLEN() const;

private:
  void manipSource(float *p, const float *source, const ShotPosition &pos, boost::function2<float, float, float> op) const;
  void recordSeis(float *seis_it, const float *p, const ShotPosition &geoPos) const;
  void removeDirectArrival(const ShotPosition &allSrcPos, const ShotPosition &allGeoPos, float* data, int nt, float t_width) const;

public:
	CPML* getCPML(int cpmlId) const;
	void initFdUtil(sf_file &vinit, Velocity *v, int nb, float dx, float dt);

private:
  const static int EXFDBNDRYLEN = 6;

private:
  const Velocity *vel;
  const Velocity *vel_real;
  const ShotPosition *allSrcPos;
  const ShotPosition *allGeoPos;
  float dt;
  float dx;
  float fm;
  int bx0, bxn;
  int bz0, bzn;
  int nt;
	int freeSurface;	//free surface
  mutable int bndrSize;
  mutable int bndrWidth;


private:
	std::vector<float> bndr;
	std::vector<float> bcoff;
	mutable Sponge *spng;
	mutable CPML **cpml;

	struct fdm2 *fd;
	struct spon *sp;
	struct abc2 *abc;
};

#endif /* SRC_FM2D_DAMP4T10D_H_ */
