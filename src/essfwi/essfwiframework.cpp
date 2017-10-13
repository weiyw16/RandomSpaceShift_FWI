/*
 * essfwiframework.cpp
 *
 *  Created on: Mar 10, 2016
 *      Author: rice
 */


extern "C"
{
#include <rsf.h>
#include <triangle.h>
}

#include <time.h>

#include <cmath>
#include <algorithm>
#include <iterator>
#include <numeric>
#include <cstdlib>
#include <functional>
#include <vector>
#include <set>

#include "logger.h"
#include "common.h"
#include "ricker-wavelet.h"
#include "sum.h"
#include "sf-velocity-reader.h"
#include "shotdata-reader.h"
#include "random-code.h"
#include "encoder.h"
#include "velocity.h"
#include "sfutil.h"
#include "parabola-vertex.h"
#include "essfwiframework.h"

#include "aux.h"
#include "ReguFactor.h"

EssFwiFramework::EssFwiFramework(ForwardModeling &method, const UpdateSteplenOp &updateSteplenOp,
    const UpdateVelOp &_updateVelOp,
    const std::vector<float> &_wlt, const std::vector<float> &_dobs) :
    FwiBase(method, _wlt, _dobs), updateStenlelOp(updateSteplenOp), updateVelOp(_updateVelOp), essRandomCodes(ESS_SEED)
{
}

void EssFwiFramework::epoch(int iter, float lambdaX, float lambdaZ, float fhi) {
  // create random codes
  const std::vector<int> encodes = essRandomCodes.genPlus1Minus1(ns);

  std::stringstream ss;
  std::copy(encodes.begin(), encodes.end(), std::ostream_iterator<int>(ss, " "));
  DEBUG() << "code is: " << ss.str();

  Encoder encoder(encodes);
  std::vector<float> encsrc  = encoder.encodeSource(wlt);
  std::vector<float> encobs_trans = encoder.encodeObsData(dobs, nt, ng);
  std::vector<float> encobs(nt * ng, 0);
	matrix_transpose(&encobs_trans[0], &encobs[0], ng, nt);
  //cbw
  int flo = 0;
  bool verb = false;
  bool phase = false;
	if(flo != -1 && fhi != -1) 
		filter(&encobs[0], nt, dt, flo, fhi, phase, verb, ng, 1);

  std::vector<float> dcal_trans(nt * ng, 0);
  std::vector<float> dcal(nt * ng, 0);
  fmMethod.EssForwardModeling(encsrc, dcal_trans);
	matrix_transpose(&dcal_trans[0], &dcal[0], ng, nt);
  //cbw
	if(flo != -1 && fhi != -1) 
		filter(&dcal[0], nt, dt, flo, fhi, phase, verb, ng, 1);

  //cbw
  /*
	bool phase = false;
	bool verb = false;
	if(fhi != -1) 
		filter(&dcal[0], nt, dt, 0, fhi, phase, verb, ng, 1);
    */

  fmMethod.removeDirectArrival(&encobs[0]);
  fmMethod.removeDirectArrival(&dcal[0]);

  std::vector<float> vsrc(nt * ng, 0);
  vectorMinus(encobs, dcal, vsrc);
  float obj1 = cal_objective(&vsrc[0], vsrc.size());
  initobj = iter == 0 ? obj1 : initobj;
  DEBUG() << format("obj: %e") % obj1;

    Velocity &exvel = fmMethod.getVelocity();
    if (!(lambdaX == 0 && lambdaZ == 0)) {
      ReguFactor fac(&exvel.dat[0], nx, nz, lambdaX, lambdaZ);
      obj1 += fac.getReguTerm();
    }


  transVsrc(vsrc, nt, ng);

  std::vector<float> g1(nx * nz, 0);
  calgradient(fmMethod, encsrc, vsrc, g1, nt, dt);

  DEBUG() << format("grad %.20f") % sum(g1);

  fmMethod.scaleGradient(&g1[0]);
  fmMethod.maskGradient(&g1[0]);

  //sf_file sf_g = sf_output("g1.rsf");
  //sf_putint(sf_g, "n1", nz);
  //sf_putint(sf_g, "n2", nx);
  //sf_floatwrite(&g1[0], nx * nz, sf_g);
//#define SMOOTH
#ifdef SMOOTH 
  INFO() << "smoothing";
  if(fhi != -1)
  {
    //int rectz = 10 - fhi;
    int rectz = 5.0 / fhi;
    int rectx = rectz;
    rectz = rectz > 1 ? rectz : 1;
    rectx = rectx > 1 ? rectx : 1;
    smooth(&g1[0], nz, nx, rectz, rectx);
  }
#endif

  updateGrad(&g0[0], &g1[0], &updateDirection[0], g0.size(), iter);

  updateStenlelOp.bindEncSrcObs(encsrc, encobs);
  float steplen;
  updateStenlelOp.calsteplen(updateDirection, obj1, iter, lambdaX, lambdaZ, steplen, updateobj);

//  Velocity &exvel = fmMethod.getVelocity();
  updateVelOp.update(exvel, exvel, updateDirection, steplen);

  fmMethod.refillBoundary(&exvel.dat[0]);
}

void EssFwiFramework::calgradient(const ForwardModeling &fmMethod,
    const std::vector<float> &encSrc,
    const std::vector<float> &vsrc,
    std::vector<float> &g0,
    int nt, float dt)
{
  int nx = fmMethod.getnx();
  int nz = fmMethod.getnz();
  int ns = fmMethod.getns();
  int ng = fmMethod.getng();
  const ShotPosition &allGeoPos = fmMethod.getAllGeoPos();
  const ShotPosition &allSrcPos = fmMethod.getAllSrcPos();

  std::vector<float> bndr = fmMethod.initBndryVector(nt);
  std::vector<float> sp0(nz * nx, 0);
  std::vector<float> sp1(nz * nx, 0);
  std::vector<float> gp0(nz * nx, 0);
  std::vector<float> gp1(nz * nx, 0);


  for(int it=0; it<nt; it++) {
    fmMethod.addSource(&sp1[0], &encSrc[it * ns], allSrcPos);
    fmMethod.stepForward(sp0,sp1);
    std::swap(sp1, sp0);
    fmMethod.writeBndry(&bndr[0], &sp0[0], it);
  }

  std::vector<float> vsrc_trans(nt * ng, 0);
  matrix_transpose(const_cast<float*>(&vsrc[0]), &vsrc_trans[0], nt, ng);

  for(int it = nt - 1; it >= 0 ; it--) {
    fmMethod.readBndry(&bndr[0], &sp0[0], it);
    std::swap(sp0, sp1);
    fmMethod.stepBackward(sp0, sp1);
    fmMethod.subEncodedSource(&sp0[0], &encSrc[it * ns]);

    /**
     * forward propagate receviers
     */
    fmMethod.addSource(&gp1[0], &vsrc_trans[it * ng], allGeoPos);
    fmMethod.stepForward(gp0,gp1);
    std::swap(gp1, gp0);

    if (dt * it > 0.4) {
      cross_correlation(&sp0[0], &gp0[0], &g0[0], g0.size(), 1.0);
    } else if (dt * it > 0.3) {
      cross_correlation(&sp0[0], &gp0[0], &g0[0], g0.size(), (dt * it - 0.3) / 0.1);
    } else {
      break;
    }
 }
}

