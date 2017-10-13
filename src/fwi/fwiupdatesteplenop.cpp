/*
 * ipdatesteplenop.cpp
 *
 *  Created on: Mar 14, 2016
 *      Author: rice
 */

#include <set>
#include <cmath>
#include "fwiupdatesteplenop.h"
#include "logger.h"
#include "common.h"
#include "parabola-vertex.h"
#include "sum.h"
#include "mpi.h"

namespace {
typedef std::pair<float, float> ParaPoint;

bool parabolicLessComp(const ParaPoint &a, const ParaPoint &b) {
  return a.second - b.second < 1e-10;
}



void calMaxAlpha2_3(const Velocity &exvel,  const float *grad, float dt, float dx, float maxdv,
                    float &ret_alpha2, float &ret_alpha3) {
  const int nx = exvel.nx;
  const int nz = exvel.nz;

  const std::vector<float> &vel = exvel.dat;
  float alpha2 = FLT_MAX;
  for (int i = 0; i < nx * nz; i++) {
    float tmpv = dx / (dt * std::sqrt(vel[i]));
    tmpv -= maxdv;
    tmpv = (dx / (dt * tmpv)) * (dx / (dt * tmpv));
    if (std::fabs(grad[i]) < 1e-10 ) {
      continue;
    }
    if (alpha2 > (tmpv - vel[i]) / std::fabs(grad[i])) {
      alpha2 = (tmpv - vel[i]) / std::fabs(grad[i]);
    }
  }

  /// return the value
  ret_alpha2 = alpha2;
  ret_alpha3 =  2 * alpha2;
}

} /// end of name space

FwiUpdateSteplenOp::FwiUpdateSteplenOp(const ForwardModeling &fmMethod, const FwiUpdateVelOp &updateVelOp,
    int max_iter_select_alpha3, float maxdv, int ns, int ng, int nt, std::vector<float> *encsrc) :
  fmMethod(fmMethod), updateVelOp(updateVelOp), encsrc(encsrc), encobs(NULL),
  max_iter_select_alpha3(max_iter_select_alpha3), maxdv(maxdv), ns(ns), ng(ng), nt(nt)
{

}

float FwiUpdateSteplenOp::calobjval(const std::vector<float>& grad,
    float steplen, int shot_id) const {
  int nx = fmMethod.getnx();
  int nz = fmMethod.getnz();
  int nt = fmMethod.getnt();

  const Velocity &oldVel = fmMethod.getVelocity();
  Velocity newVel(nx, nz);

  /*
	sf_file sf_oldvel = sf_output("oldvel_before.rsf");
	sf_putint(sf_oldvel, "n1", nz);
	sf_putint(sf_oldvel, "n2", nx);
	sf_floatwrite(const_cast<float*>(&oldVel.dat[0]), nx * nz, sf_oldvel);

	sf_file sf_grad = sf_output("grad_before.rsf");
	sf_putint(sf_grad, "n1", nz);
	sf_putint(sf_grad, "n2", nx);
	sf_floatwrite(const_cast<float*>(&grad[0]), nx * nz, sf_grad);
  */

  updateVelOp.update(newVel, oldVel, grad, steplen);

  /*
	sf_file sf_newvel = sf_output("newvel_after.rsf");
	sf_putint(sf_newvel, "n1", nz);
	sf_putint(sf_newvel, "n2", nx);
	sf_floatwrite(&newVel.dat[0], nx * nz, sf_newvel);
  exit(1);
  */

  ForwardModeling *updateMethod = const_cast<ForwardModeling*>(&fmMethod);
  updateMethod->bindVelocity(newVel);

  //forward modeling
  int ng = fmMethod.getng();
  std::vector<float> dcal(nt * ng);
  std::vector<float> dcal_trans(nt * ng);
  updateMethod->FwiForwardModeling(*encsrc, dcal_trans, shot_id);
	matrix_transpose(&dcal_trans[0], &dcal[0], ng, nt);

  /*
	sf_file sf_dcal2 = sf_output("dcal2.rsf");
	sf_putint(sf_dcal2, "n1", nt);
	sf_putint(sf_dcal2, "n2", ng);
	sf_floatwrite(&dcal[0], nt * ng, sf_dcal2);
  exit(1);
  */

  updateMethod->bindVelocity(oldVel);  //-test
  updateMethod->fwiRemoveDirectArrival(&dcal[0], shot_id);
  //updateMethod.bindVelocity(newVel);  //-test

	/*
	sf_file sf_dcal3 = sf_output("dcal3.rsf");
	sf_putint(sf_dcal3, "n1", nt);
	sf_putint(sf_dcal3, "n2", ng);
	sf_floatwrite(&dcal[0], nt * ng, sf_dcal3);
  exit(1);
	*/

  std::vector<float> vdiff(nt * ng, 0);
		INFO() << "****sum encobs: " << std::accumulate((*encobs).begin(), (*encobs).begin() + ng * nt, 0.0f);
		INFO() << "****sum2 dcal: " << std::accumulate(dcal.begin(), dcal.begin() + ng * nt, 0.0f);
	
  vectorMinus(*encobs, dcal, vdiff);
  float val = cal_objective(&vdiff[0], vdiff.size());

  DEBUG() << format("curr_alpha = %e, pure object value = %e") % steplen % val;

  return val;
}

bool FwiUpdateSteplenOp::refineAlpha(const std::vector<float> &grad, float obj_val1, float maxAlpha3,
    float& _alpha2, float& _obj_val2, float& _alpha3, float& _obj_val3, int shot_id) const {

  TRACE() << "SELECTING THE RIGHT OBJECTIVE VALUE 3";

  float alpha3 = _alpha3;
  float alpha2 = _alpha2;
  float obj_val2, obj_val3 = 0;

	fmMethod.fwiRemoveDirectArrival(&(*encobs)[0], shot_id);

  obj_val2 = calobjval(grad, alpha2, shot_id);
  obj_val3 = calobjval(grad, alpha3, shot_id);

  //DEBUG() << "BEFORE TUNNING";
  DEBUG() << __FUNCTION__ << format(" alpha1 = %e, obj_val1 = %e") % 0. % obj_val1;
  DEBUG() << __FUNCTION__ << format(" alpha2 = %e, obj_val2 = %e") % alpha2 % obj_val2;
  DEBUG() << __FUNCTION__ << format(" alpha3 = %e, obj_val3 = %e") % alpha3 % obj_val3;


	/*
  TRACE() << "maintain a set to store alpha2 that we ever tuned";
  std::set<ParaPoint, bool (*)(const ParaPoint &, const ParaPoint &) > tunedAlpha(parabolicLessComp);
  tunedAlpha.insert(std::make_pair(alpha2, obj_val2));

  DEBUG() << "BEGIN TUNING";
  /// obj_val2 might be quite large, so we should make it smaller by halfing alpha2

  int iter = 0;
  for (; iter < max_iter_select_alpha3 && obj_val2 > obj_val1; iter++) {

    /// pass the property of alpha2 to alpha3
    alpha3 = alpha2;
    obj_val3 = obj_val2;

    /// update alpha2
    alpha2 /= 2;
    obj_val2 = calobjval(grad, alpha2);

    /// store it
    tunedAlpha.insert(std::make_pair(alpha2, obj_val2));
    DEBUG() << __FUNCTION__ << format(" iter = %d, alpha2 = %e, obj_val2 = %e") % iter % alpha2 % obj_val2;
    DEBUG() << __FUNCTION__ << format(" iter = %d, alpha3 = %e, obj_val3 = %e\n") % iter % alpha3 % obj_val3;
  }

  DEBUG() << "SELECT A BETTER ALPHA2 IN " << iter << " ITERS";


  DEBUG() << "tunedAlpha size: " << tunedAlpha.size();
  for (std::set<ParaPoint, bool (*)(const ParaPoint &, const ParaPoint &) >::iterator it = tunedAlpha.begin();
      it != tunedAlpha.end(); ++it) {
    DEBUG() << format("alpha %e, obj %e") % it->first % it->second;
  }

  TRACE() << "check if we need to forward tuning";
  TRACE() << "after halfing in the previous step, obj_val2 might still be larger than obj_val1"
          "then we should stop tunting and choose a best alpha2 ever got";
  if (obj_val2 > obj_val1) {
    DEBUG() << "UNABLE TO TUNING A ALPHA2 BY HALFING";
    DEBUG() << "SELECT A BEST ALPHA2 EVER GOT";
    std::set<ParaPoint>::iterator it = tunedAlpha.begin();
    _alpha2 = it->first;
    _obj_val2 = it->second;

    _alpha3 = std::min(_alpha2 * 2, maxAlpha3);
    _obj_val3 = calobjval(grad, alpha3);


    DEBUG() << __FUNCTION__ << format(" alpha2 = %e, obj_val2 = %e") % _alpha2 % _obj_val2;
    DEBUG() << __FUNCTION__ << format(" alpha3 = %e, obj_val3 = %e") % _alpha3 % _obj_val3;

    bool toParabolicFit = false;
    return toParabolicFit;
  }

  TRACE() << "now we can make sure that obj_val2 < obj_val1";

  const float alpha1 = 0;
  float linearFitAlph3 = (obj_val2 - obj_val1) / (alpha2 - alpha1) * (alpha3 - alpha1) + obj_val1;
  DEBUG() << __FUNCTION__ << format(" linear fit alpha3 = %e ") % linearFitAlph3;

  TRACE() << "keep the alpha we tuned";
  tunedAlpha.clear();
  tunedAlpha.insert(std::make_pair(alpha3, obj_val3));

  while (obj_val3 < linearFitAlph3 && obj_val3 < obj_val1 && alpha3 < maxAlpha3) {
    TRACE() << "if in this case, we should enlarge alpha3";
    alpha2 = alpha3;
    obj_val2 = obj_val3;

    alpha3 = std::min(alpha3 * 2, maxAlpha3);
    obj_val3 = calobjval(grad, alpha3);

    tunedAlpha.insert(std::make_pair(alpha3, obj_val3));

    DEBUG() << __FUNCTION__ << format(" tune alpha3, alpha2 = %e, obj_val2 = %e") % alpha2 % obj_val2;
    DEBUG() << __FUNCTION__ << format(" tune alpha3, alpha3 = %e, obj_val3 = %e") % alpha3 % obj_val3;
  }


  TRACE() << "If we couldnot tune a good alpha3";
  if (alpha3 > maxAlpha3 + 0.1) {
    DEBUG() << "UNABLE TO TUNING A ALPHA3 BY DOUBLING";
    DEBUG() << "SELECT A BEST ALPHA3 EVER GOT";
    std::set<ParaPoint>::iterator it = tunedAlpha.begin();
    _alpha3 = it->first;
    _obj_val3 = it->second;

    _alpha2 = _alpha3 / 2;
    _obj_val2 = calobjval(grad, alpha2);

    bool toParabolicFit = false;

    DEBUG() << __FUNCTION__ << format(" alpha2 = %e, obj_val2 = %e") % _alpha2 % _obj_val2;
    DEBUG() << __FUNCTION__ << format(" alpha3 = %e, obj_val3 = %e") % _alpha3 % _obj_val3;
    return toParabolicFit;
  }

	*/
  /// return objval2 and objval3
  bool toParabolicFit = true;
  //_alpha2 = alpha2;
  //_alpha3 = alpha3;
  _obj_val2 = obj_val2;
  _obj_val3 = obj_val3;

  //DEBUG() << __FUNCTION__ << format(" alpha2 = %e, obj_val2 = %e") % _alpha2 % _obj_val2;
  //DEBUG() << __FUNCTION__ << format(" alpha3 = %e, obj_val3 = %e") % _alpha3 % _obj_val3;

  return toParabolicFit;
}

void FwiUpdateSteplenOp::parabola_fit(float x0, float x1, float x2, float y0, float y1, float y2, float max_alpha3, bool toParabolic, int iter, float &xmin, float &objval) {
  float a, b, c, k;

  k = 1. / ((((float)x0) - x1) * (x0 - x2) * (x1 - x2));
  a = k * (y0 * (x1 - x2) - y1 * (x0 - x2) + y2 * (x0 - x1));
  b =  -k * (y0 * (x1 * x1 - x2 * x2) - y1 * (x0 * x0 - x2 * x2) + y2 * (x0 * x0 - x1 * x1));
  c = k * (y0 * x1 * x2 * (x1 - x2) - y1 * x0 * x2 * (x0 - x2) + y2 * x0 * x1 * (x0 - x1));

  if (a == 0) {
    if (b > 0.) {
      xmin = x2;
    } else {
      xmin = x0 / 100;
    }
  } else if (a > 0.) {
    if ((-b / (2.*a)) < 0) {
      xmin = x1 / 100;
    } else if ((-b / (2.*a)) < (x2 / 15.) ) {
      xmin = -b / (2.*a);
    } else if ((-b / (2.*a)) > x2 ) {
      xmin = x2;
    } else {
      xmin = -b / (2.*a);
    }
  } else if (a < 0.) {
    if ((-b / (2 * a)) < 0.5 * (x0 + x2) ) {
      xmin = x2;
    }
    if ((-b / (2 * a)) > 0.5 * (x0 + x2) ) {
      xmin = x1 / 100;
    }
  }
}

/*
void FwiUpdateSteplenOp::parabola_fit(float alpha1, float alpha2, float alpha3, float obj_val1, float obj_val2, float obj_val3, float max_alpha3, bool toParabolic, int iter, float &steplen, float &objval) {
  float alpha4, obj_val4;
  if (toParabolic) {
    DEBUG() << "parabolic fit";

    TRACE() << "max_alpha3: " << max_alpha3;
    TRACE() << "alpha1: " << alpha1 << ", obj_val1: " << obj_val1 << ", alpha2: " << alpha2 << ", obj_val2: " << obj_val2
            << "alpha3: " << alpha3 << ", obj_val3: " << obj_val3;
    parabolaVertex(alpha1, obj_val1, alpha2, obj_val2, alpha3, obj_val3, max_alpha3, alpha4, obj_val4);
    TRACE() << "parabolaVertex done";

    if (alpha4 > max_alpha3) {
      DEBUG() << format("alpha4 = %e, max_alpha3 = %e") % alpha4 % max_alpha3;
      DEBUG() << format("alpha4 is greater than max_alpha3, set it to alpha3");
      alpha4 = max_alpha3;
    }
  } else {
    DEBUG() << "NO need to perform parabolic fit";
    alpha4 = alpha3;
    obj_val4 = obj_val3;
  }

  INFO() << format("iter %d  alpha  = %e total obj_val1 = %e") % iter % alpha1 % obj_val1;
  INFO() << format("iter %d  alpha2 = %e total obj_val2 = %e") % iter % alpha2 % obj_val2;
  INFO() << format("iter %d  alpha3 = %e total obj_val3 = %e") % iter % alpha3 % obj_val3;
  INFO() << format("iter %d  alpha4 = %e total obj_val4 = %e\n") % iter % alpha4 % obj_val4;

  preservedAlpha.alpha = alpha4;
  steplen = alpha4;
  objval = obj_val4;
}
*/

void FwiUpdateSteplenOp::calsteplen(const std::vector<float> &dobs, const std::vector<float>& grad,
    float obj_val1, int iter, float &steplen, float &objval, int rank, int shot_begin, int shot_end) {

  float dt = fmMethod.getdt();
  float dx = fmMethod.getdx();

  /// "calculate the initial value of alpha2 and alpha3";
  float max_alpha2, max_alpha3;

  calMaxAlpha2_3(fmMethod.getVelocity(), &grad[0], dt, dx, maxdv, max_alpha2, max_alpha3);
  DEBUG() << format("               max_alpha2 = %e,  max_alpha3: = %e") % max_alpha2 % max_alpha3;

  alpha1 = 0;
	this->obj_val1 = obj_val1;
  initAlpha23(max_alpha3, alpha2, alpha3);
  DEBUG() << format("after init alpha,  alpha2 = %e,      alpha3: = %e") % alpha2 % alpha3;
	float local_obj_val2_sum = 0.0f;
	float local_obj_val3_sum = 0.0f;
	obj_val1_sum = 0.0f;
	obj_val2_sum = 0.0f;
	obj_val3_sum = 0.0f;

	maxAlpha3 = max_alpha3;

	for(int is = shot_begin ; is < shot_end ; is ++)
	{
		std::vector<float> t_obs(ng * nt);
		std::vector<float> t_obs_trans(ng * nt);
		memcpy(&t_obs_trans[0], &dobs[is * ng * nt], sizeof(float) * ng * nt);
	  matrix_transpose(&t_obs_trans[0], &t_obs[0], ng, nt);

		encobs = &t_obs;
		INFO() << format("calculate steplen, shot id: %d") % is;
		toParabolic = refineAlpha(grad, obj_val1, max_alpha3, alpha2, obj_val2, alpha3, obj_val3, is);
		local_obj_val2_sum += obj_val2;
		local_obj_val3_sum += obj_val3;
	}
	obj_val1_sum = obj_val1;
	MPI_Allreduce(&local_obj_val2_sum, &obj_val2_sum, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&local_obj_val3_sum, &obj_val3_sum, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
	if(rank == 0)
	{
		INFO() << format("In calsteplen(): iter %d  alpha = %e total obj_val1 = %e") % iter % alpha1 % obj_val1_sum;
		INFO() << format("In calsteplen(): iter %d  alpha2 = %e total obj_val2 = %e") % iter % alpha2 % obj_val2_sum;
		INFO() << format("In calsteplen(): iter %d  alpha3 = %e total obj_val3 = %e") % iter % alpha3 % obj_val3_sum;
	}
}

void FwiUpdateSteplenOp::bindEncSrcObs(const std::vector<float>& encsrc,
    const std::vector<float>& encobs) {
  this->encsrc = &encsrc;
  //this->encobs = &encobs;
}

void FwiUpdateSteplenOp::initAlpha23(float maxAlpha3, float &initAlpha2, float &initAlpha3) {
	/*
  const float minAlpha   = 1.0E-7;
  const float resetAlpha = 1.0E-4;

  if (!preservedAlpha.init) {
    preservedAlpha.init = true;
    preservedAlpha.alpha = maxAlpha3;
  }

  initAlpha3 = preservedAlpha.alpha;
  initAlpha3 = initAlpha3 < minAlpha ? resetAlpha : initAlpha3;
  initAlpha2 = initAlpha3 * 0.5;
	*/
	initAlpha3 = maxAlpha3;
	initAlpha2 = initAlpha3 * 0.5;
}
