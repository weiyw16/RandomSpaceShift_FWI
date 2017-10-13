/*
 * fwibase.cpp
 *
 *  Created on: Nov 22, 2016
 *      Author: cbw
 */


extern "C"
{
#include <rsf.h>
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
#include "fwibase.h"

#include "aux.h"

#ifdef USE_SW
extern "C"
{
#include "fd4t10s-damp-zjh-cg.h"
}
#endif

FwiBase::FwiBase(ForwardModeling &method, const std::vector<float> &_wlt, const std::vector<float> &_dobs) :
    fmMethod(method), wlt(_wlt), dobs(_dobs),
    ns(method.getns()), ng(method.getng()), nt(method.getnt()),
    nx(method.getnx()), nz(method.getnz()), dx(method.getdx()), dt(method.getdt()),
    updateobj(0), initobj(0)
{
  g0.resize(nx*nz, 0);
  updateDirection.resize(nx*nz, 0);
}

void FwiBase::writeVel(sf_file file) const {
	fmMethod.sfWriteVel(fmMethod.getVelocity().dat, file);
}

float FwiBase::getUpdateObj() const {
	return updateobj;
}

float FwiBase::getInitObj() const {
	return initobj;
}

void FwiBase::cross_correlation(float *src_wave, float *vsrc_wave, float *image, int model_size, float scale) {
	/*
  for (int i = 0; i < model_size; i ++) {
    image[i] -= src_wave[i] * vsrc_wave[i] * scale;
  }
	*/
#ifdef USE_SW
	int nx = fmMethod.getnx();
	int nz = fmMethod.getnz();
	cross_correlation_cg(src_wave, vsrc_wave, image, nx, nz, scale);
#else
  for (int i = 0; i < model_size; i ++) {
    image[i] -= src_wave[i] * vsrc_wave[i] * scale;
  }
#endif

}

void FwiBase::transVsrc(std::vector<float> &vsrc, int nt, int ng) {
  for (int ig = 0; ig < ng; ig++) {
    second_order_virtual_source_forth_accuracy(&vsrc[ig * nt], nt);
  }
}

void FwiBase::updateGrad(float *pre_gradient, const float *cur_gradient, float *update_direction,
                           int model_size, int iter) {
  if (iter == 0) {
    std::copy(cur_gradient, cur_gradient + model_size, update_direction);
    std::copy(cur_gradient, cur_gradient + model_size, pre_gradient);
  } else {
    float beta = 0.0f;
    float a = 0.0f;
    float b = 0.0f;
    float c = 0.0f;
    int   i = 0;
    for (i = 0; i < model_size; i ++) {
      a += (cur_gradient[i] * cur_gradient[i]);
      b += (cur_gradient[i] * pre_gradient[i]);
      c += (pre_gradient[i] * pre_gradient[i]);
    }

    beta = (a - b) / c;

    if (beta < 0.0f) {
      beta = 0.0f;
    }

    for (i = 0; i < model_size; i ++) {
      update_direction[i] = cur_gradient[i] + beta * update_direction[i];
    }

    TRACE() << "Save current gradient to pre_gradient for the next iteration's computation";
    std::copy(cur_gradient, cur_gradient + model_size, pre_gradient);
  }
}

void FwiBase::one_order_virtual_source_forth_accuracy(float *vsrc, int num) {
  float *tmp_vsrc = (float *)malloc(num * sizeof(float));
  memcpy(tmp_vsrc, vsrc, num * sizeof(float));
  int i = 0;
  for (i = 0; i < num; i ++) {
    if ( i <= 1) {
      vsrc[i] = 0.0f;
      continue;
    }

    if ( (num - 1) == i || (num - 2) == i) {
      vsrc[i] = 0.0f;
      continue;
    }

    //vsrc[i] = -1. / 3 * tmp_vsrc[i - 1] -1. / 3 * tmp_vsrc[i] + 2. / 3 * tmp_vsrc[i + 1];
    vsrc[i] = (- tmp_vsrc[i - 1] + tmp_vsrc[i + 1]) / 2;
  }

  free(tmp_vsrc);
}

void FwiBase::second_order_virtual_source_forth_accuracy(float *vsrc, int num) {
  float *tmp_vsrc = (float *)malloc(num * sizeof(float));
  memcpy(tmp_vsrc, vsrc, num * sizeof(float));
  int i = 0;
  for (i = 0; i < num; i ++) {
    if ( i <= 1) {
      vsrc[i] = 0.0f;
      continue;
    }

    if ( (num - 1) == i || (num - 2) == i) {
      vsrc[i] = 0.0f;
      continue;
    }

    vsrc[i] = -1. / 12 * tmp_vsrc[i - 2] + 4. / 3 * tmp_vsrc[i - 1] -
              2.5 * tmp_vsrc[i] + 4. / 3 * tmp_vsrc[i + 1] - 1. / 12 * tmp_vsrc[i + 2];
  }

  free(tmp_vsrc);
}

