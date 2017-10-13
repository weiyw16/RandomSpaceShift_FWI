/*
 * fwibase.h
 *
 *  Created on: Nov 22, 2016
 *      Author: cbw
 */

#ifndef SRC_FWI2D_FWIBASE_H_
#define SRC_FWI2D_FWIBASE_H_

#include "forwardmodeling.h"

class FwiBase {
public:
  FwiBase(ForwardModeling &fmMethod, const std::vector<float> &wlt,
                  const std::vector<float> &dobs);
	void cross_correlation(float *src_wave, float *vsrc_wave, float *image, int model_size, float scale);
	void transVsrc(std::vector<float> &vsrc, int nt, int ng);
	void updateGrad(float *pre_gradient, const float *cur_gradient, float *update_direction, int model_size, int iter);
	void one_order_virtual_source_forth_accuracy(float *vsrc, int num);
	void second_order_virtual_source_forth_accuracy(float *vsrc, int num);
  void writeVel(sf_file file) const;
  float getUpdateObj() const;
  float getInitObj() const;

protected:
  ForwardModeling &fmMethod;
  const std::vector<float> &wlt;  /// wavelet
  const std::vector<float> &dobs; /// actual observed data (nt*ng*ns)

protected: /// propagate from other construction
  int ns;
  int ng;
  int nt;
  int nx;
  int nz;
  float dx;
  float dt;

protected:
  std::vector<float> g0;               /// gradient in previous step
  std::vector<float> updateDirection;
  float updateobj;
  float initobj;
	float obj_val4;
};

#endif /* SRC_ESS_FWI2D_ESSFWIFRAMEWORK_H_ */
