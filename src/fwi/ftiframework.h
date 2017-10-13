/*
 * ftiframework.h
 *
 *  Created on: Nov 22, 2016
 *      Author: cbw
 */

#ifndef SRC_FWI2D_FTIFRAMEWORK_H_
#define SRC_FWI2D_FTIFRAMEWORK_H_

#include "forwardmodeling.h"
#include "fwiupdatevelop.h"
#include "fwiupdatesteplenop.h"
#include "random-code.h"
#include "fwiframework.h"

class FtiFramework: public FwiFramework {
public:
  FtiFramework(ForwardModeling &fmMethod, const FwiUpdateSteplenOp &updateSteplenOp,
                  const FwiUpdateVelOp &updateVelOp, const std::vector<float> &wlt,
                  const std::vector<float> &dobs, int jsx, int jsz);
  void epoch(int iter);
	void calgradient(const ForwardModeling &fmMethod,
    const std::vector<float> &encSrc,
    const std::vector<float> &vsrc,
    std::vector<float> &I,
    std::vector<float> &gd,
    int nt, float dt,
		int shot_id, int rank, int H);
	void image_born(const ForwardModeling &fmMethod, const std::vector<float> &encSrc,
    const std::vector<float> &vsrc,
    std::vector<float> &g0,
    int nt, float dt,
		int shot_id, int rank, int H);
	void cross_correlation(float *src_wave, float *vsrc_wave, float *image, int nx, int nz, float scale, int H);

	int jsx, jsz;
};

#endif /* SRC_ESS_FWI2D_ESSFWIFRAMEWORK_H_ */
