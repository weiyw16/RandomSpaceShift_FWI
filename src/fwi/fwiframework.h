/*
 * essfwiframework.h
 *
 *  Created on: Nov 22, 2016
 *      Author: cbw
 */

#ifndef SRC_FWI2D_FWIFRAMEWORK_H_
#define SRC_FWI2D_FWIFRAMEWORK_H_

#include "forwardmodeling.h"
#include "fwibase.h"
#include "fwiupdatevelop.h"
#include "fwiupdatesteplenop.h"
#include "random-code.h"

class FwiFramework : public FwiBase {
public:
  FwiFramework(ForwardModeling &fmMethod, const FwiUpdateSteplenOp &updateSteplenOp,
                  const FwiUpdateVelOp &updateVelOp, const std::vector<float> &wlt,
                  const std::vector<float> &dobs);
	void epoch(int iter);
	void calgradient(const ForwardModeling &fmMethod,
    const std::vector<float> &encSrc,
    const std::vector<float> &vsrc,
    std::vector<float> &g0,
    int nt, float dt,
		int shot_id, int rank);


protected:
  FwiUpdateSteplenOp updateStenlelOp;
  const FwiUpdateVelOp &updateVelOp;
};

#endif /* SRC_ESS_FWI2D_ESSFWIFRAMEWORK_H_ */
