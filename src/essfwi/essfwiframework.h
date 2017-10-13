/*
 * essfwiframework.h
 *
 *  Created on: Mar 10, 2016
 *      Author: rice
 */

#ifndef SRC_ESS_FWI2D_ESSFWIFRAMEWORK_H_
#define SRC_ESS_FWI2D_ESSFWIFRAMEWORK_H_

#include "forwardmodeling.h"
#include "fwibase.h"
#include "updatevelop.h"
#include "updatesteplenop.h"
#include "random-code.h"

class EssFwiFramework : public FwiBase {
public:
  EssFwiFramework(ForwardModeling &fmMethod, const UpdateSteplenOp &updateSteplenOp,
                  const UpdateVelOp &updateVelOp, const std::vector<float> &wlt,
                  const std::vector<float> &dobs);

  void epoch(int iter, float lambdaX = 0, float lambdaZ = 0, float fhi = 0);
	void calgradient(const ForwardModeling &fmMethod, const std::vector<float> &encSrc,
    const std::vector<float> &vsrc,
    std::vector<float> &g0,
    int nt, float dt);

private:
  static const int ESS_SEED = 1;

private:
  UpdateSteplenOp updateStenlelOp;
  const UpdateVelOp &updateVelOp;
  RandomCodes essRandomCodes;
};

#endif /* SRC_ESS_FWI2D_ESSFWIFRAMEWORK_H_ */
