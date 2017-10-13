#ifndef FD_SW_H
#define FD_SW_H
#include <simd.h>
#define SIMDType floatv4 
#include "fd4t10s-damp-zjh-cg.h"

void fd4t10s_sw(struct Param *param);
#endif
