/*
 * fd4t10s-damp-zjh-born.cpp
 *
 *  Created on: Oct 7, 2017
 *      Author: wyw
 */

#include "fd4t10s-zjh-born.h"

void fd4t10s_zjh_born(float *prev_wave, const float *curr_wave, const float *born_coff, int nx, int nz) {
  float a[6];

  const int d = 6;
  int ix, iz;

  /// Zhang, Jinhai's method
  a[0] = +1.53400796;
  a[1] = +1.78858721;
  a[2] = -0.31660756;
  a[3] = +0.07612173;
  a[4] = -0.01626042;
  a[5] = +0.00216736;

  //printf("fm 1\n");
#ifdef USE_OPENMP
  #pragma omp parallel for default(shared) private(ix, iz)
#endif


#ifdef USE_SWACC
/*#pragma acc parallel loop copyin(nx,nz,d) annotate(readonly=(nx,nz,d))*/
#pragma acc parallel loop
#endif
/* First, get the second sources */
  for (ix = d - 1; ix < nx - (d - 1); ix++) {
    for (iz = d - 1; iz < nz - (d - 1); iz++) {
      int curPos = ix * nz + iz;
      prev_wave[curPos] += born_coff[curPos] * ( -4.0 * a[0] * curr_wave[curPos] +
                   a[1] * (curr_wave[curPos - 1]  +  curr_wave[curPos + 1]  +
                           curr_wave[curPos - nz]  +  curr_wave[curPos + nz])  +
                   a[2] * (curr_wave[curPos - 2]  +  curr_wave[curPos + 2]  +
                           curr_wave[curPos - 2 * nz]  +  curr_wave[curPos + 2 * nz])  +
                   a[3] * (curr_wave[curPos - 3]  +  curr_wave[curPos + 3]  +
                           curr_wave[curPos - 3 * nz]  +  curr_wave[curPos + 3 * nz])  +
                   a[4] * (curr_wave[curPos - 4]  +  curr_wave[curPos + 4]  +
                           curr_wave[curPos - 4 * nz]  +  curr_wave[curPos + 4 * nz])  +
                   a[5] * (curr_wave[curPos - 5]  +  curr_wave[curPos + 5]  +
                           curr_wave[curPos - 5 * nz]  +  curr_wave[curPos + 5 * nz]) );
    }
  }

}
