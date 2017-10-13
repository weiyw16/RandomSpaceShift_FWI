#include <slave.h>
#include <dma.h>
#include "fd4t10s-damp-zjh-cg.h"
#include "fd4t10s-slave.h"
#define SIMDSIZE 4

#define ROW 8 //change
#define COL 8

const int d = 6;

void fd4t10s_sw(struct Param *param) {
	int id = athread_get_id(-1);
	int rid = id / ROW;
	int cid = id % ROW;
	Type *prev_wave = param->prev_wave;
	Type *curr_wave = param->curr_wave;
	Type *next_wave = param->next_wave;
	Type *vel = param->vel;
	Type *u2 = param->u2;
	int nx = param->nx;
	int nz = param->nz;
	int gnxbeg = param->gnxbeg;
	int gnzbeg = param->gnzbeg;
	int ldnx = param->ldnx;
	int ldnz = param->ldnz;
	int nb = param->nb;
	int freeSurface = param->freeSurface;
	int cnx = nx - 2 * d;
	int cnz = nz - 2 * d;
	int snx = cnx / ROW + 1;
	int snz = cnz / COL + 1;
	int snxbeg = snx * cid - d + d; 
	int snzbeg = snz * rid - d + d;
	int snxend = snx * (cid + 1) + d + d;
	int snzend = snz * (rid + 1) + d + d;
	int gxbeg = gnxbeg + snxbeg;
	int gzbeg = gnzbeg + snzbeg;
	snxbeg = snxbeg > 0 ? snxbeg : 0; //can delete
	snzbeg = snzbeg > 0 ? snzbeg : 0; //can delete
	snxend = snxend < nx ? snxend : nx;
	snzend = snzend < nz ? snzend : nz;
	snx = snxend - snxbeg;
	snz = snzend - snzbeg;
	Type *local_vel = (Type *)ldm_malloc(sizeof(Type) * (snz - 2 * d) * (snx - 2 * d));
	Type *local_prev_wave = (Type *)ldm_malloc(sizeof(Type) * (snz - 2 * d) * (snx - 2 * d));
	Type *local_curr_wave = (Type *)ldm_malloc(sizeof(Type) * snz * snx);
	Type *local_u2 = (Type *)ldm_malloc(sizeof(Type) * snz * snx);
	Type *local_bndr = (Type *)ldm_malloc(sizeof(Type) * nb);

  volatile int  replyget_curr_one = 0, replyget_curr_wave = 0, replyget_prev_one = 0, replyget_prev_wave = 0, replyput_next_wave = 0, replyget_vel_one = 0, replyput_next_one = 0, replyput_curr_one = 0;
  dma_desc dma_get_curr_one, dma_get_curr_wave, dma_get_prev_one, dma_get_prev_wave, dma_put_next_wave, dma_get_vel_one, dma_put_next_one, dma_put_curr_one;

  dma_set_op(&dma_get_curr_one, DMA_GET);
  dma_set_mode(&dma_get_curr_one, PE_MODE);
  dma_set_reply(&dma_get_curr_one, &replyget_curr_one);

  dma_set_op(&dma_get_curr_wave, DMA_GET);
  dma_set_mode(&dma_get_curr_wave, PE_MODE);
  dma_set_reply(&dma_get_curr_wave, &replyget_curr_wave);

  dma_set_op(&dma_get_prev_one, DMA_GET);
  dma_set_mode(&dma_get_prev_one, PE_MODE);
  dma_set_reply(&dma_get_prev_one, &replyget_prev_one);

  dma_set_op(&dma_get_prev_wave, DMA_GET);
  dma_set_mode(&dma_get_prev_wave, PE_MODE);
  dma_set_reply(&dma_get_prev_wave, &replyget_prev_wave);

  dma_set_op(&dma_put_next_wave, DMA_PUT);
  dma_set_mode(&dma_put_next_wave, PE_MODE);
  dma_set_reply(&dma_put_next_wave, &replyput_next_wave);

  dma_set_op(&dma_get_vel_one, DMA_GET);
  dma_set_mode(&dma_get_vel_one, PE_MODE);
  dma_set_reply(&dma_get_vel_one, &replyget_vel_one);

  dma_set_op(&dma_put_next_one, DMA_PUT);
  dma_set_mode(&dma_put_next_one, PE_MODE);
  dma_set_reply(&dma_put_next_one, &replyput_next_one);

  dma_set_op(&dma_put_curr_one, DMA_PUT);
  dma_set_mode(&dma_put_curr_one, PE_MODE);
  dma_set_reply(&dma_put_curr_one, &replyput_curr_one);

	dma_set_size(&dma_get_curr_one, snz*sizeof(Type));
	dma_set_bsize(&dma_get_curr_one, snz*sizeof(Type));
	dma_set_stepsize(&dma_get_curr_one, (ldnz - snz) * sizeof(Type));

  dma_set_size(&dma_get_curr_wave, snx*snz*sizeof(Type));
  dma_set_bsize(&dma_get_curr_wave, snz*sizeof(Type));
  dma_set_stepsize(&dma_get_curr_wave, (ldnz - snz)*sizeof(Type));

  dma_set_size(&dma_get_prev_one, (snz - 2 * d) * sizeof(Type));
  dma_set_bsize(&dma_get_prev_one, (snz - 2 * d) * sizeof(Type));
  dma_set_stepsize(&dma_get_prev_one, (ldnz - (snz - 2 * d)) *sizeof(Type));

  dma_set_size(&dma_get_prev_wave, (snx - 2 * d) * (snz - 2 * d) * sizeof(Type));
  dma_set_bsize(&dma_get_prev_wave, (snz - 2 * d) * sizeof(Type));
  dma_set_stepsize(&dma_get_prev_wave, (ldnz - (snz - 2 * d)) *sizeof(Type));

  dma_set_size(&dma_put_next_wave, (snx - 2 * d) * (snz - 2 * d) * sizeof(Type));
  dma_set_bsize(&dma_put_next_wave, (snz - 2 * d) * sizeof(Type));
  dma_set_stepsize(&dma_put_next_wave, (ldnz - (snz - 2 * d)) *sizeof(Type));

  dma_set_size(&dma_get_vel_one, (snz - 2 * d) *sizeof(Type));
  dma_set_bsize(&dma_get_vel_one, (snz - 2 * d) *sizeof(Type));
  dma_set_stepsize(&dma_get_vel_one, (ldnz - (snz - 2 * d)) * sizeof(Type));

  dma_set_size(&dma_put_next_one, (snz - 2 * d) * sizeof(Type));
  dma_set_bsize(&dma_put_next_one, (snz - 2 * d) * sizeof(Type));
  dma_set_stepsize(&dma_put_next_one, (ldnz - (snz - 2 * d)) *sizeof(Type));

  dma_set_size(&dma_put_curr_one, (snz - 2 * d) * sizeof(Type));
  dma_set_bsize(&dma_put_curr_one, (snz - 2 * d) * sizeof(Type));
  dma_set_stepsize(&dma_put_curr_one, (ldnz - (snz - 2 * d)) *sizeof(Type));

	Type *g_curr_wave = curr_wave + snxbeg * ldnz + snzbeg;
	Type *g_curr_wave_write = curr_wave + (snxbeg + d) * ldnz + (snzbeg + d);
	Type *g_next_wave = next_wave + (snxbeg + d) * ldnz + (snzbeg + d);
	Type *g_prev_wave_read = prev_wave + (snxbeg + d) * ldnz + (snzbeg + d);
	Type *g_vel = vel + (snxbeg + d) * ldnz + (snzbeg + d);

	int ix, iz;
	for(ix = 0 ; ix < 2 * d - 2 + 1 ; ix ++) {
		dma(dma_get_curr_one, (long)(g_curr_wave), (long)(local_curr_wave + ix * snz));
		dma_wait(&replyget_curr_one, 1);
		replyget_curr_one = 0;
		g_curr_wave += ldnz;
	}

	float sponge_coef = 0.015;
	int ib;
  for(ib=0;ib<nb;ib++){
    float tmp=sponge_coef*(nb-ib-1);
    local_bndr[ib]=expf(-tmp*tmp);
  }

  Type *a = (Type *)ldm_malloc(sizeof(Type) * 6);
  a[0] = +1.53400796;
  a[1] = +1.78858721;
  a[2] = -0.31660756;
  a[3] = +0.07612173;
  a[4] = -0.01626042;
  a[5] = +0.00216736;

	int i;
	/*
	SIMDType tmp_u2[1];
	float tmp2_u2[4];
	SIMDType tmp_a[1];
	SIMDType tmp_curr[4];
	float tmp2_curr[16];
	*/
  for (ix = d - 1 ; ix < snx - (d - 1) ; ix++) {
		int ixnext = ix + d < snx ? ix + d : 0;
		dma(dma_get_curr_one, (long)(g_curr_wave), (long)(local_curr_wave + ixnext * snz));
		g_curr_wave = ixnext < snx ? g_curr_wave + ldnz : curr_wave + snxbeg * ldnz + snzbeg;
		//dma_wait(&replyget_curr_one, 1);
		//replyget_curr_one = 0;

		int ixnext2 = ix - (d - 1);
		ixnext2 = ixnext2 < snx - 2 * d ? ixnext2 : snx - 2 * d - 1;
		dma(dma_get_prev_one, (long)(g_prev_wave_read), (long)(local_prev_wave + ixnext2 * (snz - 2 * d)));
		g_prev_wave_read = ixnext2 != snx - 2 * d - 1 ? g_prev_wave_read  + ldnz : g_prev_wave_read;
		//dma_wait(&replyget_prev_one, 1);
		//replyget_prev_one = 0;

		dma(dma_get_vel_one, (long)(g_vel), (long)(local_vel + ixnext2 * (snz - 2 * d)));
		g_vel = ixnext2 != snx - 2 * d - 1 ? g_vel + ldnz : g_vel;

    for (iz = d - 1 ; iz < snz - (d - 1) ; iz++) {
			float tmp;
      int curPos = ix * snz + iz;
			/*
			local_u2[curPos] = -4.0 * a[0] * local_curr_wave[curPos];

			simd_load(tmp_a[0], a + 1);
			for(i = 0 ; i < 4 ; i ++) {
				int j = i + 1;
				tmp2_curr[0 * 4 + i] = local_curr_wave[curPos - j];
				tmp2_curr[1 * 4 + i] = local_curr_wave[curPos + j];
				tmp2_curr[2 * 4 + i] = local_curr_wave[curPos - j * nz];
				tmp2_curr[3 * 4 + i] = local_curr_wave[curPos + j * nz];
			}
			tmp_u2[0] = 0;
			for(i = 0 ; i < 4 ; i ++) {
				simd_load(tmp_curr[i], (tmp2_curr + i * 4));
				tmp_u2[0] += tmp_curr[i];
			}
			//tmp_u2[0] *= tmp_a[0];
			simd_store(tmp_u2[0], (tmp2_u2));
			for(i = 0 ; i < 4 ; i ++) 
				local_u2[curPos] += tmp2_u2[i];
				*/
			/*
      tmp = local_curr_wave[curPos - 1]  +  local_curr_wave[curPos + 1] + local_curr_wave[curPos - snz]  +  local_curr_wave[curPos + snz];
      local_u2[curPos] += a[1] * tmp;
      tmp = local_curr_wave[curPos - 2]  +  local_curr_wave[curPos + 2]  + local_curr_wave[curPos - 2 * snz]  +  local_curr_wave[curPos + 2 * snz];
      local_u2[curPos] += a[2] * tmp;
      tmp = local_curr_wave[curPos - 3]  +  local_curr_wave[curPos + 3]  + local_curr_wave[curPos - 3 * snz]  +  local_curr_wave[curPos + 3 * snz];
      local_u2[curPos] += a[3] * tmp;
      tmp = local_curr_wave[curPos - 4]  +  local_curr_wave[curPos + 4]  + local_curr_wave[curPos - 4 * snz]  +  local_curr_wave[curPos + 4 * snz];
      local_u2[curPos] += a[4] * tmp;
			*/

			/*
      tmp = local_curr_wave[curPos - 5]  +  local_curr_wave[curPos + 5]  + local_curr_wave[curPos - 5 * snz]  +  local_curr_wave[curPos + 5 * snz];
      local_u2[curPos] += a[5] * tmp;
			*/
			/*
      local_u2[curPos] = -4.0 * a[0] * local_curr_wave[curPos];
      tmp = local_curr_wave[curPos - 1]  +  local_curr_wave[curPos + 1] + local_curr_wave[curPos - snz]  +  local_curr_wave[curPos + snz];
      local_u2[curPos] += a[1] * tmp;
      tmp = local_curr_wave[curPos - 2]  +  local_curr_wave[curPos + 2]  + local_curr_wave[curPos - 2 * snz]  +  local_curr_wave[curPos + 2 * snz];
      local_u2[curPos] += a[2] * tmp;
      tmp = local_curr_wave[curPos - 3]  +  local_curr_wave[curPos + 3]  + local_curr_wave[curPos - 3 * snz]  +  local_curr_wave[curPos + 3 * snz];
      local_u2[curPos] += a[3] * tmp;
      tmp = local_curr_wave[curPos - 4]  +  local_curr_wave[curPos + 4]  + local_curr_wave[curPos - 4 * snz]  +  local_curr_wave[curPos + 4 * snz];
      local_u2[curPos] += a[4] * tmp;
      tmp = local_curr_wave[curPos - 5]  +  local_curr_wave[curPos + 5]  + local_curr_wave[curPos - 5 * snz]  +  local_curr_wave[curPos + 5 * snz];
      local_u2[curPos] += a[5] * tmp;
			*/
      local_u2[curPos] = -4.0 * a[0] * local_curr_wave[curPos]+
                   a[1] * (local_curr_wave[curPos - 1]  +  local_curr_wave[curPos + 1]  +
                           local_curr_wave[curPos - snz]  +  local_curr_wave[curPos + snz])  +
                   a[2] * (local_curr_wave[curPos - 2]  +  local_curr_wave[curPos + 2]  +
                           local_curr_wave[curPos - 2 * snz]  +  local_curr_wave[curPos + 2 * snz])  +
                   a[3] * (local_curr_wave[curPos - 3]  +  local_curr_wave[curPos + 3]  +
                           local_curr_wave[curPos - 3 * snz]  +  local_curr_wave[curPos + 3 * snz])  +
                   a[4] * (local_curr_wave[curPos - 4]  +  local_curr_wave[curPos + 4]  +
                           local_curr_wave[curPos - 4 * snz]  +  local_curr_wave[curPos + 4 * snz])  +
                   a[5] * (local_curr_wave[curPos - 5]  +  local_curr_wave[curPos + 5]  +
                           local_curr_wave[curPos - 5 * snz]  +  local_curr_wave[curPos + 5 * snz]);
		}
		dma_wait(&replyget_curr_one, 1);
		replyget_curr_one = 0;
		dma_wait(&replyget_prev_one, 1);
		replyget_prev_one = 0;
		dma_wait(&replyget_vel_one, 1);
		replyget_vel_one = 0;
	}

	dma(dma_put_next_one, (long)g_prev_wave_read, (long)(local_prev_wave)); //just for start, attention
	//dma(dma_put_curr_one, (long)(g_prev_wave_read + ldnz), (long)(local_prev_wave + (snz - 2 * d))); //just for start, attention
  for (ix = d ; ix < snx - d ; ix++) {
    for (iz = d ; iz < snz - d ; iz++) {
      int curPos = ix * snz + iz;
			int curPos2 = (ix - d) * (snz - 2 * d) + (iz - d);
			float curvel = local_vel[curPos2];
      local_prev_wave[curPos2] = 2. * local_curr_wave[curPos] - local_prev_wave[curPos2] +
                          (1.0f / curvel) * local_u2[curPos] + /// 2nd order
                          1.0f / 12 * (1.0f / curvel) * (1.0f / curvel) *
                          (local_u2[curPos - 1] + local_u2[curPos + 1] + local_u2[curPos - snz] + local_u2[curPos + snz] - 4 * local_u2[curPos]); /// 4th order
			/*
			int gx = gxbeg + ix;
			int gxt = ldnx - gx - 1;
			int gz = gzbeg + iz;
			int gzt = ldnz - gz - 1;
			*/
			/*
			if(gx < nb) {
				local_prev_wave[curPos2] *= local_bndr[gx];
				local_curr_wave[curPos] *= local_bndr[gx];
			}
			if(!freeSurface && gz < nb) {
				local_prev_wave[curPos2] *= local_bndr[gz];
				local_curr_wave[curPos] *= local_bndr[gz];
			}
			if(gxt < nb) {
				local_prev_wave[curPos2] *= local_bndr[gxt];
				local_curr_wave[curPos] *= local_bndr[gxt];
			}
			if(gzt < nb) {
				local_prev_wave[curPos2] *= local_bndr[gzt];
				local_curr_wave[curPos] *= local_bndr[gzt];
			}
			*/
    }
		dma_wait(&replyput_next_one, 1);
		replyput_next_one = 0;
		int ix_put_prev = ix != d ?  ix - d : 0;
		dma(dma_put_next_one, (long)g_next_wave, (long)(local_prev_wave + ix_put_prev * (snz - 2 * d)));
		g_next_wave = ix != d ? g_next_wave + ldnz : g_next_wave + ldnz;

		/*
		dma_wait(&replyput_curr_one, 1);
		replyput_curr_one = 0;
		int ix_put_curr = ix != d ?  ix : d;
		dma(dma_put_curr_one, (long)g_curr_wave_write, (long)(local_curr_wave + ix_put_curr * snz + d));
		g_curr_wave_write = ix != d ? g_curr_wave_write + ldnz : g_curr_wave_write + ldnz;
		*/
  }
	dma_wait(&replyput_next_one, 1);
	replyput_next_one = 0;
	//dma_wait(&replyput_curr_one, 1);
	//vreplyput_curr_one = 0;

	ldm_free(local_prev_wave, sizeof(Type) * (snx - 2 * d) * (snz - 2 * d));
	ldm_free(local_curr_wave, sizeof(Type) * snx * snz);
	ldm_free(local_u2, sizeof(Type) * snx * snz);
	ldm_free(local_vel, sizeof(Type) * (snx - 2 * d) * (snz - 2 * d));
	ldm_free(a, sizeof(Type) * 6);
	ldm_free(local_bndr, sizeof(Type) * nb);
}
