/*
 * fd4t10s-damp-zjh.cpp
 *
 *  Created on: Mar 4, 2016
 *      Author: cbw
 */

#include <stdio.h>
#include <athread.h>
#include "fd4t10s-damp-zjh-cg.h"
#include "fd4t10s-slave.h"

extern SLAVE_FUN(fd4t10s_sw)();
extern SLAVE_FUN(cross_correlation_sw)();

struct Param;
/**
 * please note that the velocity is transformed
 */
void master_wait() {
	pthread_mutex_lock(&mut1);
	while(g_status != 0)
	{
		pthread_cond_wait(&cond1, &mut1);
	}
	pthread_mutex_unlock(&mut1);
	g_status = CG_MAX;
}

void master_signal(struct Param** params) {
	int i;
	pthread_mutex_lock(&mut2);
	for(i = 0 ; i < CG_MAX ; i ++)
		params[i * ND]->status = 0;
	pthread_cond_broadcast(&cond2);
	pthread_mutex_unlock(&mut2);
}

void slave_wait(struct Param* param) {
	pthread_mutex_lock(&mut2);
	while(param->status)
		pthread_cond_wait(&cond2, &mut2);
	pthread_mutex_unlock(&mut2);
}

void slave_signal() {
	pthread_mutex_lock(&mut1);
	g_status --;
	pthread_cond_signal(&cond1);
	pthread_mutex_unlock(&mut1);
}

void compute_ldm(struct Param *param) {
	int nx = param->nx;
	int nz = param->nz;
	int nb = param->nb;
	int cnx = nx - 2 * d;
	int cnz = nz - 2 * d;
	int snx = ceil(cnx * 1.0 / ROW);
	int snz = ceil(cnz * 1.0 / COL);
	int cid = 0;
	int rid = 0;
	int snxbeg = snx * cid - d + d; 
	int snzbeg = snz * rid - d + d;
	int snxend = snx * (cid + 1) + d + d;
	int snzend = snz * (rid + 1) + d + d;
	snxend = snxend < nx ? snxend : nx;
	snzend = snzend < nz ? snzend : nz;
	snx = snxend - snxbeg;
	snz = snzend - snzbeg;
	int vel_size = (snx - 2 * d) * (snz - 2 * d) * sizeof(float);
	int curr_size = snx * snz * sizeof(float);
	int total = vel_size * 2 + curr_size * 2 + nb * sizeof(float) + 6 * sizeof(float);
	printf("LDM consume: vel_size = %dB, curr_size = %dB, total = %dB\n", vel_size, curr_size, total);
}

int fd4t10s_4cg(void *ptr) {
	athread_init();
	struct Param *param = (struct Param *)ptr;
	struct Param **params = param->params;
	int id = param->id;
	compute_ldm(param);
	while(1) {
		slave_wait(param);
		if(param->exit)
			break;
		param->status = 1;
		int i, j;
		for(i = 0 ; i < ND ; i ++) {
			param = params[id + i];
			if(param->task == 0) {
				athread_spawn(fd4t10s_sw, param);
				athread_join();
			}
			else if(param->task == 1) {
				athread_spawn(cross_correlation_sw, param);
				athread_join();
				/*
				int ix, iz;
				for(ix = 0 ; ix < param->crnx ; ix ++) {
					for(iz = 0 ; iz < param->crnz ; iz ++) {
						j = ix * param->ldnz + iz;
						param->image[j] -= param->src_wave[j] * param->vsrc_wave[j] * param->scale;
					}
				}
				*/
			}
			param = params[id];
		}
		slave_signal();
	}
}

int id, icg, ix, iz;
int cnx, cnz;
int snx, snz;
int snxbeg, snzbeg;
int snxend, snzend;

pthread_t pt[CG_MAX];
struct Param *params[X_MAX * Z_MAX];
void fd4t10s_nobndry_zjh_2d_vtrans_cg(const float *prev_wave, const float *curr_wave, float *next_wave, const float *vel, float *u2, int nx, int nz, int nb, int nt, int freeSurface) {
	static int init = 1;
	if(init) {
		if(pthread_mutex_init(&mut1, NULL) != 0)
			printf("mutex init error\n");

		if(pthread_cond_init(&cond1, NULL) != 0)
			printf("cond init error\n");

		for(ix = 0 ; ix < X_MAX ; ix ++) {
			cnx = nx - 2 * d;
			snx = ceil(cnx * 1.0 / X_MAX);
			snxbeg = snx * ix - d + d;
			snxend = snx * (ix + 1) + d + d;
			snxend = snxend < nx ? snxend : nx;
			snx = snxend - snxbeg;
			for(iz = 0 ; iz < Z_MAX ; iz ++) {
				cnz = nz - 2 * d;
				snz = ceil(cnz * 1.0 / Z_MAX);
				snzbeg = snz * iz - d + d;
				snzend = snz * (iz + 1) + d + d;
				snzend = snzend < nz ? snzend : nz;
				snz = snzend - snzbeg;
				id = ix * Z_MAX + iz;
				params[id] = (struct Param *)malloc(sizeof(struct Param));
				params[id]->nx = snx;
				params[id]->nz = snz;
				params[id]->ldnx = nx;
				params[id]->ldnz = nz;
				params[id]->nb = nb;
				params[id]->freeSurface = freeSurface;
				params[id]->id = id;
				params[id]->nd = ND;
				params[id]->icg = id / ND;
				params[id]->status = 1;	//stop
				params[id]->exit = 0;		//can not exit
				params[id]->gnxbeg = snxbeg;	
				params[id]->gnzbeg = snzbeg;
				params[id]->params = params;
			}
		}
		g_status = CG_MAX;
	}

	if(init) {
		for(icg = 0 ; icg < CG_MAX ; icg ++) {
			int id_beg = icg * ND;
			pthread_create(&pt[icg], NULL, fd4t10s_4cg, (void *)params[id_beg]);
		}
		init = 0;
	}
	
	for(ix = 0 ; ix < X_MAX ; ix ++) {
		cnx = nx - 2 * d;
		snx = ceil(cnx * 1.0 / X_MAX);
		snxbeg = snx * ix - d + d;
		for(iz = 0 ; iz < Z_MAX ; iz ++) {
			cnz = nz - 2 * d;
			snz = ceil(cnz * 1.0 / Z_MAX);
			snzbeg = snz * iz - d + d;
			id = ix * Z_MAX + iz;
			params[id]->prev_wave = prev_wave + snxbeg * nz + snzbeg;
			params[id]->curr_wave = curr_wave + snxbeg * nz + snzbeg;
			params[id]->next_wave = next_wave + snxbeg * nz + snzbeg;
			params[id]->vel = vel + snxbeg * nz + snzbeg;
			params[id]->u2 = u2 + snxbeg * nz + snzbeg;
			params[id]->task = 0;
		}
	}

	master_signal(params);
	master_wait();
}

void cross_correlation_cg(float *src_wave, float *vsrc_wave, float *image, int nx, int nz, float scale) {
	for(ix = 0 ; ix < X_MAX ; ix ++) {
		snx = ceil(nx * 1.0 / X_MAX);
		snxbeg = snx * ix;
		snxend = snx * (ix + 1);
		snxend = snxend < nx ? snxend : nx;
		snx = snxend - snxbeg;
		for(iz = 0 ; iz < Z_MAX ; iz ++) {
			snz = ceil(nz * 1.0 / Z_MAX);
			snzbeg = snz * iz;
			snzend = snz * (iz + 1);
			snzend = snzend < nz ? snzend : nz;
			snz = snzend - snzbeg;
			id = ix * Z_MAX + iz;
			params[id]->src_wave = src_wave + snxbeg * nz + snzbeg;
			params[id]->vsrc_wave = vsrc_wave + snxbeg * nz + snzbeg;
			params[id]->image = image + snxbeg * nz + snzbeg;
			params[id]->scale = scale;
			params[id]->crnx = snx;
			params[id]->crnz = snz;
			params[id]->task = 1;
		}
	}

	master_signal(params);
	master_wait();
}

void fd4t10s_4cg_exit() {
	int icg;
	for(icg = 0 ; icg < CG_MAX ; icg ++)
		params[icg * ND]->exit = 1;
	master_signal(params);

	for(icg = 0 ; icg < CG_MAX ; icg ++) {
		pthread_join(pt[icg], NULL);
	}

	int ix, iz;
	for(ix = 0 ; ix < X_MAX ; ix ++)
		for(iz = 0 ; iz < Z_MAX ; iz ++)
			free(params[ix * Z_MAX + iz]);
}
