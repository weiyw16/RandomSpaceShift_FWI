/*
 * fd4t10s-damp-zjh.h
 *
 *  Created on: Mar 4, 2016
 *      Author: cbw
 */

#ifndef SRC_MDLIB_FD4T10S_DAMP_ZJH_CG_H_
#define SRC_MDLIB_FD4T10S_DAMP_ZJH_CG_H_
#include <string.h>
#include <pthread.h>

//pthread
#define CG_MAX 4
#define X_MAX 4
#define Z_MAX 1
#define ND (X_MAX*Z_MAX/CG_MAX)

//athread
#define ROW 8 //change
#define COL 8 

const int d = 6;

int g_status;
pthread_mutex_t mut1 = PTHREAD_MUTEX_INITIALIZER, mut2 = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t cond1 = PTHREAD_COND_INITIALIZER, cond2 = PTHREAD_COND_INITIALIZER;

typedef float Type;
struct Param {
		int task;

		//forward modeling
		Type *prev_wave;
		Type *curr_wave;
		Type *next_wave;
		Type *vel;
		Type *u2;
		int nx;
		int nz;
		int ldnx;
		int ldnz;
		int nb;
		int freeSurface;
		int id;
		int nd;
		int icg;
		int status;
		int exit;
		int gnxbeg;
		int gnzbeg;
		struct Param **params;

		//cross correlation
		int crnx;
		int crnz;
		float scale;
		float *src_wave;
		float *vsrc_wave;
		float *image;
};

void cross_correlation_cg(float *src_wave, float *vsrc_wave, float *image, int nx, int nz, float scale);
void fd4t10s_nobndry_zjh_2d_vtrans_cg(const float *prev_wave, const float *curr_wave, float *next_wave, const float *vel, float *u2, int nx, int nz, int nb, int nt, int freeSurface);
void fd4t10s_4cg_exit();

#endif /* SRC_MDLIB_FD4T10S_DAMP_ZJH_H_ */
