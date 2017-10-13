/*
 * ftiframework.cpp
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
#include "ftiframework.h"

FtiFramework::FtiFramework(ForwardModeling &method, const FwiUpdateSteplenOp &updateSteplenOp,
    const FwiUpdateVelOp &_updateVelOp,
    const std::vector<float> &_wlt, const std::vector<float> &_dobs, int _jsx, int _jsz) :
    FwiFramework(method, updateSteplenOp, _updateVelOp, _wlt, _dobs), jsx(_jsx), jsz(_jsz)
{
}

void FtiFramework::epoch(int iter) {
	int nwx = 200;
	std::vector<float> tap = taper(ng, nwx);
	std::vector<float> encobs(ng * nt, 0);
	int rank, np, k, ntask, shot_begin, shot_end;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &np);
	k = std::ceil(ns * 1.0 / np);
	ntask = std::min(k, ns - rank*k);
	shot_begin = rank * k;
	shot_end = shot_begin + ntask;
	float local_obj1 = 0.0f, obj1 = 0.0f;
	int H = 60;
	std::vector<float> img((2 * H + 1) * nx * nz, 0);
	std::vector<float> g2((2 * H + 1) * nx * nz, 0);

	///*
	sf_file sf_g2;
	if(rank == 0 && iter == 0)
	{
		sf_g2 = sf_output("g2.rsf");
		sf_putint(sf_g2, "n1", nz);
		sf_putint(sf_g2, "n2", nx);
		sf_putint(sf_g2, "n3", 2 * H + 1);
	}

	
	for(int is = shot_begin ; is < shot_end ; is ++) {
		std::vector<float> encobs_trans(nt * ng, 0.0f);
		INFO() << format("calculate image, shot id: %d") % is;
		memcpy(&encobs_trans[0], &dobs[is * ng * nt], sizeof(float) * ng * nt);
		for(int it = 0 ; it < nt ; it ++) {
			for(int ig = 0 ; ig < ng ; ig ++) {
				encobs_trans[it * ng + ig] *= tap[ig];
			}
		}
		matrix_transpose(&encobs_trans[0], &encobs[0], ng, nt);

		sf_file shots2 = sf_output("dshots2.rsf");
		sf_putint(shots2, "n1", nt);
		sf_putint(shots2, "n2", ng);
		sf_floatwrite(&encobs[0], nt * ng, shots2);

		//fmMethod.fwiRemoveDirectArrival(&encobs[0], is);
		img.assign((2 * H + 1) * nx * nz, 0.0f);
		image_born(fmMethod, wlt, encobs, img, nt, dt, is, rank, H);
		DEBUG() << ("sum grad: ") << std::accumulate(&img[H * nx * nz], &img[(H + 1) * nx * nz], 0.0f);
		fmMethod.bornMaskGradient(&img[0], H);
		std::transform(g2.begin(), g2.end(), img.begin(), g2.begin(), std::plus<float>());
		DEBUG() << ("global grad: ") << std::accumulate(&g2[H * nx * nz], &g2[(H + 1) * nx * nz], 0.0f);
	}

	img.assign((2 * H + 1) * nx * nz, 0.0f);
	MPI_Allreduce(&g2[0], &img[0], g2.size(), MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&local_obj1, &obj1, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);

	if(rank == 0)
	{
		DEBUG() << ("****** global grad: ") << std::accumulate(&img[H * nx * nz], &img[(H + 1) * nx * nz], 0.0f);
		DEBUG() << format("****** sum obj: %.20f") % obj1;
	}

	if(rank == 0 && iter == 0)
	{
		sf_floatwrite(&img[0], (2 * H + 1) * nx * nz, sf_g2);
	}

	MPI_Barrier(MPI_COMM_WORLD);
	exit(1);
	//*/

	sf_file sf_img0 = sf_input("g3.rsf");
	sf_floatread(&img[0], (2 * H + 1) * nx * nz, sf_img0);

	std::vector<float> gd(nx * nz, 0);
	std::vector<float> grad(nx * nz, 0);
	for(int is = shot_begin ; is < shot_end ; is ++) {
		INFO() << format("************Calculating gradient %d:") % is;
		std::vector<float> encobs_trans(nt * ng, 0.0f);
		memcpy(&encobs_trans[0], &dobs[is * ng * nt], sizeof(float) * ng * nt);
		matrix_transpose(&encobs_trans[0], &encobs[0], ng, nt);	//removeDirectArrival?
		calgradient(fmMethod, wlt, encobs, img, gd, nt, dt, is, rank, H);
	}
	MPI_Allreduce(&gd[0], &grad[0], gd.size(), MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
	fmMethod.maskGradient(&grad[0]);

	if(rank == 0 && iter == 0) {
		 sf_file sf_g2 = sf_output("grad.rsf");
		 sf_putint(sf_g2, "n1", nz);
		 sf_putint(sf_g2, "n2", nx);
		 sf_floatwrite(&grad[0], nx * nz, sf_g2);
	}
	MPI_Barrier(MPI_COMM_WORLD);

	exit(1);

	updateGrad(&g0[0], &img[0], &updateDirection[0], g0.size(), iter);

	/*
		 if(rank == 0 && iter == 0)
		 {
		 sf_file sf_ud = sf_output("updateDirection.rsf");
		 sf_putint(sf_ud, "n1", nz);
		 sf_putint(sf_ud, "n2", nx);
		 sf_floatwrite(&updateDirection[0], nx * nz, sf_ud);
		 }
		 MPI_Barrier(MPI_COMM_WORLD);
		 exit(1);
		 */

	float steplen;
	float obj_val1 = 0, obj_val2 = 0, obj_val3 = 0;

	updateStenlelOp.calsteplen(dobs, updateDirection, obj1, iter, steplen, updateobj, rank, shot_begin, shot_end);


	float alpha1 = updateStenlelOp.alpha1;
	float alpha2 = updateStenlelOp.alpha2;
	float alpha3 = updateStenlelOp.alpha3;
	float obj_val1_sum = updateStenlelOp.obj_val1_sum;
	float obj_val2_sum = updateStenlelOp.obj_val2_sum;
	float obj_val3_sum = updateStenlelOp.obj_val3_sum;
	float maxAlpha3 = updateStenlelOp.maxAlpha3;
	bool	toParabolic = updateStenlelOp.toParabolic;
	updateStenlelOp.parabola_fit(alpha1, alpha2, alpha3, obj_val1_sum, obj_val2_sum, obj_val3_sum, maxAlpha3, toParabolic, iter, steplen, updateobj);

	if(rank == 0)
	{
		INFO() << format("In calculate_steplen(): iter %d  steplen (alpha4) = %e") % iter % steplen;

		INFO() << format("steplen = %.20f") % steplen;
	}

	Velocity &exvel = fmMethod.getVelocity();

	/*
		 if(iter == 1)
		 {
		 sf_file sf_exvel = sf_output("exvel_before.rsf");
		 sf_putint(sf_exvel, "n1", nz);
		 sf_putint(sf_exvel, "n2", nx);
		 sf_floatwrite(&exvel.dat[0], nx * nz, sf_exvel);
		 exit(1);
		 }
		 */

	if(rank == 0)
		INFO() << format("sum vel %f") % sum(exvel.dat);

	updateVelOp.update(exvel, exvel, updateDirection, steplen);

	if(rank == 0)
		INFO() << format("sum vel2 %f") % sum(exvel.dat);

	/*
	if(rank == 0)
	{
		char f_name[64];
		sprintf(f_name, "exvel_after%02d.rsf", iter);
		sf_file sf_exvel2 = sf_output(f_name);
		sf_putint(sf_exvel2, "n1", nz);
		sf_putint(sf_exvel2, "n2", nx);
		sf_floatwrite(&exvel.dat[0], nx * nz, sf_exvel2);
		if(iter == 3)
		exit(1);
	}
	*/

	//fmMethod.refillBoundary(&exvel.dat[0]);

}

void FtiFramework::calgradient(const ForwardModeling &fmMethod,
    const std::vector<float> &wlt,
    const std::vector<float> &vsrc,
		std::vector<float> &img,
    std::vector<float> &gd,
    int nt, float dt,
		int shot_id, int rank, int H)
{
  int nx = fmMethod.getnx();
  int nz = fmMethod.getnz();
  int ns = fmMethod.getns();
  int ng = fmMethod.getng();
	int nb = fmMethod.getbx0();
	std::vector<float> gd0(nx * nz, 0.0f);
  const ShotPosition &allGeoPos = fmMethod.getAllGeoPos();
  const ShotPosition &allSrcPos = fmMethod.getAllSrcPos();

  //std::vector<float> bndr = fmMethod.initBndryVector(nt);
  std::vector<float> sp0(nz * nx, 0);
  std::vector<float> sp1(nz * nx, 0);
  std::vector<float> gp0(nz * nx, 0);
  std::vector<float> gp1(nz * nx, 0);

	std::vector<float> ps(nt * nx * nz, 0);
	std::vector<float> pg(nt * nx * nz, 0);

  ShotPosition curSrcPos = allSrcPos.clipRange(shot_id, shot_id);
	std::vector<float> src = wlt;
	one_order_virtual_source_forth_accuracy(&src[0], nt);
	
	int dn = 10;
	/*
	sf_file fullwv1 = sf_output("fullwv1.rsf");
	sf_putint(fullwv1, "n1" , nz);
	sf_putint(fullwv1, "n2" , nx);
	sf_putint(fullwv1, "n3" , nt / dn);

	sf_file fullwv2 = sf_output("fullwv2.rsf");
	sf_putint(fullwv2, "n1" , nz);
	sf_putint(fullwv2, "n2" , nx);
	sf_putint(fullwv2, "n3" , nt / dn);
	*/

	sf_file fullwv3 = sf_output("fullwv3.rsf");
	sf_putint(fullwv3, "n1" , nz);
	sf_putint(fullwv3, "n2" , nx);
	sf_putint(fullwv3, "n3" , nt / dn);

  for(int it=0; it<nt; it++) {
    //fmMethod.addSource(&sp1[0], &wlt[it], curSrcPos);
    fmMethod.addSource(&sp1[0], &src[it], curSrcPos);
    fmMethod.stepForward(sp0,sp1,0);
    std::swap(sp1, sp0);
		/*
		if(it % dn == 0)
			sf_floatwrite(&sp0[0], nx * nz, fullwv1);
			*/
#pragma omp parallel for
		for(int ix = 0 ; ix < nx ; ix ++)	//lack of the last timestep?
			for(int iz = 0 ; iz < nz ; iz ++)
				//ps[it * nx * nz + ix * nz + iz] = sp1[ix * nz + iz] - sp0[ix * nz + iz];
				ps[it * nx * nz + ix * nz + iz] = sp0[ix * nz + iz];
  }

	printf("1\n");
  std::vector<float> vsrc_trans(ng * nt, 0.0f);
	for(int ig = 0 ; ig < ng ; ig ++)
		one_order_virtual_source_forth_accuracy(const_cast<float*>(&vsrc[ig * nt]), nt);
  matrix_transpose(const_cast<float*>(&vsrc[0]), &vsrc_trans[0], nt, ng);

  for(int it = nt - 1; it >= 0 ; it--) {
		/*
    //fmMethod.readBndry(&bndr[0], &sp0[0], it);	-test
    //std::swap(sp0, sp1); -test
    fmMethod.stepBackward(sp0, sp1);
    //fmMethod.subEncodedSource(&sp0[0], &wlt[it]);
    std::swap(sp0, sp1);	//-test
    fmMethod.subSource(&sp0[0], &wlt[it], curSrcPos);
		*/

    /**
     * forward propagate receviers
     */
    fmMethod.addSource(&gp1[0], &vsrc_trans[it * ng], allGeoPos);
    fmMethod.stepForward(gp0,gp1,0);
    std::swap(gp1, gp0);
		/*
		if(it % dn == 0)
			sf_floatwrite(&gp0[0], nx * nz, fullwv2);
			*/
#pragma omp parallel for
		for(int ix = 0 ; ix < nx ; ix ++)
			for(int iz = 0 ; iz < nz ; iz ++)
				//pg[it * nx * nz + ix * nz + iz] = gp1[ix * nz + iz] - gp0[ix * nz + iz];
				pg[it * nx * nz + ix * nz + iz] = gp0[ix * nz + iz];
	}

	const Velocity &exvel = fmMethod.getVelocity();

	std::vector<float> ggg(nx * nz, 0.0f);
	for(int it=0; it<nt; it++) {
#pragma omp parallel for 
		for(int ix = 0 ; ix < nx ; ix ++) 
			for(int iz = 0 ; iz < nz ; iz ++) 
				ggg[ix * nz + iz] += ps[it * nx * nz + ix * nz + iz] * pg[it * nx * nz + ix * nz + iz];
	}
	/*
	sf_file sf_g2 = sf_output("ggg_test.rsf");
	sf_putint(sf_g2, "n1", nz);
	sf_putint(sf_g2, "n2", nx);
	sf_floatwrite(&ggg[0], nx * nz, sf_g2);
	*/
	printf("2\n");

	sp0.assign(nx * nz, 0);
	sp1.assign(nx * nz, 0);

	/*
	for(int h = -H ; h <= H ; h ++) {
		int ind = h + H;
		for(int ix = nb ; ix < nx - nb ; ix ++) {
			//printf("h = %d, H = %d, ix = %d\n", h, H, ix);
			for(int iz = nb ; iz < nz - nb ; iz ++) { 
				if(h != 0)
					img[ind * nx * nz + ix * nz + iz] = 0.0f;
			}
		}
	}
	*/

	sf_file shots = sf_output("shots_test.rsf");
  sf_putint(shots,"n1",nt);
  sf_putint(shots,"n2",ng);
  std::vector<float> dobs_trans(nt * ng, 0);
  std::vector<float> dobs(nt * ng, 0);
  std::vector<float> record(nx * nz, 0);
	float ps_t = 0, pg_t = 0, img_t = 0;
	for(int it=0; it<nt; it++) {
				for(int h = -H ; h <= H ; h ++) {
					int ind = h + H;
#pragma omp parallel for private(ps_t, pg_t, img_t)
		for(int ix = nb ; ix < nx - nb ; ix ++) {
			//printf("h = %d, H = %d, ix = %d\n", h, H, ix);
			for(int iz = nb ; iz < nz - nb ; iz ++) { 
				//for(int h = -H ; h <= H ; h ++) {
					//int ind = h + H;
					ps_t = ix + 2 * h >= 0 && ix + 2 * h < nx ? ps[it * nx * nz + (ix + 2 * h) * nz + iz] : 0;
					img_t = ix + h >= 0 && ix + h < nx ? img[ind * nx * nz + (ix + h) * nz + iz] : 0; 
					//img_t = ix + h >= nb && ix + h < nx - nb ? img[(ix + h) * nz * (2 * H + 1) + iz * (2 * H + 1) + ind] : 0; 
					sp1[ix * nz + iz] += ps_t * h * h * img_t;
					//sp1[ix * nz + iz] += ps_t * img_t;
				}
			}
		}
		fmMethod.stepForward(sp0,sp1,0);
		std::swap(sp1, sp0);
		if(it % dn == 0)
			sf_floatwrite(&sp0[0], nx * nz, fullwv3);
    fmMethod.recordSeis(&dobs_trans[it*ng], &sp0[0]);
#pragma omp parallel for 
		for(int ix = 0 ; ix < nx ; ix ++) 
			for(int iz = 0 ; iz < nz ; iz ++) 
				gd0[ix * nz + iz] += 2 * sp0[ix * nz + iz] * pg[it * nx * nz + ix * nz + iz] * exvel.dat[ix * nz + iz];
	}
  matrix_transpose(&dobs_trans[0], &dobs[0], ng, nt);
	for(int ig = 0 ; ig < ng ; ig ++)
		for(int it = 0 ; it < nt ; it ++)
			if(it == 0 || it == 1)
				dobs[ig * nt + it] = 0.0f;
			else if(it == nt - 1 || it == nt - 2)
				dobs[ig * nt + it] = 0.0f;
			else
				dobs[ig * nt + it] = (-dobs[ig * nt + it - 1] + dobs[ig * nt + it + 1]) / 2;
	sf_floatwrite(&dobs[0], ng*nt, shots);

	/*
	sf_file gd1 = sf_output("gd1_test.rsf");
  sf_putint(gd1,"n1",nz);
  sf_putint(gd1,"n2",nx);
	sf_floatwrite(&gd0[0], nx*nz, gd1);
	*/

	printf("3\n");
	gp0.assign(nx * nz, 0);
	gp1.assign(nx * nz, 0);

	for(int it = nt - 1; it >= 0 ; it--) {
				for(int h = -H ; h <= H ; h ++) {
					int ind = h + H;
#pragma omp parallel for private(ps_t, pg_t, img_t)
		for(int ix = nb ; ix < nx - nb ; ix ++) 
			for(int iz = nb ; iz < nz - nb ; iz ++) {
				//for(int h = -H ; h <= H ; h ++) {
					//int ind = h + H;
					pg_t = ix - 2 * h >= 0 && ix - 2 * h < nx ? pg[it * nx * nz + (ix - 2 * h) * nz + iz] : 0;
					img_t = ix - h >= 0 && ix - h < nx ? img[ind * nx * nz + (ix - h) * nz + iz] : 0;
					//img_t = ix - h >= nb && ix - h < nx - nb ? img[(ix - h) * nz * (2 * H + 1) + iz * (2 * H + 1) + ind] : 0;
					gp1[ix * nz + iz] +=  pg_t * h * h * img_t;
					//gp1[ix * nz + iz] +=  pg_t * img_t;
			}
		}
    fmMethod.stepForward(gp0,gp1,0);
    std::swap(gp1, gp0);
#pragma omp parallel for 
		for(int ix = 0 ; ix < nx ; ix ++) 
			for(int iz = 0 ; iz < nz ; iz ++) 
				gd0[ix * nz + iz] += 2 * gp0[ix * nz + iz] * ps[it * nx * nz + ix * nz + iz] * exvel.dat[ix * nz + iz];
	}
	printf("4\n");
	char filename[20];
	sprintf(filename, "grad_%d.rsf", shot_id);
	sf_file sf_g2 = sf_output(filename);
	sf_putint(sf_g2, "n1", nz);
	sf_putint(sf_g2, "n2", nx);
	sf_floatwrite(&gd[0], nx * nz, sf_g2);

	std::transform(gd.begin(), gd.end(), gd0.begin(), gd.begin(), std::plus<float>());
}

void FtiFramework::image_born(const ForwardModeling &fmMethod,
    const std::vector<float> &wlt,
    const std::vector<float> &vsrc,
    std::vector<float> &g0,
    int nt, float dt,
		int shot_id, int rank, int H)
{
  int nx = fmMethod.getnx();
  int nz = fmMethod.getnz();
  int ns = fmMethod.getns();
  int ng = fmMethod.getng();
  const ShotPosition &allGeoPos = fmMethod.getAllGeoPos();
  const ShotPosition &allSrcPos = fmMethod.getAllSrcPos();

  //std::vector<float> bndr = fmMethod.initBndryVector(nt);
  std::vector<float> sp0(nz * nx, 0);
  std::vector<float> sp1(nz * nx, 0);
  std::vector<float> gp0(nz * nx, 0);
  std::vector<float> gp1(nz * nx, 0);


  ShotPosition curSrcPos = allSrcPos.clipRange(shot_id, shot_id);

	std::vector<float> ps(nt * nx * nz, 0);

  for(int it=0; it<nt; it++) {
    fmMethod.addSource(&sp1[0], &wlt[it], curSrcPos);
    //fmMethod.stepForward(sp0,sp1);
    fmMethod.stepForward(sp0,sp1,0);
    std::swap(sp1, sp0);
#pragma omp parallel for
		for(int ix = 0 ; ix < nx ; ix ++)	//lack of the last timestep?
			for(int iz = 0 ; iz < nz ; iz ++)
				ps[it * nx * nz + ix * nz + iz] = sp0[ix * nz + iz];
    //fmMethod.writeBndry(&bndr[0], &sp0[0], it); //-test
		/*
		const int check_step = 5;
    if ((it > 0) && (it != (nt - 1)) && !(it % check_step)) {
      char check_file_name1[64];
      char check_file_name2[64];
      sprintf(check_file_name1, "./rank_%d_check_time_%d_1.su", rank, it);
      sprintf(check_file_name2, "./rank_%d_check_time_%d_2.su", rank, it);
			FILE *f1 = fopen(check_file_name1, "wb");
			FILE *f2 = fopen(check_file_name2, "wb");
			fwrite(&sp0[0], sizeof(float), nx * nz, f1);
			fwrite(&sp1[0], sizeof(float), nx * nz, f2);
			fclose(f1);
			fclose(f2);
    }
		*/
  }
	/*
	char check_file_name1[64];
	char check_file_name2[64];
	sprintf(check_file_name1, "./rank_%d_check_time_last_1.su", rank);
	sprintf(check_file_name2, "./rank_%d_check_time_last_2.su", rank);
	FILE *f1 = fopen(check_file_name1, "wb");
	FILE *f2 = fopen(check_file_name2, "wb");
	fwrite(&sp0[0], sizeof(float), nx * nz, f1);
	fwrite(&sp1[0], sizeof(float), nx * nz, f2);
	fclose(f1);
	fclose(f2);
	*/

  std::vector<float> vsrc_trans(ng * nt, 0.0f);
  matrix_transpose(const_cast<float*>(&vsrc[0]), &vsrc_trans[0], nt, ng);

  for(int it = nt - 1; it >= 0 ; it--) {
		/*
		const int check_step = 5;
		if(it == nt - 1)
		{
			char check_file_name1[64];
			char check_file_name2[64];
			sprintf(check_file_name1, "./rank_%d_check_time_last_1.su", rank);
			sprintf(check_file_name2, "./rank_%d_check_time_last_2.su", rank);
			FILE *f1 = fopen(check_file_name1, "rb");
			FILE *f2 = fopen(check_file_name2, "rb");
			fread(&sp1[0], sizeof(float), nx * nz, f1);
			fread(&sp0[0], sizeof(float), nx * nz, f2);
			fclose(f1);
			fclose(f2);
		}
		else if ((check_step > 0) && !(it % check_step) && (it != 0)) {
			char check_file_name1[64];
			char check_file_name2[64];
			sprintf(check_file_name1, "./rank_%d_check_time_%d_1.su", rank, it);
			sprintf(check_file_name2, "./rank_%d_check_time_%d_2.su", rank, it);
			FILE *f1 = fopen(check_file_name1, "rb");
			FILE *f2 = fopen(check_file_name2, "rb");
			fread(&sp1[0], sizeof(float), nx * nz, f1);
			fread(&sp0[0], sizeof(float), nx * nz, f2);
			fclose(f1);
			fclose(f2);
		}
    //fmMethod.readBndry(&bndr[0], &sp0[0], it);	//-test
    //std::swap(sp0, sp1); //-test
    fmMethod.stepBackward(sp0, sp1);
    //fmMethod.subEncodedSource(&sp0[0], &wlt[it]);
    std::swap(sp0, sp1); //-test
    fmMethod.subSource(&sp0[0], &wlt[it], curSrcPos);
		*/

    /**
     * forward propagate receviers
     */
    fmMethod.addSource(&gp1[0], &vsrc_trans[it * ng], allGeoPos);
    fmMethod.stepForward(gp0,gp1,0);
    std::swap(gp1, gp0);

    cross_correlation(&ps[it * nx * nz], &gp0[0], &g0[0], nx, nz, 1.0, H);
 }
}

void FtiFramework::cross_correlation(float *src_wave, float *vsrc_wave, float *image, int nx, int nz, float scale, int H) {
	float t_src_wave, t_vsrc_wave;
	int nb = fmMethod.getbx0();
#pragma omp parallel for private(t_src_wave, t_vsrc_wave)
	for(int h = -H ; h <= H ; h ++) {
		int ind = h + H;
		/*
		for (int i = 0 ; i < nx ; i ++) {
			for (int j = 0 ; j < nz ; j ++) {
				t_src_wave = i + h < nx && i + h >= 0 ? src_wave[(i + h) * nz + j] : 0;
				t_vsrc_wave = i - h < nx && i - h >= 0 ? vsrc_wave[(i - h) * nz + j] : 0;
				image[ind * nx * nz + i * nz + j] += t_src_wave * t_vsrc_wave * scale;
			}
		}
		*/

		for (int i = nb ; i < nx - nb ; i ++) {
			for (int j = nb ; j < nz - nb ; j ++) {
				t_src_wave = i + h < nx - nb && i + h >= nb ? src_wave[(i + h) * nz + j] : 0;
				t_vsrc_wave = i - h < nx - nb && i - h >= nb ? vsrc_wave[(i - h) * nz + j] : 0;
				image[ind * nx * nz + i * nz + j] += t_src_wave * t_vsrc_wave * scale;
			}
		}
	}
}

