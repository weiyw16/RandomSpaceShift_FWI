
extern "C" {
#include <rsf.h>
}

#ifdef _OPENMP
#include <omp.h>
#endif

#include "mpi.h"
#include "logger.h"
#include "common.h"
#include "sum.h"
#include "ricker-wavelet.h"
#include "velocity.h"
#include "sf-velocity-reader.h"
#include "common.h"
#include "shot-position.h"
#include "forwardmodeling.h"
#include "sfutil.h"
#include "timer.h"
#include "environment.h"
#include "math.h"

#define OFFSET 1000000
namespace {
class Params {
public:
  Params();
  ~Params();

private:
  Params(const Params &);
  void operator=(const Params &);
  void check();

public:
  sf_file vinit;
  sf_file vreal;
  sf_file shots_rf;
  sf_file shots_bg;
  int nb;
  int nz;
  int nx;
  float dz;
  float dx;
  int nt;
  int ng;
  int ns;
  float dt;
  float amp;
  float fm;
  int sxbeg;
  int szbeg;
  int gxbeg;
  int gzbeg;
  int jsx;
  int jsz;
  int jgx;
  int jgz;
	int freeSurface;

public:
  int rank;
  int k;
  int np;
  int ntask; /// exactly the # of task each process owns
};


Params::Params() {
  /*< set up I/O files >*/
  vinit=sf_input ("vinit");   /* initial velocity model, unit=m/s */
  vreal=sf_input ("vreal");   /* initial velocity model, unit=m/s */
  shots_rf=sf_output("shots_rf");
  shots_bg=sf_output("shots_bg");

  /* get parameters for forward modeling */
  if (!sf_histint(vinit,"n1",&nz)) sf_error("no n1");
  if (!sf_histint(vinit,"n2",&nx)) sf_error("no n2");
  if (!sf_histfloat(vinit,"d1",&dz)) sf_error("no d1");
  if (!sf_histfloat(vinit,"d2",&dx)) sf_error("no d2");

  if (!sf_getfloat("amp",&amp)) amp=1000;
  /* maximum amplitude of ricker */
  if (!sf_getfloat("fm",&fm)) fm=10;
  /* dominant freq of ricker */
  if (!sf_getint("nb",&nb))   nb=30;
  /* thickness of sponge ABC  */
  if (!sf_getfloat("dt",&dt)) sf_error("no dt");
  /* time interval */
  if (!sf_getint("nt",&nt))   sf_error("no nt");
  /* total modeling time steps */
  if (!sf_getint("ns",&ns))   sf_error("no ns");
  /* total shots */
  if (!sf_getint("ng",&ng))   sf_error("no ng");
  /* total receivers in each shot */
  if (!sf_getint("jsx",&jsx))   sf_error("no jsx");
  /* source x-axis  jump interval  */
  if (!sf_getint("jsz",&jsz))   jsz=0;
  /* source z-axis jump interval  */
  if (!sf_getint("jgx",&jgx))   jgx=1;
  /* receiver x-axis jump interval */
  if (!sf_getint("jgz",&jgz))   jgz=0;
  /* receiver z-axis jump interval */
  if (!sf_getint("sxbeg",&sxbeg))   sf_error("no sxbeg");
  /* x-begining index of sources, starting from 0 */
  if (!sf_getint("szbeg",&szbeg))   sf_error("no szbeg");
  /* z-begining index of sources, starting from 0 */
  if (!sf_getint("gxbeg",&gxbeg))   sf_error("no gxbeg");
  /* x-begining index of receivers, starting from 0 */
  if (!sf_getint("gzbeg",&gzbeg))   sf_error("no gzbeg");
  /* z-begining index of receivers, starting from 0 */
	if (!sf_getint("free", &freeSurface)) sf_error("no freeSurface");
	/* whether it is freeSurface */

  sf_putint(shots_rf,"n1",nt);
  sf_putint(shots_rf,"n2",ng);
  sf_putint(shots_rf,"n3",ns);
  sf_putfloat(shots_rf,"d1",dt);
  sf_putfloat(shots_rf,"d2",jgx*dx);
  sf_putfloat(shots_rf,"o1",0);
  sf_putfloat(shots_rf,"o2",0);
  sf_putstring(shots_rf,"label1","Time");
  sf_putstring(shots_rf,"label2","Lateral");
  sf_putstring(shots_rf,"label3","Shot");
  sf_putstring(shots_rf,"unit1","sec");
  sf_putstring(shots_rf,"unit2","m");
  sf_putfloat(shots_rf,"amp",amp);
  sf_putfloat(shots_rf,"fm",fm);
  sf_putint(shots_rf,"ng",ng);
  sf_putint(shots_rf,"szbeg",szbeg);
  sf_putint(shots_rf,"sxbeg",sxbeg);
  sf_putint(shots_rf,"gzbeg",gzbeg);
  sf_putint(shots_rf,"gxbeg",gxbeg);
  sf_putint(shots_rf,"jsx",jsx);
  sf_putint(shots_rf,"jsz",jsz);
  sf_putint(shots_rf,"jgx",jgx);
  sf_putint(shots_rf,"jgz",jgz);
  sf_putint(shots_rf, "nb", nb);
  sf_putint(shots_rf, "free", freeSurface);

  sf_putint(shots_bg,"n1",nt);
  sf_putint(shots_bg,"n2",ng);
  sf_putint(shots_bg,"n3",ns);
  sf_putfloat(shots_bg,"d1",dt);
  sf_putfloat(shots_bg,"d2",jgx*dx);
  sf_putfloat(shots_bg,"o1",0);
  sf_putfloat(shots_bg,"o2",0);
  sf_putstring(shots_bg,"label1","Time");
  sf_putstring(shots_bg,"label2","Lateral");
  sf_putstring(shots_bg,"label3","Shot");
  sf_putstring(shots_bg,"unit1","sec");
  sf_putstring(shots_bg,"unit2","m");
  sf_putfloat(shots_bg,"amp",amp);
  sf_putfloat(shots_bg,"fm",fm);
  sf_putint(shots_bg,"ng",ng);
  sf_putint(shots_bg,"szbeg",szbeg);
  sf_putint(shots_bg,"sxbeg",sxbeg);
  sf_putint(shots_bg,"gzbeg",gzbeg);
  sf_putint(shots_bg,"gxbeg",gxbeg);
  sf_putint(shots_bg,"jsx",jsx);
  sf_putint(shots_bg,"jsz",jsz);
  sf_putint(shots_bg,"jgx",jgx);
  sf_putint(shots_bg,"jgz",jgz);
  sf_putint(shots_bg, "nb", nb);
  sf_putint(shots_bg, "free", freeSurface);

	//!!!!!You should put the vmin and vmax of vreal not vinit to the shots, because fti will use it as input!!!!!!
  Velocity v = SfVelocityReader::read(vreal, nx, nz);
  float vmin = *std::min_element(v.dat.begin(), v.dat.end());
  float vmax = *std::max_element(v.dat.begin(), v.dat.end());
  sf_putfloat(shots_rf, "vmin", vmin);
  sf_putfloat(shots_rf, "vmax", vmax);

  sf_putfloat(shots_bg, "vmin", vmin);
  sf_putfloat(shots_bg, "vmax", vmax);

  MPI_Comm_size(MPI_COMM_WORLD, &np);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  k = std::ceil(ns * 1.0 / np);
  ntask = std::min(k, ns - rank*k);

  check();
}

Params::~Params() {
  sf_close();
}

void Params::check() {
  if (!(sxbeg >= 0 && szbeg >= 0 && sxbeg + (ns - 1)*jsx < nx && szbeg + (ns - 1)*jsz < nz)) {
    sf_warning("sources exceeds the computing zone!\n");
    exit(1);
  }

  if (!(gxbeg >= 0 && gzbeg >= 0 && gxbeg + (ng - 1)*jgx < nx && gzbeg + (ng - 1)*jgz < nz)) {
    sf_warning("geophones exceeds the computing zone!\n");
    exit(1);
  }

}

} /// end of name space

int main(int argc, char* argv[]) {
  MPI_Init(&argc, &argv);
	MPI_Status status;
	MPI_Request request;
  /* initialize Madagascar */
  sf_init(argc,argv);
  Environment::setDatapath();

  Params params;
  Timer totalTimer;

  /// configure logger
	char logfile[64];
	sprintf(logfile, "born-damp-%02d.log", params.rank);
  FILELog::setLogFile(logfile);
	printGitInfo();

  int nz = params.nz;
  int nx = params.nx;
  int nb = params.nb;
  int ng = params.ng;
  int nt = params.nt;
  int ns = params.ns;
  float dt = params.dt;
  float fm = params.fm;
  int np = params.np;
  int rank = params.rank;
  int k = params.k;
  int ntask = params.ntask;

  ShotPosition allSrcPos(params.szbeg, params.sxbeg, params.jsz, params.jsx, ns, nz);
  ShotPosition allGeoPos(params.gzbeg, params.gxbeg, params.jgz, params.jgx, ng, nz);
  ForwardModeling fmMethod(allSrcPos, allGeoPos, dt, params.dx, params.fm, nb, nt, params.freeSurface);

  Velocity exvel = fmMethod.expandDomain(SfVelocityReader::read(params.vinit, nx, nz));
  Velocity exvel_real = fmMethod.expandDomain(SfVelocityReader::read(params.vreal, nx, nz));

	nx = exvel.nx;
	nz = exvel.nz;

#ifdef dot_product_test
  fmMethod.bindVelocity(exvel);

  std::vector<float> wlt(nt);
  rickerWavelet(&wlt[0], nt, fm, dt, params.amp);

  std::vector<float> dobs(params.ntask * params.nt * params.ng, 0);
  std::vector<float> fullwv(nt * nz * nx, 0);

	std::vector<float> p0(nz * nx, 0);
	std::vector<float> p1(nz * nx, 0);
	std::vector<float> dobs_trans(params.nt * params.ng, 0);
	ShotPosition curSrcPos = allSrcPos.clipRange(is, is);

	std::vector<float> rand1(params.nt * params.ns, 0);
	std::vector<float> rand2(params.nt * params.ng, 0);
	
	for(int it=0; it<nt; it++) {
		fmMethod.addSource(&p1[0], &rand1[it * ns], allSrcPos);
		fmMethod.stepForward(p0,p1);
		std::swap(p1, p0);
		fmMethod.recordSeis(&dobs_trans[it*ng], &p0[0]);
	}

	p0.assign(nz * nx, 0);
	p1.assign(nz * nx, 0);
	dobs_trans.assign(params.nt * params.ng, 0);
	for(int it=0; it<nt; it++) {
		fmMethod.addSource(&p1[0], &rand2[it * ng], allGeoPos);
		fmMethod.stepForward(p0,p1);
		std::swap(p1, p0);
		fmMethod.recordSeis(&dobs_trans[it*ng], &p0[0]);
	}
	
#endif 

#define gradient_test
#ifdef gradient_test
	std::vector<float> exvel_m(nx * nz, 0);
	vectorMinus(exvel_real.dat, exvel.dat, exvel_m);

  fmMethod.bindVelocity(exvel);

  std::vector<float> wlt(nt);
  rickerWavelet(&wlt[0], nt, fm, dt, params.amp);

  std::vector<float> dobs(params.ntask * params.nt * params.ng, 0);
  std::vector<float> dobs_t(params.ntask * params.nt * params.ng, 0);
  std::vector<float> fullwv(3 * nz * nx, 0);
	// fullwv_t0: t - 1 timestep; fulllwv_t1: t timestep; fullwv_t2: t + 1 timestep
	float *fullwv_t0, *fullwv_t1, *fullwv_t2;	
	fullwv_t0 = &fullwv[0];
	fullwv_t1 = &fullwv[nz * nx];
	fullwv_t2 = &fullwv[2 * nz * nx];
  for(int is=rank*k; is<rank*k+ntask; is++) {
    int local_is = is - rank * k;
    Timer timer;
    std::vector<float> dobs_trans(params.nt * params.ng, 0);
    std::vector<float> dobs_trans_t(params.nt * params.ng, 0);
		//fmMethod.BornForwardModeling(exvel_m, wlt, dobs_trans, is);
    std::vector<float> p0(nz * nx, 0);
    std::vector<float> p1(nz * nx, 0);
    std::vector<float> rp0(nz * nx, 0);
    std::vector<float> rp1(nz * nx, 0);
    ShotPosition curSrcPos = allSrcPos.clipRange(is, is);

    for(int it0 = 0 ; it0 < nt + 1 ; it0 ++) {
      fmMethod.addSource(&p1[0], &wlt[it0], curSrcPos);
      fmMethod.stepForward(p0,p1,0);
      std::swap(p1, p0);
			if(it0 < nt)
				fmMethod.recordSeis(&dobs_trans_t[it0*ng], &p0[0]);

			swap3(fullwv_t0, fullwv_t1, fullwv_t2);
			std::copy(p0.begin(), p0.end(), fullwv_t2);

			int it = it0 - 1;
			if(it < 0)
				continue;
			fmMethod.addBornwv(fullwv_t0, fullwv_t1, fullwv_t2, &exvel_m[0], dt, it, &rp1[0]);
      fmMethod.stepForward(rp0,rp1,1);
      std::swap(rp1, rp0);
      fmMethod.recordSeis(&dobs_trans[it*ng], &rp0[0]);
    }

    matrix_transpose(&dobs_trans[0], &dobs[local_is * ng * nt], ng, nt);
		if(np == 1) {
			sf_floatwrite(&dobs[local_is * ng * nt], ng*nt, params.shots_rf);
		}
		else {
			if(rank == 0) {
				sf_floatwrite(&dobs[local_is * ng * nt], ng*nt, params.shots_rf);
				if(is == rank * k + ntask - 1) {
					for(int other_is = rank * k + ntask ; other_is < ns ; other_is ++) {
						MPI_Recv(&dobs[0], ng*nt, MPI_FLOAT, other_is / k, other_is, MPI_COMM_WORLD, &status);
						sf_floatwrite(&dobs[0], ng*nt, params.shots_rf);
					}
				}
			}
			else {
				MPI_Isend(&dobs[local_is * ng * nt], ng*nt, MPI_FLOAT, 0, is, MPI_COMM_WORLD, &request);
			}
		}

    matrix_transpose(&dobs_trans_t[0], &dobs_t[local_is * ng * nt], ng, nt);
		if(np == 1) {
			sf_floatwrite(&dobs_t[local_is * ng * nt], ng*nt, params.shots_bg);
		}
		else {
			if(rank == 0) {
				sf_floatwrite(&dobs_t[local_is * ng * nt], ng*nt, params.shots_bg);
				if(is == rank * k + ntask - 1) {
					for(int other_is = rank * k + ntask ; other_is < ns ; other_is ++) {
						MPI_Recv(&dobs_t[0], ng*nt, MPI_FLOAT, other_is / k, other_is + OFFSET, MPI_COMM_WORLD, &status);
						sf_floatwrite(&dobs_t[0], ng*nt, params.shots_bg);
					}
				}
			}
			else {
				MPI_Isend(&dobs_t[local_is * ng * nt], ng*nt, MPI_FLOAT, 0, is + OFFSET, MPI_COMM_WORLD, &request);
			}
		}
    INFO() << format("shot %d, elapsed time %fs") % is % timer.elapsed();
  }

  INFO() << format("total elapsed time %fs") % totalTimer.elapsed();
#endif

  MPI_Finalize();
  return 0;
}

