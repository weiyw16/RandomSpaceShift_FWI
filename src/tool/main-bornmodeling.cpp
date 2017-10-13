
extern "C" {
#include <rsf.h>
#include "fdutil.h"
#include "Mbandpass.h"
}

#ifdef _OPENMP
#include <omp.h>
#endif

#include "mpi.h"
#include "logger.h"
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
  sf_file shots;
  sf_file shotsborn;
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
  vreal=sf_input ("vreal");   /* real velocity model, unit=m/s */
  shots=sf_output("shots");
  shotsborn=sf_output("shotsborn");

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

  sf_putint(shots,"n1",nt);
  sf_putint(shots,"n2",ng);
  sf_putint(shots,"n3",ns);
  sf_putfloat(shots,"d1",dt);
  sf_putfloat(shots,"d2",jgx*dx);
  sf_putfloat(shots,"o1",0);
  sf_putfloat(shots,"o2",0);
  sf_putstring(shots,"label1","Time");
  sf_putstring(shots,"label2","Lateral");
  sf_putstring(shots,"label3","Shot");
  sf_putstring(shots,"unit1","sec");
  sf_putstring(shots,"unit2","m");
  sf_putfloat(shots,"amp",amp);
  sf_putfloat(shots,"fm",fm);
  sf_putint(shots,"ng",ng);
  sf_putint(shots,"szbeg",szbeg);
  sf_putint(shots,"sxbeg",sxbeg);
  sf_putint(shots,"gzbeg",gzbeg);
  sf_putint(shots,"gxbeg",gxbeg);
  sf_putint(shots,"jsx",jsx);
  sf_putint(shots,"jsz",jsz);
  sf_putint(shots,"jgx",jgx);
  sf_putint(shots,"jgz",jgz);
  sf_putint(shots, "nb", nb);
  sf_putint(shots, "free", freeSurface);

  sf_putint(shotsborn,"n1",nt);
  sf_putint(shotsborn,"n2",ng);
  sf_putint(shotsborn,"n3",ns);
  sf_putfloat(shotsborn,"d1",dt);
  sf_putfloat(shotsborn,"d2",jgx*dx);
  sf_putfloat(shotsborn,"o1",0);
  sf_putfloat(shotsborn,"o2",0);
  sf_putstring(shotsborn,"label1","Time");
  sf_putstring(shotsborn,"label2","Lateral");
  sf_putstring(shotsborn,"label3","Shot");
  sf_putstring(shotsborn,"unit1","sec");
  sf_putstring(shotsborn,"unit2","m");
  sf_putfloat(shotsborn,"amp",amp);
  sf_putfloat(shotsborn,"fm",fm);
  sf_putint(shotsborn,"ng",ng);
  sf_putint(shotsborn,"szbeg",szbeg);
  sf_putint(shotsborn,"sxbeg",sxbeg);
  sf_putint(shotsborn,"gzbeg",gzbeg);
  sf_putint(shotsborn,"gxbeg",gxbeg);
  sf_putint(shotsborn,"jsx",jsx);
  sf_putint(shotsborn,"jsz",jsz);
  sf_putint(shotsborn,"jgx",jgx);
  sf_putint(shotsborn,"jgz",jgz);
  sf_putint(shotsborn, "nb", nb);
  sf_putint(shotsborn, "free", freeSurface);
  Velocity v = SfVelocityReader::read(vinit, nx, nz);

  float vmin = *std::min_element(v.dat.begin(), v.dat.end());
  float vmax = *std::max_element(v.dat.begin(), v.dat.end());
  sf_putfloat(shots, "vmin", vmin);
  sf_putfloat(shots, "vmax", vmax);
  sf_putfloat(shotsborn, "vmin", vmin);
  sf_putfloat(shotsborn, "vmax", vmax);

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
	sprintf(logfile, "fm-damp-%02d.log", params.rank);
  FILELog::setLogFile(logfile);
#ifndef USE_SW
	printGitInfo();
#endif

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
  std::vector<float> bcoff(exvel.nz * exvel.nx, 0.0);
  //bcoff = fmMethod.getBornCoff();
  //bcoff = fmMethod.getBornCoff(exvel, exvel_real, params.dx, params.dt);

  INFO() << "init velocity and coff " ;

  fmMethod.bindVelocity(exvel);
  fmMethod.bindRealVelocity(exvel_real);
  //bcoff = fmMethod.getBornCoff();
  //bcoff = fmMethod.getBornCoff(exvel, exvel_real, params.dx);
  bcoff = fmMethod.getBornCoff(exvel, exvel_real, params.dx, params.dt);
  fmMethod.bindBornCoff(bcoff);
  INFO() << "get born coff" ;

  std::vector<float> wlt(nt);
  rickerWavelet(&wlt[0], nt, fm, dt, params.amp);

  std::vector<float> dobs(params.ntask * params.nt * params.ng, 0);
  std::vector<float> dobs_born(params.ntask * params.nt * params.ng, 0);
  for(int is=rank*k; is<rank*k+ntask; is++) {
    int local_is = is - rank * k;
    Timer timer;
    std::vector<float> p0(exvel.nz * exvel.nx, 0);
    std::vector<float> p1(exvel.nz * exvel.nx, 0);
    std::vector<float> p2(exvel.nz * exvel.nx, 0);
    std::vector<float> pborn0(exvel.nz * exvel.nx, 0);
    std::vector<float> pborn1(exvel.nz * exvel.nx, 0);
    std::vector<float> dobs_trans(params.nt * params.ng, 0);
    std::vector<float> dobs_trans_born(params.nt * params.ng, 0);
	INFO() << "define p0 p1 pborn0 pborn1" ;
    ShotPosition curSrcPos = allSrcPos.clipRange(is, is);
		/*
		int dn = 10;
		sf_file fullwv = sf_output("fullwv.rsf");
		sf_putint(fullwv, "n1" , exvel.nz);
		sf_putint(fullwv, "n2" , exvel.nx);
		sf_putint(fullwv, "n3" , nt / dn);
		*/

		//fmMethod.initFdUtil(params.vinit, &exvel, nb, params.dx, dt);
    for(int it=0; it<nt; it++) {
		INFO() << format("start itetration %d ") % it ;
      fmMethod.addSource(&pborn1[0], &wlt[it], curSrcPos);
      fmMethod.stepbornForward(pborn0, pborn1);
      //fmMethod.stepForward(p0, p1, 0);
      std::swap(pborn1, pborn0);
      fmMethod.recordSeis(&dobs_trans_born[it*ng], &pborn0[0]);
      //fmMethod.addSource(&p1[0], &pborn1[0], curSrcPos);
	  vectorMinus(p1, pborn1, p2);
      fmMethod.stepForward(p0, p2);
      //fmMethod.stepForward(p0, p1, 0);
      swap3(p1, p0, p2);
      fmMethod.recordSeis(&dobs_trans[it*ng], &p1[0]);
			/*
			if(it % dn == 0)
				sf_floatwrite(&p0[0], exvel.nx * exvel.nz, fullwv);
				*/
    }
		//exit(1);
    matrix_transpose(&dobs_trans[0], &dobs[local_is * ng * nt], ng, nt);
    matrix_transpose(&dobs_trans_born[0], &dobs_born[local_is * ng * nt], ng, nt);

		//fmMethod.fwiRemoveDirectArrival(&dobs[local_is * ng * nt], local_is);
		if(np == 1) {
			sf_floatwrite(&dobs[local_is * ng * nt], ng*nt, params.shots);
			sf_floatwrite(&dobs_born[local_is * ng * nt], ng*nt, params.shotsborn);
		}
		else {
			if(rank == 0) {
				sf_floatwrite(&dobs[local_is * ng * nt], ng*nt, params.shots);
				sf_floatwrite(&dobs_born[local_is * ng * nt], ng*nt, params.shotsborn);
				if(is == rank * k + ntask - 1) {
					for(int other_is = rank * k + ntask ; other_is < ns ; other_is ++) {
						MPI_Recv(&dobs[0], ng*nt, MPI_FLOAT, other_is / k, other_is, MPI_COMM_WORLD, &status);
						sf_floatwrite(&dobs[0], ng*nt, params.shots);
						MPI_Recv(&dobs_born[0], ng*nt, MPI_FLOAT, other_is / k, other_is, MPI_COMM_WORLD, &status);
						sf_floatwrite(&dobs_born[0], ng*nt, params.shotsborn);
					}
				}
			}
			else {
				MPI_Isend(&dobs[local_is * ng * nt], ng*nt, MPI_FLOAT, 0, is, MPI_COMM_WORLD, &request);
				MPI_Isend(&dobs_born[local_is * ng * nt], ng*nt, MPI_FLOAT, 0, is, MPI_COMM_WORLD, &request);
			}
		}
    INFO() << format("shot %d, elapsed time %fs") % is % timer.elapsed();
  }

  INFO() << format("total elapsed time %fs") % totalTimer.elapsed();

  MPI_Finalize();
  return 0;
}

