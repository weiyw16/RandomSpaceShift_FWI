/* this is for dot product test 
 * written by weiyw at 2017/09/27
 * contact at weiyw17@gmail.com
 */
extern "C"{
#include <rsf.h>
}

#include <new>

float *getMemory(unsigned long size){
	float *p = new float[size]();
	if(!p) printf("get p memory failed! \n");
	return p;
}

int main(int argc, char *argv[]){
	/* initivate */
	sf_init(argc, argv);

	/* parameters */
	sf_file Filem1;
	sf_file Filem2;
	sf_file Filed1;
	sf_file Filed2;
	int nx, nz, nt, ng;
	int i, j, k;
	//float d_result, m_result;
	float sum_buff=0, sum=0;
	/*
	//float **m1buff, **m2buff, **d1buff, **d2buff, **mresult, **dresult;	
	//float *m1buff, *m2buff, *d1buff, *d2buff, *mresult, *dresult;	
	float *m1buff = new float[nz*nx]();
	float *m2buff = new float[nz*nx]();
	float *d1buff = new float[nt*ng]();
	float *d2buff = new float[nt*ng]();
	float *mresult = new float[nz*nx]();
	float *dresult = new float[nt*ng]();
	*/
	
	Filem1 = sf_input("m1");
	Filem2 = sf_input("m2");
	Filed1 = sf_input("d1");
	Filed2 = sf_input("d2");
	
	if(!sf_histint(Filem1, "n1", &nz))  sf_error("No n1= in input");
	if(!sf_histint(Filem1, "n2", &nx))  sf_error("No n2= in input");
	if(!sf_histint(Filed1, "n1", &nt))  sf_error("No n1= in input");
	if(!sf_histint(Filed1, "n2", &ng))  sf_error("No n2= in input");
	
	printf("get nx, nz, nt, ng done \n");
	printf("get nx, nz, nt, ng done: nx = %d \n", nx);
	printf("get nx, nz, nt, ng done: nz = %d \n", nz);
	printf("get nx, nz, nt, ng done: nt = %d \n", nt);
	printf("get nx, nz, nt, ng done: ng = %d \n", ng);

	/*
	float *m1buff = new float[nz*nx]();
	float *m2buff = new float[nz*nx]();
	float *d1buff = new float[nt*ng]();
	float *d2buff = new float[nt*ng]();
	float *mresult = new float[nz*nx]();
	float *dresult = new float[nt*ng]();
	*/
	
	float *m1buff = getMemory(nz*nx);
	float *m2buff = getMemory(nz*nx);
	float *d1buff = getMemory(nt*ng);
	float *d2buff = getMemory(nt*ng);
	float *mresult = getMemory(nx*nx);
	float *dresult = getMemory(ng*ng);


	//float *m1_buff[nz*nx], *m2_buff[nz*nx], *d1_buff[nt*ng], *d2_buff[nt*ng];
	/*
	m1buff = sf_floatalloc2(nz,nx);
	m2buff = sf_floatalloc2(nz,nx);
	d1buff = sf_floatalloc2(nt,ng);
	d2buff = sf_floatalloc2(nt,ng);
	mresult = sf_floatalloc2(nx,nx);
	dresult = sf_floatalloc2(ng,ng);	// ng should be equal with nx */ 
	//sf_floatread(m1buff, nz*nx, Filem1);
	//sf_floatread(m2buff, nz*nx, Filem2);
	//sf_floatread(d1buff, nt*ng, Filed1);
	//sf_floatread(d2buff, nt*ng, Filed2);
	
	for(i = 0; i < nx*nz; i++){
		m1buff[i] = 0;
		m2buff[i] = 0;
	}

	for(i = 0; i < nt*ng; i++){
		d1buff[i] = 0;
		d2buff[i] = 0;
	}

	for(i = 0; i < ng*ng; i++){
		dresult[i] = 0;
	}

	for(i = 0; i < nx*nx; i++){
		mresult[i] = 0;
	}

	sf_floatread(m1buff, nz*nx, Filem1);
	sf_floatread(m2buff, nz*nx, Filem2);

	printf("read form sf-file done \n");

	for(i = 0; i < nx; i++){
		for(j = 0; j < nx; j++){
			for(k = 0; k < nz; k++){
				mresult[j * nx + i] += m1buff[i * nz + k] * m2buff[j * nz + k];
				//mresult[i][j] += m1buff[k][j] * m2buff[k][i];
			}}}
	printf("get mresult done \n");

	delete []m1buff;
	delete []m2buff;
	sf_floatread(d1buff, nt*ng, Filed1);
	sf_floatread(d2buff, nt*ng, Filed2);

	for(i = 0; i < ng; i++){
		for(j = 0; j < ng; j++){
			for(k = 0; k < nt; k++){
				dresult[j * ng + i] += d1buff[i * nt + k] * d2buff[j * nt + k];
				//dresult[i][j] += d1buff[k][j] * d2buff[k][i];
			}}}
	printf("get dresult done \n");

	for(i = 0; i < ng; i++){
	for(j = 0; j < ng; j++){
		sum_buff =  dresult[i * ng + j] - mresult[i * ng + j];
		if(sum_buff < 0 ) sum_buff = 0 - sum_buff;
		sum += sum_buff;
	}
	}

	printf("sum = %f \n", sum);

	//delete []m1buff;
	//delete []m2buff;
	delete []d1buff;
	delete []d2buff;
	delete []mresult;
	delete []dresult;

	sf_close();
	return 0;

}
