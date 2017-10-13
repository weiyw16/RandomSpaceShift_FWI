/* Triangle smoothing */
/*
  Copyright (C) 2004 University of Texas at Austin
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include "rsf.h"
#include "triangle.h"

#include "_bool.h"
/*^*/

#include "alloc.h"
#include "blas.h"

#ifndef _sf_triangle_h

typedef struct sf_Triangle *sf_triangle;
/* abstract data type */
/*^*/

#endif

struct sf_Triangle {
    float *tmp;
    int np, nb, nx;
};

static void fold (int o, int d, int nx, int nb, int np, 
		  const float *x, float* tmp);
static void fold2 (int o, int d, int nx, int nb, int np, 
		   float *x, const float* tmp);
static void doubint (int nx, float *x, bool der);
static void triple (int o, int d, int nx, int nb, 
		    float* x, const float* tmp, bool box);
static void triple2 (int o, int d, int nx, int nb, const float* x, float* tmp, bool box);

sf_triangle sf_triangle_init (int nbox /* triangle length */, 
			      int ndat /* data length */)
/*< initialize >*/
{
    sf_triangle tr;

    tr = (sf_triangle) sf_alloc(1,sizeof(*tr));

    tr->nx = ndat;
    tr->nb = nbox;
    tr->np = ndat + 2*nbox;
    
    tr->tmp = sf_floatalloc(tr->np);

    return tr;
}

static void fold (int o, int d, int nx, int nb, int np, 
		  const float *x, float* tmp)
{
    int i, j;

    /* copy middle */
    for (i=0; i < nx; i++) 
	tmp[i+nb] = x[o+i*d];
    
    /* reflections from the right side */
    for (j=nb+nx; j < np; j += nx) {
	for (i=0; i < nx && i < np-j; i++)
	    tmp[j+i] = x[o+(nx-1-i)*d];
	j += nx;
	for (i=0; i < nx && i < np-j; i++)
	    tmp[j+i] = x[o+i*d];
    }
    
    /* reflections from the left side */
    for (j=nb; j >= 0; j -= nx) {
	for (i=0; i < nx && i < j; i++)
	    tmp[j-1-i] = x[o+i*d];
	j -= nx;
	for (i=0; i < nx && i < j; i++)
	    tmp[j-1-i] = x[o+(nx-1-i)*d];
    }
}

static void fold2 (int o, int d, int nx, int nb, int np, 
		   float *x, const float* tmp)
{
    int i, j;

    /* copy middle */
    for (i=0; i < nx; i++) 
	x[o+i*d] = tmp[i+nb];

    /* reflections from the right side */
    for (j=nb+nx; j < np; j += nx) {
	for (i=0; i < nx && i < np-j; i++)
	    x[o+(nx-1-i)*d] += tmp[j+i];
	j += nx;
	for (i=0; i < nx && i < np-j; i++)
	    x[o+i*d] += tmp[j+i];
    }
    
    /* reflections from the left side */
    for (j=nb; j >= 0; j -= nx) {
	for (i=0; i < nx && i < j; i++)
	    x[o+i*d] += tmp[j-1-i];
	j -= nx;
	for (i=0; i < nx && i < j; i++)
	    x[o+(nx-1-i)*d] += tmp[j-1-i];
    }
}
    
static void doubint (int nx, float *xx, bool der)
{
    int i;
    float t;

    /* integrate backward */
    t = 0.;
    for (i=nx-1; i >= 0; i--) {
	t += xx[i];
	xx[i] = t;
    }

    if (der) return;

    /* integrate forward */
    t=0.;
    for (i=0; i < nx; i++) {
	t += xx[i];
	xx[i] = t;
    }
}

static void doubint2 (int nx, float *xx, bool der)
{
    int i;
    float t;


    /* integrate forward */
    t=0.;
    for (i=0; i < nx; i++) {
	t += xx[i];
	xx[i] = t;
    }

    if (der) return;

    /* integrate backward */
    t = 0.;
    for (i=nx-1; i >= 0; i--) {
	t += xx[i];
	xx[i] = t;
    }
}

static void triple (int o, int d, int nx, int nb, float* x, const float* tmp, bool box)
{
    int i;
    const float *tmp1, *tmp2;
    float wt;
    
    if (box) {
	tmp2 = tmp + 2*nb;

	wt = 1.0/(2*nb-1);
	for (i=0; i < nx; i++) {
	    x[o+i*d] = (tmp[i+1] - tmp2[i])*wt;
	}
    } else {
	tmp1 = tmp + nb;
	tmp2 = tmp + 2*nb;

	wt = 1.0/(nb*nb);
	for (i=0; i < nx; i++) {
	    x[o+i*d] = (2.*tmp1[i] - tmp[i] - tmp2[i])*wt;
	}
    }
}

static void triple2 (int o, int d, int nx, int nb, const float* x, float* tmp, bool box)
{
    int i;
    float wt;

    for (i=0; i < nx + 2*nb; i++) {
	tmp[i] = 0;
    }

    if (box) {
	wt = 1.0/(2*nb-1);

	cblas_saxpy(nx,  +wt,x+o,d,tmp+1   ,1);
	cblas_saxpy(nx,  -wt,x+o,d,tmp+2*nb,1);
    } else {
	wt = 1.0/(nb*nb);
    
	cblas_saxpy(nx,  -wt,x+o,d,tmp     ,1);
	cblas_saxpy(nx,2.*wt,x+o,d,tmp+nb  ,1);
	cblas_saxpy(nx,  -wt,x+o,d,tmp+2*nb,1);
    }
}

void sf_smooth (sf_triangle tr  /* smoothing object */, 
		int o, int d    /* trace sampling */, 
		bool der        /* if derivative */, 
		bool box        /* if box filter */,
		float *x        /* data (smoothed in place) */)
/*< apply triangle smoothing >*/
{
    fold (o,d,tr->nx,tr->nb,tr->np,x,tr->tmp);
    doubint (tr->np,tr->tmp,(bool) (box || der));
    triple (o,d,tr->nx,tr->nb,x,tr->tmp,box);
}

void sf_smooth2 (sf_triangle tr  /* smoothing object */, 
		 int o, int d    /* trace sampling */, 
		 bool der        /* if derivative */,
		 bool box        /* if box filter */,
		 float *x        /* data (smoothed in place) */)
/*< apply adjoint triangle smoothing >*/
{
    triple2 (o,d,tr->nx,tr->nb,x,tr->tmp,box);
    doubint2 (tr->np,tr->tmp,(bool) (box || der));
    fold2 (o,d,tr->nx,tr->nb,tr->np,x,tr->tmp);
}

void  sf_triangle_close(sf_triangle tr)
/*< free allocated storage >*/
{
    free (tr->tmp);
    free (tr);
}

/* 	$Id: triangle.c 7107 2011-04-10 02:04:14Z ivlad $	 */
int sf_first_index (int i          /* dimension [0...dim-1] */, 
		    int j        /* line coordinate */, 
		    int dim        /* number of dimensions */, 
		    const int *n /* box size [dim] */, 
		    const int *s /* step [dim] */)
/*< Find first index for multidimensional transforms >*/
{
    int i0, n123, ii;
    int k;

    n123 = 1;
    i0 = 0;
    for (k=0; k < dim; k++) {
	if (k == i) continue;
	ii = (j/n123)%n[k]; /* to cartesian */
	n123 *= n[k];	
	i0 += ii*s[k];      /* back to line */
    }

    return i0;
}

int sf_filedims (sf_file file, /*@out@*/ int *n) 
/*< Find file dimensions.
--- 
Outputs the number of dimensions dim and a dimension array n[dim] >*/
{
    int i, dim;
    char key[3];

    dim = 1;
    for (i=0; i < SF_MAX_DIM; i++) {
	(void) snprintf(key,3,"n%d",i+1);
	if (!sf_histint(file,key,n+i)) {
	    n[i]=1;
	} else if (n[i] > 1) {
	    dim=i+1;
	}
    }
    return dim;
}

int smooth(float *origin, int nz, int nx, int rectz, int rectx)
{
    int dim, dim1, i, j, n[SF_MAX_DIM], rect[SF_MAX_DIM], s[SF_MAX_DIM];
    int nrep, irep, n1, n2, i2, i0;
    bool adj, diff[SF_MAX_DIM], box[SF_MAX_DIM];
    char key[6];
    float* data;
    sf_triangle tr;
    sf_file in, out;

    //cbw
    //in  = sf_input ("g2.rsf");
    //out = sf_output ("out.rsf");
    //sf_putint(out, "n1", nz);
    //sf_putint(out, "n2", nx);

    //if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");

    //cbw
    //dim = sf_filedims (in,n);
    dim = 2;
    n[0] = nz;
    n[1] = nx;

    dim1 = -1;
    for (i=0; i < dim; i++) {
	snprintf(key,6,"rect%d",i+1);
	if (!sf_getint(key,rect+i)) rect[i]=1;
  rect[0] = rectz;
  rect[1] = rectx;
	/*( rect#=(1,1,...) smoothing radius on #-th axis )*/ 
	if (rect[i] > 1) dim1 = i;
	snprintf(key,6,"diff%d",i+1);
	if (!sf_getbool(key,diff+i)) diff[i]=false;
	/*( diff#=(n,n,...) differentiation on #-th axis )*/
	snprintf(key,5,"box%d",i+1);
	if (!sf_getbool(key,box+i)) box[i]=false;
	/*( box#=(n,n,...) box (rather than triangle) on #-th axis )*/
    }


    n1 = n2 = 1;
    for (i=0; i < dim; i++) {
	if (i <= dim1) {
	    s[i] = n1;
	    n1 *= n[i];
	} else {
	    n2 *= n[i];
	}
    }

    data = sf_floatalloc (n1);

    if (!sf_getint("repeat",&nrep)) nrep=1;
    /* repeat filtering several times */

    if (!sf_getbool("adj",&adj)) adj=false;
    /* run in the adjoint mode */

    for (i2=0; i2 < n2; i2++) {
      //cbw
    //sf_floatread(data,n1,in);
    //memcpy(data, &origin[n1 * i2], sizeof(float) * n1);
    data = &origin[n1 * i2];

	for (i=0; i <= dim1; i++) {
	    if (rect[i] <= 1) continue;
	    tr = sf_triangle_init (rect[i],n[i]);
	    for (j=0; j < n1/n[i]; j++) {
		i0 = sf_first_index (i,j,dim1+1,n,s);
		for (irep=0; irep < nrep; irep++) {
		    if (adj) {
			sf_smooth (tr,i0,s[i],diff[i],box[i],data);
		    } else {
			sf_smooth2 (tr,i0,s[i],diff[i],box[i],data);
		    }
		}
	    }
	    sf_triangle_close(tr);
	}
	
	//sf_floatwrite(data,n1,out);
    }    
}

