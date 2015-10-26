//
//  ftr.c
//  FTR
//
//  Created by Alexander Rudy on 2015-08-13.
//  Copyright 2015 Alexander Rudy. All rights reserved.
//

#include "ftr.h"
#include <stdlib.h>
#include "dbg.h"

#ifndef FTR_PRECOMUTE
#define FTR_PRECOMUTE FFTW_PATIENT
#endif

struct ftr_plan_s {
  int nx, ny, nn, nft, nf;
  int *ift, *ifs;
  double *sx, *sy, *est;
  fftw_complex *sx_ft, *sy_ft, *est_ft;
  fftw_complex *gx_ft, *gy_ft, *gd_ft;
  fftw_plan p_sx;  // Forward x slope transform.
  fftw_plan p_sy;  // Forward y slope transform.
  fftw_plan p_est; // Inverse phase transform.
};

int ftr_init(int nthreads)
{
  int rc;
  rc = fftw_init_threads();
  check(rc != 0, "FFTW Threads did not initialize properly.");
  fftw_plan_with_nthreads(nthreads);
  
  return 0;
error:
  return 1;
}

ftr_halfcomplex ftr_halfcomplex_map(const int ny, const int nx)
{
    ftr_halfcomplex hcmap;
    int nn, nft;
    hcmap = calloc(1, sizeof(struct ftr_halfcomplex_s));
    check_mem(hcmap);
    
    hcmap->ny = ny;
    hcmap->nx = nx;
    hcmap->nf = (nx / 2) + 1;
    nn = ny * nx;
    nft = ny * hcmap->nf;
    
    hcmap->f2hc = calloc(nn, sizeof(int));
    check_mem(hcmap->f2hc);
    hcmap->hc2f = calloc(nft, sizeof(int));
    check_mem(hcmap->hc2f);
    
    ftr_map_half_complex(ny, nx, hcmap->f2hc, hcmap->hc2f);
    return hcmap;
    
error:
    ftr_halfcomplex_destroy(hcmap);
    return NULL;
}

void ftr_halfcomplex_destroy(ftr_halfcomplex hcmap)
{
    if(hcmap){
        if(hcmap->f2hc) free(hcmap->f2hc);
        if(hcmap->hc2f) free(hcmap->hc2f);
        free(hcmap);
    }
    return;
}

void ftr_map_half_complex(int ny, int nx, int * map, int * imap)
{
  int x, y;
  int nn, nf, nft;
  int i_full, i_half;
  size_t i;
  
  // Compute dimensions.
  nn = nx * ny;
  nf = (nx / 2) + 1;
  nft = nf * ny;
  
  for(i = 0; i < nft; ++i)
  {
    y = i / nf;
    x = i % nf;
    
    i_half = (y * nf) + x;
    i_full = (y * nx) + x;
    
    map[i_full] = i_half;
    imap[i_half] = i_full;
    if(y == 0 && x > 0 && x < nf - 1)
    {
      i_full = (nx - x);
      map[i_full] = i_half;
    }
    if(y > 0 && x > 0 && x < nf - 1)
    {
      i_full = (nx - x) + ((ny - y) * nx);
      map[i_full] = i_half;
    }
  }
}

ftr_plan ftr_plan_reconstructor(int ny, int nx, double *sx, double *sy, double *est) {
  int nfo;
  ftr_plan recon;
  recon = calloc(1, sizeof(struct ftr_plan_s));
  check_mem(recon);
  
  /* Dimensions of the arrays. */
  recon->nx = nx;
  recon->ny = ny;
  recon->nn = nx * ny;
  recon->nf = (nx / 2) + 1;
  recon->nft = ny * recon->nf;
  nfo = recon->nx % 2 == 0 ? recon->nf - 1 : recon->nf;
  
  
  /* Index array */
  recon->ift = calloc(recon->nn, sizeof(int));
  check_mem(recon->ift);
  recon->ifs = calloc(recon->nft, sizeof(int));
  check_mem(recon->ifs);
  
  /* Input and output arrays. */
  recon->sx = sx;
  recon->sy = sy;
  recon->est = est;
  
  /* Complex data arrays which are allocated but not initialized. */
  recon->est_ft = fftw_malloc(sizeof(fftw_complex) * recon->nn);
  check_mem(recon->est_ft);
  
  recon->sx_ft  = fftw_malloc(sizeof(fftw_complex) * recon->nn);
  check_mem(recon->sx_ft);
  
  recon->sy_ft  = fftw_malloc(sizeof(fftw_complex) * recon->nn);
  check_mem(recon->sy_ft);
  
  recon->gx_ft  = fftw_malloc(sizeof(fftw_complex) * recon->nn);
  check_mem(recon->gx_ft);
  
  recon->gy_ft  = fftw_malloc(sizeof(fftw_complex) * recon->nn);
  check_mem(recon->gy_ft);
  
  recon->gd_ft  = fftw_malloc(sizeof(fftw_complex) * recon->nn);
  check_mem(recon->gd_ft);
  
  /* Prepare the mapping for the half-complex problem. */
  ftr_map_half_complex(ny, nx, recon->ift, recon->ifs);
  
  /* Allocate FFTW plans. */
  recon->p_sx  = fftw_plan_dft_r2c_2d(recon->ny, recon->nx, recon->sx, recon->sx_ft, FTR_PRECOMUTE);
  check(recon->p_sx, "Failed to plan sx FFT r2c 2d");
  recon->p_sy  = fftw_plan_dft_r2c_2d(recon->ny, recon->nx, recon->sy, recon->sy_ft, FTR_PRECOMUTE);
  check(recon->p_sx, "Failed to plan sy FFT r2c 2d");
  recon->p_est = fftw_plan_dft_c2r_2d(recon->ny, recon->nx, recon->est_ft, recon->est, FTR_PRECOMUTE);
  check(recon->p_est, "Failed to plan est FFT c2r 2d");
  
  return recon;

error:
  ftr_destroy(recon);
  return NULL;
}

void
ftr_set_filter(ftr_plan recon, const fftw_complex *gx, const fftw_complex *gy) {
  
  size_t i;
  double denom;
  memcpy(recon->gx_ft, gx, sizeof(fftw_complex) * recon->nn);
  memcpy(recon->gy_ft, gy, sizeof(fftw_complex) * recon->nn);
  for(i = 0; i < recon->nn; ++i)
  {
    denom = (pow(cabs(recon->gx_ft[i]), 2) + pow(cabs(recon->gy_ft[i]), 2)) * (fftw_complex)(recon->nn / 2.0);
    // Note the normalization factor on the end here. The factor of two accounts for the two FFTs involved
    // in the forward transform, with only one FFT involved in the backward transform.
    
    // This ensures that we never divide by zero, so long as the filter sums to zero.
    recon->gd_ft[i] = (denom > 0.0) ? (1.0 / denom) : 1.0 + 0.0 * I;
  }
}

#ifndef FTR_RECONSTRUCT_AS_MACRO
void
ftr_reconstruct(ftr_plan recon)
{
    ftr_reconstruct_with_callback(recon, NULL, NULL);
    return;
}
#endif

void
ftr_forward_transform(ftr_plan recon)
{
  // Forward slope transforms.
  fftw_execute(recon->p_sx);
  fftw_execute(recon->p_sy);
}

void
ftr_reconstruct_with_callback(ftr_plan recon, ftr_estimate_callback callback, void * data)
{
    // Forward slope transforms.
    fftw_execute(recon->p_sx);
    fftw_execute(recon->p_sy);
  
    // Estimation using the gx/gy filters.
    ftr_apply_filter(recon);
    
    // Filtering callback for the estimate.
    if(callback) callback(data, recon->ny, recon->nx, recon->est_ft);
    
    // Inverse transform to compute results.
    fftw_execute(recon->p_est);
    return; 
}

void
ftr_apply_callback(ftr_plan recon, ftr_estimate_callback callback, void * data)
{
    callback(data, recon->ny, recon->nx, recon->est_ft);
    return;
}

void
ftr_backward_transform(ftr_plan recon)
{
  fftw_execute(recon->p_est);
}

void ftr_apply_filter(ftr_plan recon) {
  size_t i;
  for(i = 0; i < recon->nft; ++i)
  {
      recon->est_ft[i] = (conj(recon->gx_ft[recon->ifs[i]]) * recon->sx_ft[i]
        + conj(recon->gy_ft[recon->ifs[i]]) * recon->sy_ft[i])
        * recon->gd_ft[recon->ifs[i]];
  }
  return;
}

void
ftr_destroy(ftr_plan recon)
{
    if(recon){
        if(recon->p_sx) fftw_destroy_plan(recon->p_sx);
        if(recon->p_sy) fftw_destroy_plan(recon->p_sy);
        if(recon->p_est) fftw_destroy_plan(recon->p_est);
    
        if(recon->est_ft) fftw_free(recon->est_ft);
        if(recon->sx_ft) fftw_free(recon->sx_ft);
        if(recon->sy_ft) fftw_free(recon->sy_ft);
        if(recon->gx_ft) fftw_free(recon->gx_ft);
        if(recon->gy_ft) fftw_free(recon->gy_ft);
        if(recon->gd_ft) fftw_free(recon->gd_ft);
    
        if(recon->ift) free(recon->ift);
        if(recon->ifs) free(recon->ifs);
        
        free(recon);
    }
    return;
}

