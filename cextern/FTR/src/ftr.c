//
//  ftr.c
//  FTR
//
//  Created by Alexander Rudy on 2015-08-13.
//  Copyright 2015 Alexander Rudy. All rights reserved.
//

#include "ftr.h"
#include <stdlib.h>
#define FTR_PRECOMUTE FFTW_ESTIMATE

struct ftr_plan_s {
  int nx, ny, nn, nft, nf;
  int *ift, *iconj, *ifs;
  double *sx, *sy, *est;
  fftw_complex *sx_ft, *sy_ft, *est_ft;
  fftw_complex *gx_ft, *gy_ft, *gd_ft;
  fftw_plan p_sx;  // Forward x slope transform.
  fftw_plan p_sy;  // Forward y slope transform.
  fftw_plan p_est; // Inverse phase transform.
};

/* 
This is a private method, decalared here for use inside
ftr_plan_reconstructor.
*/
void ftr_allocate_fftw_plans(ftr_plan recon);

ftr_plan ftr_plan_reconstructor(int ny, int nx, double *sx, double *sy, double *est) {
  size_t i;
  int x, y, io, nfo;
  ftr_plan recon;
  recon = malloc(sizeof(struct ftr_plan_s));
  /* Dimensions of the arrays. */
  recon->nx = nx;
  recon->ny = ny;
  recon->nn = nx * ny;
  recon->nf = (nx / 2) + 1;
  recon->nft = ny * recon->nf;
  nfo = recon->nx % 2 == 0 ? recon->nf - 1 : recon->nf;
  
  
  /* Index array */
  recon->ift = malloc(sizeof(int) * recon->nn);
  recon->ifs = malloc(sizeof(int) * recon->nft);
  recon->iconj = malloc(sizeof(int) * recon->nn);
  
  /* Input and output arrays. */
  recon->sx = sx;
  recon->sy = sy;
  recon->est = est;
  
  /* Complex data arrays which are allocated but not initialized. */
  recon->est_ft = fftw_malloc(sizeof(fftw_complex) * recon->nn);
  recon->sx_ft  = fftw_malloc(sizeof(fftw_complex) * recon->nn);
  recon->sy_ft  = fftw_malloc(sizeof(fftw_complex) * recon->nn);

  recon->gx_ft  = fftw_malloc(sizeof(fftw_complex) * recon->nn);
  recon->gy_ft  = fftw_malloc(sizeof(fftw_complex) * recon->nn);
  recon->gd_ft  = fftw_malloc(sizeof(fftw_complex) * recon->nn);
  
  for(i = 0; i < recon->nn; ++i)
  {
    recon->ift[i] = -1;
  }
  /*
  Create the indexing arrays for unpacking complex transform results (which are hermetian)
  TODO: This might not matter, since the other direction proceeds to ignore it anyways!
  */
  for(i = 0; i < recon->nft; ++i)
  {
    y = i / recon->nf;
    x = i % recon->nf;
    io = (y * recon->nx) + x;
    
    recon->ift[io] = i;
    recon->ifs[i] = io;
    recon->iconj[io] = 0;
    
    if(y == 0 && x > 0 && x < nfo)
    {
      io = (recon->nx - x);
      recon->ift[io] = i;
      recon->iconj[io] = 1;
    }
    if(y > 0 && x > 0 && x < nfo)
    {
      io = (recon->nx - x) + ((recon->ny - y) * recon->nx);
      recon->ift[io] = i;
      recon->iconj[io] = 1;
      
    }
  }
  
  ftr_allocate_fftw_plans(recon);
  return recon;
}

void 
ftr_allocate_fftw_plans(ftr_plan recon) {
  recon->p_sx  = fftw_plan_dft_r2c_2d(recon->nx, recon->ny, recon->sx, recon->sx_ft, FTR_PRECOMUTE);
  recon->p_sy  = fftw_plan_dft_r2c_2d(recon->nx, recon->ny, recon->sy, recon->sy_ft, FTR_PRECOMUTE);
  recon->p_est = fftw_plan_dft_c2r_2d(recon->nx, recon->ny, recon->est_ft, recon->est, FTR_PRECOMUTE);
}

void
ftr_set_filter(ftr_plan recon, fftw_complex *gx, fftw_complex *gy) {
  
  size_t i;
  int x, y;
  double denom;
  for(i = 0; i < recon->nn; ++i)
  {
    recon->gx_ft[i] = gx[i];
    recon->gy_ft[i] = gy[i];
    denom = (pow(cabs(recon->gx_ft[i]), 2) + pow(cabs(recon->gy_ft[i]), 2)) * (fftw_complex)(recon->nn / 2.0);
    recon->gd_ft[i] = (denom > 0.0) ? (1.0 / denom) : 1.0 + 0.0 * I;
  }
}

void
ftr_reconstruct(ftr_plan recon) {
  
  fftw_execute(recon->p_sx);
  fftw_execute(recon->p_sy);
  ftr_estimate(recon);
  fftw_execute(recon->p_est);
}

void ftr_estimate(ftr_plan recon) {
  size_t i, j;
  int x, y;
  fftw_complex numerx, numery, numer, est;
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
    fftw_destroy_plan(recon->p_sx);
    fftw_destroy_plan(recon->p_sy);
    fftw_destroy_plan(recon->p_est);
    
    fftw_free(recon->est_ft);
    fftw_free(recon->sx_ft);
    fftw_free(recon->sy_ft);
    fftw_free(recon->gx_ft);
    fftw_free(recon->gy_ft);
    fftw_free(recon->gd_ft);
    
    free(recon->ift);
    free(recon->ifs);
    free(recon->iconj);
    
    free(recon);
    return;
}

