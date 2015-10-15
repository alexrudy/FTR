//
//  ftr.c
//  FTR
//
//  Created by Alexander Rudy on 2015-08-13.
//  Copyright 2015 Alexander Rudy. All rights reserved.
//

#include "ftr.h"
#include <stdlib.h>
#define FTR_PRECOMUTE FFTW_MEASURE

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

/* 
This is a private method, decalared here for use inside
ftr_plan_reconstructor.
*/
void ftr_allocate_fftw_plans(ftr_plan recon);

//TODO: This is probably the broken piece when nx or ny are odd.
//I'm not sure why it is dependent on even or odd behavior, but hey?
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
  
  /* Prepare the mapping for the half-complex problem. */
  ftr_map_half_complex(ny, nx, recon->ift, recon->ifs);
  
  /* Allocate FFTW plans. */
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
    // Note the normalization factor on the end here. The factor of two accounts for the two FFTs involved
    // in the forward transform, with only one FFT involved in the backward transform.
    
    // This ensures that we never divide by zero, so long as the filter sums to zero.
    recon->gd_ft[i] = (denom > 0.0) ? (1.0 / denom) : 1.0 + 0.0 * I;
  }
}

void
ftr_reconstruct(ftr_plan recon)
{
    ftr_reconstruct_with_callback(recon, NULL, NULL);
    return;
}  


void
ftr_reconstruct_with_callback(ftr_plan recon, ftr_estimate_callback callback, void * data)
{
    // Forward slope transforms.
    fftw_execute(recon->p_sx);
    fftw_execute(recon->p_sy);
  
    // Estimation using the gx/gy filters.
    ftr_estimate(recon);
    if(callback) callback(data, recon->est_ft);
    fftw_execute(recon->p_est);
    return; 
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
        
    free(recon);
    return;
}

