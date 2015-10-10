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

ftr_plan ftr_plan_reconstructor(int nx, int ny, double *sx, double *sy, double *est) {
  ftr_plan recon;
  recon = malloc(sizeof(struct ftr_plan_s));
  /* Dimensions of the arrays. */
  recon->nx = nx;
  recon->ny = ny;
  
  /* Input and output arrays. */
  recon->sx = sx;
  recon->sy = sy;
  recon->est = est;
  
  /* Complex data arrays which are allocated but not initialized. */
  recon->est_ft = fftw_malloc(sizeof(fftw_complex) * nx * ny);
  recon->sx_ft  = fftw_malloc(sizeof(fftw_complex) * nx * ny);
  recon->sy_ft  = fftw_malloc(sizeof(fftw_complex) * nx * ny);
  recon->gd_ft  = fftw_malloc(sizeof(fftw_complex) * nx * ny);
  
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
  double denom;
  recon->gx_ft = gx;
  recon->gy_ft = gy;
  
  for(i = 0; i < (recon->nx * recon->ny); ++i)
  {
    denom = pow(cabs(recon->gx_ft[i]), 2.0) + pow(cabs(recon->gy_ft[i]),2.0);
    recon->gd_ft[i] = (denom > 0.0) ? (1.0 / denom) : 1.0;
  }
}

void
ftr_reconstruct(ftr_plan recon) {
  int nn;
  size_t i;
  nn = (recon->nx * recon->ny);
  
  fftw_execute(recon->p_sx);
  fftw_execute(recon->p_sy);
  
  for(i = 0; i < nn; ++i)
  {
    recon->est_ft[i] = (conj(recon->gx_ft[i]) * recon->sx_ft[i] + conj(recon->gy_ft[i]) * recon->sy_ft[i]) * recon->gd_ft[i];
  }

  fftw_execute(recon->p_est);
}

