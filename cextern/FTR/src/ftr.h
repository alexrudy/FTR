//
//  ftr.h
//  FTR
//
//  Created by Alexander Rudy on 2015-08-13.
//  Copyright 2015 Alexander Rudy. All rights reserved.
//

#ifndef FTR_H_D3963E38
#define FTR_H_D3963E38

// Inclue FFTW, and Complex to use native complex data types.
#include <math.h>
#include <complex.h>
#include <fftw3.h>


// Type definitions
typedef struct ftr_plan ftr_plan;

// Functions
ftr_plan
ftr_plan_reconstructor(int nx, int ny, double *sx, double *sy, double *est);

void
ftr_allocate_fftw_plans(ftr_plan recon);

void
ftr_set_filter(ftr_plan recon, fftw_complex *gx, fftw_complex *gy);

void
ftr_reconstruct(ftr_plan recon);

#endif /* end of include guard: FTR_H_D3963E38 */

