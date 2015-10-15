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

// Type definition for the FTR Plan.
// FTR Plans are opaque structures managed
// by the functions below.
typedef struct ftr_plan_s * ftr_plan;


// Functions

// Create an FTR plan.
ftr_plan
ftr_plan_reconstructor(int nx, int ny, double *sx, double *sy, double *est);

// Set the filter associated with an FTR plan.
void
ftr_set_filter(ftr_plan recon, fftw_complex *gx, fftw_complex *gy);

// Execute the plan, acting on the arrays specified when this plan
// was constructed.
void
ftr_reconstruct(ftr_plan recon);

// Destroy the plan, cleaning up memory.
void
ftr_destroy(ftr_plan recon);

// Perform only the FTR estimation step on the transformed arrays.
void 
ftr_estimate(ftr_plan recon);

// Utilties
void ftr_map_half_complex(int ny, int nx, int * map, int * imap);

#endif /* end of include guard: FTR_H_D3963E38 */

