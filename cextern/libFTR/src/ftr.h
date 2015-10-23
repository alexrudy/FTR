//
//  ftr.h
//  FTR
//
//  Created by Alexander Rudy on 2015-08-13.
//  Copyright 2015 Alexander Rudy. All rights reserved.
//

#ifndef FTR_H_D3963E38
#define FTR_H_D3963E38

#define FTR_PRECOMUTE FFTW_MEASURE
#define FTR_RECONSTRUCT_AS_MACRO

// Inclue FFTW, and Complex to use native complex data types.
#include <math.h>
#include <complex.h>
#include <fftw3.h>

// Type definition for the FTR Plan.
// FTR Plans are opaque structures managed
// by the functions below.
typedef struct ftr_plan_s * ftr_plan;

// Define the callback function.
typedef void (*ftr_estimate_callback)(void * data, fftw_complex * est_ft);


// Functions
// Initialize multithreaded things.
int ftr_init(int nthreads);


// Create an FTR plan.
ftr_plan
ftr_plan_reconstructor(int nx, int ny, double *sx, double *sy, double *est);

// Set the filter associated with an FTR plan.
void
ftr_set_filter(ftr_plan recon, const fftw_complex *gx, const fftw_complex *gy);

// Execute the plan, acting on the arrays specified when this plan
// was constructed.
#ifdef FTR_RECONSTRUCT_AS_MACRO
#define ftr_reconstruct(R) ftr_reconstruct_with_callback(R, NULL, NULL)
#else
void
ftr_reconstruct(ftr_plan recon);
#endif

void
ftr_reconstruct_with_callback(ftr_plan recon, ftr_estimate_callback callback, void * data);

// Destroy the plan, cleaning up memory.
void
ftr_destroy(ftr_plan recon);

// Perform only the FTR estimation step on the transformed arrays.
void 
ftr_apply_filter(ftr_plan recon);

void
ftr_forward_transform(ftr_plan recon);

void
ftr_backward_transform(ftr_plan recon);

void
ftr_apply_callback(ftr_plan recon, ftr_estimate_callback callback, void * data);

// Utilties
void
ftr_map_half_complex(int ny, int nx, int * map, int * imap);

#endif /* end of include guard: FTR_H_D3963E38 */

