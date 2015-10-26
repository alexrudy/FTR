//
//  lqgd.h
//  FTR
//
//  Created by Alexander Rudy on 2015-10-26.
//  Copyright 2015 Alexander Rudy. All rights reserved.
//
#ifndef LQGD_H
#define LQGD_H

#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include "lqg.h"

typedef struct lqgd_differentiator_s * lqgd_differentiator;

/*
Define an LQG in-place differentiator
*/
lqgd_differentiator lqgd_new_differentiator(const int ny, const int nx, double leak);
void lqgd_set_lqg(lqgd_differentiator diff, lqg_filter filter);

// Reset the differntiator
void lqgd_reset_diff(lqgd_differentiator diff);

void lqgd_destroy(lqgd_differentiator diff);

void lqgd_apply(lqgd_differentiator diff, fftw_complex * est_ft);

void lqgd_callback(void * diff, const int ny, const int nx, fftw_complex * est_ft);


#endif /* end of include guard: LQGD_H */
