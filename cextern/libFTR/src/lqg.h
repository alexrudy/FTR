//
//  lgq.h
//  FTR
//
//  linear-quadratic-guassian filters implemented in c.
//

#ifndef LGQ_H
#define LGQ_H

#include <math.h>
#include <complex.h>
#include <fftw3.h>

typedef struct lqg_filter_s *lqg_filter;

/*
Define a new LQG Filter.
*/
lqg_filter lqg_new_filter(const int nl, const int ny, const int nx, 
    fftw_complex * gains, fftw_complex * alphas, fftw_complex * hp_coefficients);

/*
Apply the LQG filter to an estimate of the phase, modifying the phase in place.
*/
void lqg_apply_filter(lqg_filter filter, fftw_complex * est_ft);

// Callback function which correctly applies the filter using the above technique.
void lqg_filter_callback(void * filter, const int ny, const int nx, fftw_complex * est_ft);

// Reset the interal state of the filter.
void lqg_reset(lqg_filter filter);

/*
Destroy an LQG Filter.
*/
void lqg_destroy(lqg_filter filter);

#endif /* end of include guard: LGQ_H */


