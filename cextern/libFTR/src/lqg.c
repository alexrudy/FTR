//
//  lqg.c
//  libFTR
//  
//  linear-quadratic-guassian filters implemented in c.
//

#include "lqg.h"
#include "dbg.h"
#include "ftr.h"
#include <stdlib.h>


struct lqg_filter_s {
    int nl, ny, nx;
    int nn, nf, nft; // Auxiliary derived dimensions.
    int *ift, *ifs;
    fftw_complex *gains, *alphas, *hp_coefficients;
    fftw_complex *holder, *past, *end; // Internal variables used for Kalman integration.
};

lqg_filter lqg_new_filter(const int nl, const int ny, const int nx, 
    fftw_complex * gains, fftw_complex * alphas, fftw_complex * hp_coefficients)
{
    lqg_filter filter;
    filter = calloc(1, sizeof(struct lqg_filter_s));
    check_mem(filter);
    
    filter->nl = nl;
    filter->ny = ny;
    filter->nx = nx;
    filter->nn = nx * ny;
    filter->nf = (nx / 2) + 1;
    filter->nft = ny * filter->nf;
    
    /* Index array */
    filter->ift = calloc(filter->nn, sizeof(int));
    check_mem(filter->ift);
    filter->ifs = calloc(filter->nft, sizeof(int));
    check_mem(filter->ifs);
    
    /* Integrator state arrays */
    filter->holder = fftw_malloc(sizeof(fftw_complex) * filter->nn);
    check_mem(filter->holder);
    filter->past = fftw_malloc(sizeof(fftw_complex) * filter->nn * filter->nl);
    check_mem(filter->past);
    filter->end = fftw_malloc(sizeof(fftw_complex) * filter->nn);
    check_mem(filter->end);
    
    filter->gains = gains;
    filter->alphas = alphas;
    filter->hp_coefficients = hp_coefficients;
    
    ftr_map_half_complex(ny, nx, filter->ift, filter->ifs);
    
error:
    return NULL;
}

void lqg_reset(lqg_filter filter)
{
    memset(filter->holder, (fftw_complex) 0.0, sizeof(fftw_complex) * filter->nn);
    memset(filter->past, (fftw_complex) 0.0, sizeof(fftw_complex) * filter->nn * filter->nl);
    memset(filter->end, (fftw_complex) 0.0, sizeof(fftw_complex) * filter->nn);
}

void lqg_apply_filter(lqg_filter filter, fftw_complex * est_ft)
{
    size_t l, i;
    fftw_complex *gains, *alphas, *past;
    memset(filter->holder, (fftw_complex) 0.0, sizeof(fftw_complex) * filter->nn);
    
    // Compute the integrator for each layer.
    //TODO: Could flip these loops around to be more efficient...
    for(l = 0; l < filter->nl; ++l)
    {
        gains = filter->gains + (l * filter->nn);
        alphas = filter->alphas + (l * filter->nn);
        past = filter->past + (l * filter->nn);
        
        //TODO: Could this be nf?
        for(i = 0; i < filter->nn; ++i)
        {
            past[i] = est_ft[i] * gains[i] + past[i] * alphas[i];
            filter->holder[i] += past[i];
        }
    }
    for(i = 0; i < filter->nn; ++i)
    {
        est_ft[i] = filter->holder[i] - filter->end[i] * filter->hp_coefficients[i];
        filter->end[i] = est_ft[i];
    }
    return;
}

void lqg_filter_callback(void * filter, const int ny, const int nx, fftw_complex * est_ft)
{
    lqg_apply_filter((lqg_filter) filter, est_ft);
    return;
}