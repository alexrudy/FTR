//
//  lqgd.c
//  FTR
//
//  Created by Alexander Rudy on 2015-10-26.
//  Copyright 2015 Alexander Rudy. All rights reserved.
//

#include "lqgd.h"
#include "dbg.h"
#include "ftr.h"
#include <stdlib.h>


struct lqgd_differentiator_s {
    int ny, nx;
    int nn; // Auxiliary derived dimensions.
    double leak;
    fftw_complex *past;
    lqg_filter filter;
};

lqgd_differentiator lqgd_new_differentiator(const int ny, const int nx, double leak)
{
    lqgd_differentiator diff;
    diff = calloc(1, sizeof(struct lqgd_differentiator_s));
    check_mem(diff);
    
    diff->ny = ny;
    diff->nx = nx;
    diff->nn = ny * nx;
    diff->leak = leak;
    diff->filter = NULL;

    diff->past = malloc(diff->nn * sizeof(fftw_complex));
    check_mem(diff->past);
    lqgd_reset_diff(diff);
    
    return diff;
error:
    lqgd_destroy(diff);
    return NULL;
}

void lqgd_destroy(lqgd_differentiator diff)
{
    if(diff){
        if(diff->past) free(diff->past);
        free(diff);
    }
}

void lqgd_set_lqg(lqgd_differentiator diff, lqg_filter filter)
{
    diff->filter = filter;
}

void lqgd_reset_diff(lqgd_differentiator diff)
{
    memset(diff->past, (fftw_complex) 0.0, sizeof(fftw_complex) * diff->nn);
}

void lqgd_apply(lqgd_differentiator diff, fftw_complex * est_ft)
{
    size_t i;
    if(diff->filter){
        lqg_apply_filter(diff->filter, est_ft);
    }
    for(i = 0; i < diff->nn; ++i)
    {
        est_ft[i] = est_ft[i] - diff->past[i] * diff->leak;
        diff->past[i] = est_ft[i];
    }
    return;
}

void lqgd_callback(void * diff, const int ny, const int nx, fftw_complex * est_ft)
{
    lqgd_apply((lqgd_differentiator) diff, est_ft);
    return;
}