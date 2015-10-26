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

struct lqgd_differentiator_s {
    int ny, nx;
    int nn, nf, nft; // Auxiliary derived dimensions.
    
    fftw_complex *gains, *alphas, *hp_coefficients;
    fftw_complex *holder, *past, *end; // Internal variables used for Kalman integration.
};
