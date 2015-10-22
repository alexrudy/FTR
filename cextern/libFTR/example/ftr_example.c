//
//  ftr_example.c
//  FTR
//
//  Created by Alexander Rudy on 2015-10-09.
//  Copyright 2015 Alexander Rudy. All rights reserved.
//

#include "ftr.h"
#include "slopemanage.h"
#include "clock.h"
#include "dbg.h"
#include "aperture.h"

int main (int argc, char const *argv[])
{
    int i;
    struct timespec t_start, t_slopemanage, t_ftr;
    double sm_avg, ftr_avg;
    int navg = 100;
    int nx = 10;
    int ny = 10;
    int nn = nx * ny;
    int iters;
    double *sx, *sy, *est;
    aperture ap;
    fftw_complex *gx, *gy;
    ftr_plan plan;
    sm_plan manage_plan;
    
    if(argc > 1)
    {
      iters = (int)atof(argv[1]);
    }else{
      iters = 1e5;
    }
    
    printf("Conducting %d iterations.\n", iters);
    printf("Allocating arrays. ");
    
    sx = calloc(nn, sizeof(double));
    sy = calloc(nn, sizeof(double));
    est = calloc(nn, sizeof(double));
    ap = aperture_create(ny, nx);
    
    gx = fftw_malloc(sizeof(fftw_complex) * nn);
    memset(gx, 0.0, sizeof(fftw_complex) * nn);
    gy = fftw_malloc(sizeof(fftw_complex) * nn);
    memset(gy, 0.0, sizeof(fftw_complex) * nn);
    
    
    // This only needs to be done once to allocate memory, etc.
    printf("Planning reconstruction.\n");
    plan = ftr_plan_reconstructor(nx, ny, sx, sy, est);
    ftr_set_filter(plan, gx, gy);
    manage_plan = slope_management_plan(ny, nx, ap->ap);
    
    for(i = 0; i < iters; ++i)
    {
        // This part you might do many many times in a for-loop.
        clock_gettime (CLOCK_REALTIME, &t_start); 
        if(i == 0)
        {
            printf("Slope Management ");
        }
        slope_management_execute(manage_plan, sy, sx);
        clock_gettime (CLOCK_REALTIME, &t_slopemanage); 
        
        if(i == 0)
        {
            printf("Reconstructing ");
        }
        ftr_reconstruct(plan);
        clock_gettime (CLOCK_REALTIME, &t_ftr);
        if(i == 0)
        {
            sm_avg = nanoseconds(t_start, t_slopemanage);
            ftr_avg = nanoseconds(t_slopemanage, t_ftr);
            printf("Timing ");
        }
        sm_avg = moving_average(sm_avg, nanoseconds(t_start, t_slopemanage), navg);
        ftr_avg = moving_average(ftr_avg, nanoseconds(t_slopemanage, t_ftr), navg);
    }
    printf("\n");
    printf("Averaged %.0f ns for slope management.\n", sm_avg);
    printf("Averaged %.0f ns for ftr.\n", ftr_avg);
    printf("Total iteration time %.0f ns corresponds to a rate of %.0f kHz\n", (sm_avg + ftr_avg), 1e6 / (sm_avg + ftr_avg));
    // Free the memory once at the end.
    printf("Freeing...\n");
    slope_management_destroy(manage_plan);
    ftr_destroy(plan);
    aperture_destroy(ap);
    
    printf("Done!\n");
    return 0;
error:
    return 1;
}