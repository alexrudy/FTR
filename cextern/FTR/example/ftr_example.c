//
//  ftr_example.c
//  FTR
//
//  Created by Alexander Rudy on 2015-10-09.
//  Copyright 2015 Alexander Rudy. All rights reserved.
//

#include "ftr.h"
#include "slopemanage.h"
#include "ftr_example.h"

void zero(double * pt, int nn)
{
    size_t i;
    for(i = 0; i < nn; ++i)
    {   
        pt[i] = 0.0;
    }
}

void make_aperture(int * ap, int nx, int ny)
{
    size_t i, j;
    for(i = 0; i < nx; ++i)
    {
        for(j = 0; j < ny; ++j)
        {
            if(i == 0 || j == 0 || i == nx - 1 || j == ny - 1) {
                ap[(i * nx) + j] = 0;
            }else{
                ap[(i * nx) + j] = 1;
            }
        }
    }
}

double nanoseconds(const struct timespec start, const struct timespec end)
{
    double total;
    total = (end.tv_sec - start.tv_sec) * 1e9;
    total += (end.tv_nsec - start.tv_nsec);
    return total;
}

double moving_average(double average, double duration, int navg)
{
    return ((average * (navg - 1)) + duration) / navg;
}

int main (int argc, char const *argv[])
{
    size_t i;
    struct timespec t_start, t_slopemanage, t_ftr;
    double slopemanage, ftr;
    int navg = 100;
    
    int nx = 10;
    int ny = 10;
    int nn = nx * ny;
    int iters;
    double *sx, *sy, *est;
    int *ap;
    fftw_complex *gx, *gy;
    ftr_plan plan;
    sm_plan manage_plan;
    
    if(argc > 0)
    {
      iters = (int)atof(argv[1]);
    }else{
      iters = 1e5;
    }
    printf("Conducting %d iterations.\n", iters);
    printf("Allocating...\n");
    
    sx = malloc(sizeof(double) * nn);
    sy = malloc(sizeof(double) * nn);
    est = malloc(sizeof(double) * nn);
    ap = malloc(sizeof(int) * nn);
    
    gx = fftw_malloc(sizeof(fftw_complex) * nn);
    gy = fftw_malloc(sizeof(fftw_complex) * nn);
    
    for(i = 0; i < nn; ++i)
    {
        gx[i] = 0.0;
        gy[i] = 0.0;
    }
    
    zero(sx, nn);
    zero(sy, nn);
    zero(est, nn);
    make_aperture(ap, nx, ny);
    
    // This only needs to be done once to allocate memory, etc.
    printf("Planning...\n");
    plan = ftr_plan_reconstructor(nx, ny, sx, sy, est);
    ftr_set_filter(plan, gx, gy);
    manage_plan = slope_management_plan(nx, ny, ap);
    zero(sx, nn);
    zero(sy, nn);
    
    for(i = 0; i < iters; ++i)
    {
        // This part you might do many many times in a for-loop.
        clock_gettime (CLOCK_REALTIME, &t_start); 
        if(i == 0)
        {
            printf("Slope Management...\n");
        }
        slope_management_execute(manage_plan, sx, sy);
        clock_gettime (CLOCK_REALTIME, &t_slopemanage); 
        
        if(i == 0)
        {
            printf("Reconstructing...\n");
        }
        ftr_reconstruct(plan);
        clock_gettime (CLOCK_REALTIME, &t_ftr);
        if(i == 0)
        {
            slopemanage = nanoseconds(t_start, t_slopemanage);
            ftr = nanoseconds(t_slopemanage, t_ftr);
            printf("Timing...\n");
        }
        slopemanage = moving_average(slopemanage, nanoseconds(t_start, t_slopemanage), navg);
        ftr = moving_average(ftr, nanoseconds(t_slopemanage, t_ftr), navg);
    }
    
    printf("Averaged %.0f ns for slope management.\n", slopemanage);
    printf("Averaged %.0f ns for ftr.\n", ftr);
    printf("Total iteration time %.0f ns corresponds to a rate of %.0f kHz\n", (slopemanage + ftr), 1e6 / (slopemanage + ftr));
    // Free the memory once at the end.
    printf("Freeing...\n");
    slope_management_destroy(manage_plan);
    ftr_destroy(plan);
    
    printf("Done!\n");
    return 0;
}