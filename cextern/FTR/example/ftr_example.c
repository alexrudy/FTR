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
    int i;
    for(i = 0; i < nn; ++i)
    {   
        pt[i] = 0.0;
    }
}

void make_aperture(int * ap, int nx, int ny)
{
    int i, j;
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
    return ((average * (double)(navg - 1)) + duration) / (double)navg;
}

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
    manage_plan = slope_management_plan(ny, nx, ap);
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
        slope_management_execute(manage_plan, sy, sx);
        clock_gettime (CLOCK_REALTIME, &t_slopemanage); 
        
        if(i == 0)
        {
            printf("Reconstructing...\n");
        }
        ftr_reconstruct(plan);
        clock_gettime (CLOCK_REALTIME, &t_ftr);
        if(i == 0)
        {
            sm_avg = nanoseconds(t_start, t_slopemanage);
            ftr_avg = nanoseconds(t_slopemanage, t_ftr);
            printf("Timing...\n");
        }
        sm_avg = moving_average(sm_avg, nanoseconds(t_start, t_slopemanage), navg);
        ftr_avg = moving_average(ftr_avg, nanoseconds(t_slopemanage, t_ftr), navg);
        if(i == 0)
        {
          printf("Iter1 Done...\n");
        }
    }
    printf("Done with loop...\n");
    printf("Averaged %.0f ns for slope management.\n", sm_avg);
    printf("Averaged %.0f ns for ftr.\n", ftr_avg);
    printf("Total iteration time %.0f ns corresponds to a rate of %.0f kHz\n", (sm_avg + ftr_avg), 1e6 / (sm_avg + ftr_avg));
    // Free the memory once at the end.
    printf("Freeing...\n");
    slope_management_destroy(manage_plan);
    ftr_destroy(plan);
    
    printf("Done!\n");
    return 0;
}