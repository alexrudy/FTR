#include "ftr.h"
#include "slopemanage.h"
#include "ftr_example.h"
#include <cblas.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

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

double pure_ftr_reconstructor(int nx, int ny, int navg, int iters)
{
  double * sx, * sy, *est;
  struct timespec t_start, t_slopemanage, t_ftr, t_stop;
  double sm_avg, ftr_avg, total;
  int nn = nx * ny;
  int i;
  int *ap;
  
  fftw_complex *gx, *gy;
  ftr_plan plan;
  sm_plan manage_plan;
  
  printf("FTR.\n");
  
  // Allocating plans and memory.
  clock_gettime (CLOCK_REALTIME, &t_start); 
  
  sx = fftw_malloc(nn * sizeof(double));
  memset(sx, (double)1.0, sizeof(double) * nn);
  sy = fftw_malloc(nn * sizeof(double));
  memset(sy, (double)1.0, sizeof(double) * nn);
  est = fftw_malloc(nn * sizeof(double));
  memset(est, (double)0.0, sizeof(double) * nn);
  ap = calloc(nn, sizeof(int));
  
  gx = fftw_malloc(sizeof(fftw_complex) * nn);
  memset(gx, (fftw_complex) 1.0, sizeof(fftw_complex) * nn);
  gy = fftw_malloc(sizeof(fftw_complex) * nn);
  memset(gy, (fftw_complex) 1.0, sizeof(fftw_complex) * nn);
  clock_gettime(CLOCK_REALTIME, &t_stop);
  
  make_aperture(ap, nx, ny);
  
  // Zeroing arrays. We don't time this becasue it isn't really relevant.
  for(i = 0; i < nn; ++i)
  {
      sx[i] = 1.0;
      sy[i] = 1.0;
      est[i] = 0.0;
  }
  
  clock_gettime (CLOCK_REALTIME, &t_start);
  plan = ftr_plan_reconstructor(nx, ny, sx, sy, est);
  ftr_set_filter(plan, gx, gy);
  manage_plan = slope_management_plan(ny, nx, ap);
  clock_gettime(CLOCK_REALTIME, &t_stop);
  
  // Re-zero arrays which might have been touched by the planners.
  for(i = 0; i < nn; ++i)
  {
      sx[i] = (double)rand()/100.0;
      sy[i] = (double)rand()/100.0;
      est[i] = 0.0;
  }
  
  for(i = 0; i < iters; ++i)
  {
      // This part you might do many many times in a for-loop.
      clock_gettime (CLOCK_REALTIME, &t_start); 
      slope_management_execute(manage_plan, sy, sx);
      clock_gettime (CLOCK_REALTIME, &t_slopemanage); 
      ftr_reconstruct(plan);
      clock_gettime (CLOCK_REALTIME, &t_ftr);
      if(i == 1)
      {
          sm_avg = nanoseconds(t_start, t_slopemanage);
          ftr_avg = nanoseconds(t_slopemanage, t_ftr);
      }
      if(i > 1){
        sm_avg = moving_average(sm_avg, nanoseconds(t_start, t_slopemanage), navg);
        ftr_avg = moving_average(ftr_avg, nanoseconds(t_slopemanage, t_ftr), navg);
      }
  }
  
  printf("Averaged %.0f ns for slope management.\n", sm_avg);
  printf("Averaged %.0f ns for ftr.\n", ftr_avg);
  printf("Total iteration time %.0f ns corresponds to a rate of %.0f kHz\n", (sm_avg + ftr_avg), 1e6 / (sm_avg + ftr_avg));
  return (sm_avg + ftr_avg);
}

double pure_vmm_reconstructor(int nx, int ny, int navg, int iters)
{
  double * s, * a, * m;
  int nn = nx * ny;
  int ns = nn * 2;
  int na = nn;
  struct timespec t_start, t_stop;
  double recon_avg;
  int i;
  
  
  printf("VMM.\n");
  
  // Allocating plans and memory.
  clock_gettime (CLOCK_REALTIME, &t_start); 
  
  s = malloc(sizeof(double) * ns);
  memset(s, 1.0, sizeof(double) * ns);
  a = calloc(na, sizeof(double));
  m = malloc(sizeof(double) * na * ns);
  for(i = 0; i < (na * ns); ++i){
    m[i] = (double)rand()/100.0;
  }
  clock_gettime(CLOCK_REALTIME, &t_stop);

  
  for(i = 0; i < iters; ++i)
  {
      // This part you might do many many times in a for-loop.
      clock_gettime (CLOCK_REALTIME, &t_start);
      cblas_dgemv(CblasRowMajor, CblasNoTrans, na, ns, 1.0, m, ns, s, 1, 0.0, a, 1);
      clock_gettime (CLOCK_REALTIME, &t_stop);
      if(i == 1)
      {
          recon_avg = nanoseconds(t_start, t_stop);
      }
      if(i > 1){
        recon_avg = moving_average(recon_avg, nanoseconds(t_start, t_stop), navg);
      }
  }
  printf("Averaged %.0f ns for matrix multiply.\n", recon_avg);
  printf("Total iteration time %.0f ns corresponds to a rate of %.0f kHz\n", recon_avg, 1e6 / recon_avg);
  return recon_avg;
}

double dual_vmm_reconstructor(int nx, int ny, int navg, int iters)
{
  double * s, *p, * a, * m, * ma;
  int nn = nx * ny;
  int ns = nn * 2;
  int na = nn;
  struct timespec t_start, t_first, t_stop;
  double recon_avg, apply_avg, total;
  int i;
  
  
  printf("Dual VMM.\n");
  
  // Allocating plans and memory.
  clock_gettime (CLOCK_REALTIME, &t_start); 
  
  s = malloc(sizeof(double) * ns);
  memset(s, 1.0, sizeof(double) * ns);
  p = calloc(na, sizeof(double));
  a = calloc(na, sizeof(double));
  m = malloc(sizeof(double) * na * ns);
  for(i = 0; i < (na * ns); ++i){
    m[i] = (double)rand()/100.0;
  }  
  ma = malloc(sizeof(double) * na * na);
  for(i = 0; i < (na * na); ++i){
    ma[i] = (double)rand()/100.0;
  }  
  clock_gettime(CLOCK_REALTIME, &t_stop);

  
  for(i = 0; i < iters; ++i)
  {
      // This part you might do many many times in a for-loop.
      clock_gettime (CLOCK_REALTIME, &t_start);
      cblas_dgemv(CblasRowMajor, CblasNoTrans, na, ns, 1.0, m, ns, s, 1, 0.0, p, 1);
      clock_gettime (CLOCK_REALTIME, &t_first);
      cblas_dgemv(CblasRowMajor, CblasNoTrans, na, na, 1.0, ma, na, p, 1, 0.0, a, 1);
      clock_gettime (CLOCK_REALTIME, &t_stop);
      if(i == 1)
      {
          recon_avg = nanoseconds(t_start, t_first);
          apply_avg = nanoseconds(t_first, t_stop);
      }
      if(i > 1){
        recon_avg = moving_average(recon_avg, nanoseconds(t_start, t_first), navg);
        apply_avg = moving_average(apply_avg, nanoseconds(t_first, t_stop), navg);
      }
  }
  printf("Averaged %.0f ns for reconstruction matrix multiply.\n", recon_avg);
  printf("Averaged %.0f ns for apply matrix multiply.\n", apply_avg);
  total = apply_avg + recon_avg;
  printf("Total iteration time %.0f ns corresponds to a rate of %.0f kHz\n", total, 1e6 / total);
  return recon_avg;
}

double hybrid_reconstructor(int nx, int ny, int navg, int iters)
{
  double * sx, * sy, * est;
  double * a, * m;
  struct timespec t_start, t_slopemanage, t_ftr, t_stop;
  double sm_avg, ftr_avg, mat_avg;
  int nn = nx * ny;
  int i;
  int *ap;
  int ns = 2 * nn;
  int na = nn;
  
  fftw_complex *gx, *gy;
  ftr_plan plan;
  sm_plan manage_plan;
  
  printf("FTR + VMM.\n");
  
  // Allocating plans and memory.
  clock_gettime (CLOCK_REALTIME, &t_start); 
  
  sx = fftw_malloc(nn * sizeof(double));
  memset(sx, (double)1.0, sizeof(double) * nn);
  sy = fftw_malloc(nn * sizeof(double));
  memset(sy, (double)1.0, sizeof(double) * nn);
  est = fftw_malloc(nn * sizeof(double));
  memset(est, (double)0.0, sizeof(double) * nn);
  ap = calloc(nn, sizeof(int));
  
  gx = fftw_malloc(sizeof(fftw_complex) * nn);
  memset(gx, (fftw_complex) 1.0, sizeof(fftw_complex) * nn);
  gy = fftw_malloc(sizeof(fftw_complex) * nn);
  memset(gy, (fftw_complex) 1.0, sizeof(fftw_complex) * nn);
  
  a = calloc(na, sizeof(double));
  m = malloc(sizeof(double) * na * ns);
  for(i = 0; i < (na * ns); ++i){
    m[i] = (double)rand()/100.0;
  }
  
  clock_gettime(CLOCK_REALTIME, &t_stop);
  make_aperture(ap, nx, ny);
  
  // Zeroing arrays. We don't time this becasue it isn't really relevant.
  for(i = 0; i < nn; ++i)
  {
      gx[i] = (double)rand()/100.0;
      gy[i] = (double)rand()/100.0;
  }
  
  clock_gettime (CLOCK_REALTIME, &t_start); 
  plan = ftr_plan_reconstructor(nx, ny, sx, sy, est);
  ftr_set_filter(plan, gx, gy);
  manage_plan = slope_management_plan(ny, nx, ap);
  clock_gettime(CLOCK_REALTIME, &t_stop);
  
  // Re-zero arrays which might have been touched by the planners.
  for(i = 0; i < nn; ++i)
  {
      sx[i] = (double)rand()/100.0;
      sy[i] = (double)rand()/100.0;
      est[i] = 0.0;
  }
  
  for(i = 0; i < iters; ++i)
  {
      // This part you might do many many times in a for-loop.
      clock_gettime (CLOCK_REALTIME, &t_start); 
      slope_management_execute(manage_plan, sy, sx);
      clock_gettime (CLOCK_REALTIME, &t_slopemanage); 
      ftr_reconstruct(plan);
      clock_gettime (CLOCK_REALTIME, &t_ftr);
      cblas_dgemv(CblasRowMajor, CblasNoTrans, nn, nn, 1.0, m, nn, est, 1, 0.0, a, 1);
      clock_gettime (CLOCK_REALTIME, &t_stop);
      if(i == 1){
        sm_avg = nanoseconds(t_start, t_slopemanage);
        ftr_avg = nanoseconds(t_slopemanage, t_ftr);
        mat_avg = nanoseconds(t_ftr, t_stop);
      }
      if(i > 1){
        sm_avg = moving_average(sm_avg, nanoseconds(t_start, t_slopemanage), navg);
        ftr_avg = moving_average(ftr_avg, nanoseconds(t_slopemanage, t_ftr), navg);
        mat_avg = moving_average(mat_avg, nanoseconds(t_ftr, t_stop), navg);
      }
  }
  
  printf("Averaged %.0f ns for slope management.\n", sm_avg);
  printf("Averaged %.0f ns for ftr.\n", ftr_avg);
  printf("Averaged %.0f ns for vmm.\n", mat_avg);
  printf("Total iteration time %.0f ns corresponds to a rate of %.0f kHz\n", (sm_avg + ftr_avg + mat_avg), 1e6 / (sm_avg + ftr_avg + mat_avg));
  return (sm_avg + ftr_avg + mat_avg);
}

int main (int argc, char const *argv[])
{
  int navg = 1e3;
  double ftr, vmm, hybrid, dual;
  int nx = 36;
  int ny = 36;
  int iters;
  srand((unsigned) time(NULL));
  // ftr_init(12);
  if(argc > 1)
  {
    iters = (int)atof(argv[1]);
  }else{
    iters = 1e5;
  }
  if(argc > 3)
  {
    nx = (int)atoi(argv[2]);
    ny = (int)atoi(argv[3]);
  }
  printf("Simulating grid of %d x %d\n", nx, ny);
  printf("Testing %d iterations\n", iters);
  
  hybrid = hybrid_reconstructor(nx, ny, navg, iters);
  ftr = pure_ftr_reconstructor(nx, ny, navg, iters);
  vmm = pure_vmm_reconstructor(nx, ny, navg, iters);
  dual = dual_vmm_reconstructor(nx, ny, navg, iters);
  
  printf("Summary:\n");
  printf("FTR:      %6.0f ns %.1f kHz\n", ftr, 1e6 / ftr);
  printf("VMM:      %6.0f ns %.1f kHz\n", vmm, 1e6 / vmm);
  printf("Hybrid:   %6.0f ns %.1f kHz\n", hybrid, 1e6 / hybrid);
  printf("Dual VMM: %6.0f ns %.1f kHz\n", dual, 1e6 / dual);
  
  return 1;
}