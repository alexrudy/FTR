#include "ftr.h"
#include "slopemanage.h"
#include "ftr_example.h"
#include "dbg.h"
#include "aperture.h"
#include <cblas.h>
#include <openblas_config.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define cfree(N) if(N) free(N)
#define ffree(N) if(N) fftw_free(N)

double pure_ftr_reconstructor(int nx, int ny, int navg, int iters)
{
  double * sx, * sy, *est;
  struct timespec t_start, t_slopemanage, t_ftr, t_stop;
  double sm_avg = 0.0, ftr_avg = 0.0, total = 0.0;
  int nn = nx * ny;
  int i;
  int *ap;
  
  fftw_complex *gx, *gy;
  ftr_plan plan;
  sm_plan manage_plan;
  
  printf("--> FTR.\n");
  
  // Allocating plans and memory.
  clock_gettime (CLOCK_REALTIME, &t_start); 
  
  sx = fftw_malloc(nn * sizeof(double));
  check_mem(sx);
  memset(sx, (double)1.0, sizeof(double) * nn);
  sy = fftw_malloc(nn * sizeof(double));
  check_mem(sy);
  memset(sy, (double)1.0, sizeof(double) * nn);
  est = fftw_malloc(nn * sizeof(double));
  check_mem(est);
  memset(est, (double)0.0, sizeof(double) * nn);
  ap = calloc(nn, sizeof(int));
  check_mem(ap);
  
  gx = fftw_malloc(sizeof(fftw_complex) * nn);
  check_mem(gx);
  memset(gx, (fftw_complex) 1.0, sizeof(fftw_complex) * nn);
  gy = fftw_malloc(sizeof(fftw_complex) * nn);
  check_mem(gy);
  memset(gy, (fftw_complex) 1.0, sizeof(fftw_complex) * nn);
  clock_gettime(CLOCK_REALTIME, &t_stop);
  
  make_aperture(ny, nx, ap);
  
  // Zeroing arrays. We don't time this becasue it isn't really relevant.
  for(i = 0; i < nn; ++i)
  {
      sx[i] = 1.0;
      sy[i] = 1.0;
      est[i] = 0.0;
  }
  
  clock_gettime (CLOCK_REALTIME, &t_start);
  plan = ftr_plan_reconstructor(nx, ny, sx, sy, est);
  check(plan, "allocating plan");
  ftr_set_filter(plan, gx, gy);
  manage_plan = slope_management_plan(ny, nx, ap);
  check(manage_plan, "allocating slope management");
  clock_gettime(CLOCK_REALTIME, &t_stop);
  
  // Re-zero arrays which might have been touched by the planners.
  for(i = 0; i < nn; ++i)
  {
      sx[i] = (double)rand()/100.0;
      sy[i] = (double)rand()/100.0;
      est[i] = 0.0;
  }
  
  // Warm up.
  clock_gettime(CLOCK_REALTIME, &t_start);
  for(i = 0; i < iters; ++i)
  {
      slope_management_execute(manage_plan, sy, sx);
      ftr_reconstruct(plan);
  }
  clock_gettime(CLOCK_REALTIME, &t_stop);
  printf("Long-Average %.0f ns.\n", nanoseconds(t_start, t_stop) / iters);
  
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
error:
  if(plan) ftr_destroy(plan);
  if(manage_plan) slope_management_destroy(manage_plan);
  ffree(sx);
  ffree(sy);
  ffree(est);
  ffree(gx);
  ffree(gy);
  ffree(ap);
  return (sm_avg + ftr_avg);
}

double pure_vmm_reconstructor(int nx, int ny, int navg, int iters)
{
  double * s, * a, * m;
  int nn = nx * ny;
  int ns, *ap;
  int na = 32 * 32;
  struct timespec t_start, t_stop;
  double recon_avg = 0.0;
  int i;
  
  
  printf("--> VMM.\n");
  
  // Allocating plans and memory.
  
  // Make an aperture so we know how many slopes to use.
  ap = malloc(sizeof(int) * nn);
  check_mem(ap);
  ns = 2*make_aperture(ny, nx, ap);
  
  // Allocate Slopes
  s = malloc(sizeof(double) * ns);
  check_mem(s);
  memset(s, 1.0, sizeof(double) * ns);
  
  // Allocate Actuators
  a = calloc(na, sizeof(double));
  check_mem(a);
  
  // Allocate Matrix
  m = malloc(sizeof(double) * na * ns);
  check_mem(m);
  for(i = 0; i < (na * ns); ++i){
    m[i] = (double)rand()/100.0;
  }
    
  // Warm up
  clock_gettime(CLOCK_REALTIME, &t_start);
  for(i = 0; i < iters; ++i)
  {
      cblas_dgemv(CblasRowMajor, CblasNoTrans, na, ns, 1.0, m, ns, s, 1, 0.0, a, 1);
  }
  clock_gettime(CLOCK_REALTIME, &t_stop);
  printf("Long-Average %.0f ns.\n", nanoseconds(t_start, t_stop) / iters);
  
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
  printf("Averaged %.0f ns for matrix multiply (%dx%d).\n", recon_avg, na, ns);
  printf("Total iteration time %.0f ns corresponds to a rate of %.0f kHz\n", recon_avg, 1e6 / recon_avg);
error:
  cfree(s);
  cfree(a);
  cfree(m);
  return recon_avg;
}

double dual_vmm_reconstructor(int nx, int ny, int navg, int iters)
{
  double * s, *p, * a, * m, * ma;
  int nn = nx * ny;
  int na = 32 * 32;
  int ns, *ap;
  int nm = (nx / 2 + 1) * ny;
  struct timespec t_start, t_first, t_stop;
  double recon_avg = 0.0, apply_avg = 0.0, total = 0.0;
  int i;
  
  
  printf("--> Dual VMM.\n");
  
  // Allocating plans and memory.
  clock_gettime (CLOCK_REALTIME, &t_start); 
  
  // Make an aperture so we know how many slopes to use.
  ap = malloc(sizeof(int) * nn);
  check_mem(ap);
  ns = 2*make_aperture(ny, nx, ap);
  
  s = malloc(sizeof(double) * ns);
  check_mem(s);
  for(i = 0; i < ns; ++i){
    s[i] = (double)rand()/100.0;
  }
  p = calloc(nm, sizeof(double));
  check_mem(p);
  a = calloc(na, sizeof(double));
  check_mem(a);
  m = malloc(sizeof(double) * nm * ns);
  check_mem(m);
  for(i = 0; i < (nm * ns); ++i){
    m[i] = (double)rand()/100.0;
  }  
  ma = malloc(sizeof(double) * na * nm);
  check_mem(ma);
  for(i = 0; i < (na * nm); ++i){
    ma[i] = (double)rand()/100.0;
  }
  clock_gettime(CLOCK_REALTIME, &t_stop);
  clock_gettime (CLOCK_REALTIME, &t_start);
  for(i = 0; i < iters; ++i)
  {
      cblas_dgemv(CblasRowMajor, CblasNoTrans, nm, ns, 1.0, m, ns, s, 1, 0.0, p, 1);
      cblas_dgemv(CblasRowMajor, CblasNoTrans, na, nm, 1.0, ma, nm, p, 1, 0.0, a, 1);
  }
  clock_gettime(CLOCK_REALTIME, &t_stop);
  printf("Long-Average %.0f ns.\n", nanoseconds(t_start, t_stop) / iters);
  
  for(i = 0; i < iters; ++i)
  {
      // This part you might do many many times in a for-loop.
      clock_gettime (CLOCK_REALTIME, &t_start);
      cblas_dgemv(CblasRowMajor, CblasNoTrans, nm, ns, 1.0, m, ns, s, 1, 0.0, p, 1);
      clock_gettime (CLOCK_REALTIME, &t_first);
      cblas_dgemv(CblasRowMajor, CblasNoTrans, na, nm, 1.0, ma, nm, p, 1, 0.0, a, 1);            
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
  printf("Averaged %.0f ns for reconstruction matrix multiply (%dx%d).\n", recon_avg, nm, ns);
  printf("Averaged %.0f ns for apply matrix multiply (%dx%d).\n", apply_avg, na, nm);
  total = apply_avg + recon_avg;
  printf("Total iteration time %.0f ns corresponds to a rate of %.0f kHz\n", total, 1e6 / total);
error:
  cfree(s);
  cfree(p);
  cfree(a);
  cfree(m);
  cfree(ma);
  return total;
}

double hybrid_reconstructor(int nx, int ny, int navg, int iters)
{
  double * sx, * sy, * est;
  double * a, * m;
  struct timespec t_start, t_slopemanage, t_ftr, t_stop;
  double sm_avg = 0.0, ftr_avg = 0.0, mat_avg = 0.0;
  int nn = nx * ny;
  int i;
  int ns, *ap;
  int na = 32 * 32;
  int nm = (nx / 2 + 1) * ny;
  
  fftw_complex *gx, *gy;
  ftr_plan plan;
  sm_plan manage_plan;
  
  printf("--> FTR + VMM.\n");
  
  // Allocating plans and memory.
  clock_gettime (CLOCK_REALTIME, &t_start); 
  
  sx = fftw_malloc(nn * sizeof(double));
  check_mem(sx);
  memset(sx, (double)1.0, sizeof(double) * nn);
  sy = fftw_malloc(nn * sizeof(double));
  check_mem(sy);
  memset(sy, (double)1.0, sizeof(double) * nn);
  est = fftw_malloc(nn * sizeof(double));
  check_mem(est);
  memset(est, (double)0.0, sizeof(double) * nn);
  ap = calloc(nn, sizeof(int));
  check_mem(ap);
  ns = 2*make_aperture(ny, nx, ap);
  print_aperture(ny, nx, ap);
  
  gx = fftw_malloc(sizeof(fftw_complex) * nn);
  memset(gx, (fftw_complex) 1.0, sizeof(fftw_complex) * nn);
  check_mem(gx);
  gy = fftw_malloc(sizeof(fftw_complex) * nn);
  memset(gy, (fftw_complex) 1.0, sizeof(fftw_complex) * nn);
  check_mem(gy);
  
  a = calloc(na, sizeof(double));
  m = malloc(sizeof(double) * na * nm);
  for(i = 0; i < (na * nm); ++i){
    m[i] = (double)rand()/100.0;
  }
  
  clock_gettime(CLOCK_REALTIME, &t_stop);
  
  // Zeroing arrays. We don't time this becasue it isn't really relevant.
  for(i = 0; i < nn; ++i)
  {
      gx[i] = (double)rand()/100.0;
      gy[i] = (double)rand()/100.0;
  }
  
  clock_gettime (CLOCK_REALTIME, &t_start); 
  plan = ftr_plan_reconstructor(nx, ny, sx, sy, est);
  check(plan, "Allocating FTR Plan");
  ftr_set_filter(plan, gx, gy);
  manage_plan = slope_management_plan(ny, nx, ap);
  check(manage_plan, "Allocating SM Plan");
  clock_gettime(CLOCK_REALTIME, &t_stop);
  
  // Re-zero arrays which might have been touched by the planners.
  for(i = 0; i < nn; ++i)
  {
      sx[i] = (double)rand()/100.0;
      sy[i] = (double)rand()/100.0;
      est[i] = 0.0;
  }
  clock_gettime (CLOCK_REALTIME, &t_start); 
  for(i = 0; i < iters; ++i)
  {
      slope_management_execute(manage_plan, sy, sx);
      ftr_reconstruct(plan);
      cblas_dgemv(CblasRowMajor, CblasNoTrans, na, nm, 1.0, m, nm, est, 1, 0.0, a, 1);
  }
  clock_gettime (CLOCK_REALTIME, &t_stop);
  printf("Long-Average %.0f ns.\n", nanoseconds(t_start, t_stop) / iters);
  
  for(i = 0; i < iters; ++i)
  {
      // This part you might do many many times in a for-loop.
      clock_gettime (CLOCK_REALTIME, &t_start); 
      slope_management_execute(manage_plan, sy, sx);
      clock_gettime (CLOCK_REALTIME, &t_slopemanage); 
      ftr_reconstruct(plan);
      clock_gettime (CLOCK_REALTIME, &t_ftr);
      cblas_dgemv(CblasRowMajor, CblasNoTrans, na, nm, 1.0, m, nm, est, 1, 0.0, a, 1);
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
  printf("Averaged %.0f ns for vmm (%dx%d).\n", mat_avg, na, nm);
  printf("Total iteration time %.0f ns corresponds to a rate of %.0f kHz\n", (sm_avg + ftr_avg + mat_avg), 1e6 / (sm_avg + ftr_avg + mat_avg));
error:
  if(plan) ftr_destroy(plan);
  if(manage_plan) slope_management_destroy(manage_plan);
  ffree(sx);
  ffree(sy);
  ffree(est);
  ffree(ap);
  ffree(gx);
  ffree(gy);
  cfree(a);
  cfree(m);
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
  dual = dual_vmm_reconstructor(nx, ny, navg, iters);  
  hybrid = hybrid_reconstructor(nx, ny, navg, iters);
  ftr = pure_ftr_reconstructor(nx, ny, navg, iters);
  vmm = pure_vmm_reconstructor(nx, ny, navg, iters);
  
  printf("Summary:\n");
  printf("FTR:      %6.0f ns %.1f kHz\n", ftr, 1e6 / ftr);
  printf("VMM:      %6.0f ns %.1f kHz\n", vmm, 1e6 / vmm);
  printf("Hybrid:   %6.0f ns %.1f kHz\n", hybrid, 1e6 / hybrid);
  printf("Dual VMM: %6.0f ns %.1f kHz\n", dual, 1e6 / dual);
  
  return 1;
}