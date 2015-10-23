#include "ftr.h"
#include "slopemanage.h"
#include "clock.h"
#include "dbg.h"
#include "aperture.h"
// #include <cblas.h>
#include <cblas_openblas.h>
#include <openblas_config.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <math.h>

#define cfree(N) if(N) free(N)
#define ffree(N) if(N) fftw_free(N)

#define MAX_STR_LEN 80
#define NACTUATORS 1024 + 52
#define NCOEFF 14

/* Function to create apertures that are appropriately sized.
Figured this out by trial and error. Had to special case 8.
Might also have to special case 30, but I'm not sure, becasue
Don's documents don't line up with Sri's subaperture extraction.
(Especially w/r/t partially illuminated subapertures).
*/
aperture ShaneAO_aperture(const int ng, const int nacross)
{
  double outer, inner;
  switch (nacross) {
  case 8:
    outer = 3.6;
    break;
  default:
    outer = ((double)nacross / 2.0) - 1.0;
  }
  inner = 0.3 * outer;
  return aperture_create_with_radii(ng, ng, outer, inner);
}

///////////////
// Result types
///////////////

// Typedefs for reconstructor results.
typedef struct reconstructor_result * recon_result;
typedef struct reconstructor_result_part * recon_result_part;
typedef recon_result (*reconstructor)(aperture ap, int navg, int iters);

/*
Individual parts have some timing information. We group them in structs here.
*/
struct reconstructor_result_part {
    double duration;
    char description[MAX_STR_LEN];
};

/*
Reconstructors also have a collection of parts for timing.
*/
struct reconstructor_result {
    int n_parts;
    recon_result_part *parts;
};

recon_result_part rr_new_part(double duration, char * description)
{
    recon_result_part part;
    part = malloc(sizeof(struct reconstructor_result_part));
    check_mem(part);
    part->duration = duration;
    strncpy(part->description, description, MAX_STR_LEN);
error:
    return part;
}

recon_result rr_new(int n_parts, ...)
{
    recon_result result;
    result = malloc(sizeof(struct reconstructor_result));
    result->n_parts = n_parts;
    result->parts = malloc(n_parts * sizeof(struct reconstructor_result_part));
    va_list argp;
    va_start(argp, n_parts);
    for(size_t i = 0; i < n_parts; ++i)
    {
        result->parts[i] = va_arg(argp, recon_result_part);
    }
    va_end(argp);
    return result;

}

void rr_print_result(recon_result result, char * description)
{
    double duration, total;
    size_t j;
    
    total = 0.0;
    for(j = 0; j < result->n_parts; ++j)
    {
        duration = result->parts[j]->duration;
        total += duration;
        printf("* %-40s %6.0f ns %6.1f kHz\n", result->parts[j]->description, duration, 1e6 / duration);
        
    }
    printf("                                           --------- ----------\n");
    
    printf("Total %-36s %6.0f ns %6.1f kHz\n", description, total, 1e6 / total);
    printf("========================================== ========= ==========\n");
    
}

void rr_del_result(recon_result result)
{
    size_t i;
    recon_result_part part;
    
    if(result){
        for(i = 0; i < result->n_parts; ++i)
        {
            part = result->parts[i];
            if(part) free(part);
        }
        free(result);
    }
    return;
}

void shuffle(int *array, size_t n) {    
    srand((unsigned) time(NULL));
    if (n > 1) {
        size_t i;
        for (i = n - 1; i > 0; i--) {
            size_t j = (unsigned int) (((double)rand()/(double)RAND_MAX)*(i+1));
            int t = array[j];
            array[j] = array[i];
            array[i] = t;
        }
    }
}

////////////////////////////
// Reconstruction functions:
////////////////////////////

/*
Pure FTR, with no accounting for alignment etc.
*/
recon_result pure_ftr_reconstructor(aperture ap, int navg, int iters)
{
  double * sx, * sy, *est;
  struct timespec t_start, t_slopemanage, t_ftr, t_stop;
  double sm_avg = 0.0, ftr_avg = 0.0, total = 0.0;
  int nn = ap->nx * ap->ny;
  int i;
  recon_result result;
  recon_result_part sm, ftr;
  fftw_complex *gx, *gy;
  ftr_plan plan;
  sm_plan manage_plan;
  
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
  
  gx = fftw_malloc(sizeof(fftw_complex) * nn);
  check_mem(gx);
  memset(gx, (fftw_complex) 1.0, sizeof(fftw_complex) * nn);
  gy = fftw_malloc(sizeof(fftw_complex) * nn);
  check_mem(gy);
  memset(gy, (fftw_complex) 1.0, sizeof(fftw_complex) * nn);
  clock_gettime(CLOCK_REALTIME, &t_stop);
    
  // Zeroing arrays. We don't time this becasue it isn't really relevant.
  memset(sx, (double) 1.0, sizeof(double) * nn);
  memset(sy, (double) 1.0, sizeof(double) * nn);
  memset(est, (double) 0.0, sizeof(double) * nn);
  
  clock_gettime (CLOCK_REALTIME, &t_start);
  
  // Planning FTR
  plan = ftr_plan_reconstructor(ap->ny, ap->nx, sx, sy, est);
  check(plan, "allocating plan");
  ftr_set_filter(plan, gx, gy);
  manage_plan = slope_management_plan_from_aperture(ap);
  check(manage_plan, "allocating slope management");
  clock_gettime(CLOCK_REALTIME, &t_stop);
  
  // Seed the slopes with random numbers.
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
  log_info("Warm-up average %.0f ns.", nanoseconds(t_start, t_stop) / iters);
  
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
  
  sm = rr_new_part(sm_avg, "Slope Management");
  ftr = rr_new_part(ftr_avg, "FTR");
  result = rr_new(2, sm, ftr);
error:
  if(plan) ftr_destroy(plan);
  if(manage_plan) slope_management_destroy(manage_plan);
  ffree(sx);
  ffree(sy);
  ffree(est);
  ffree(gx);
  ffree(gy);
  return result;
}

recon_result pure_vmm_reconstructor(aperture ap, int navg, int iters)
{
  double * s, * a, * m;
  int nn = ap->nx * ap->ny;
  int ns = 2*ap->ni + NCOEFF;
  int na = NACTUATORS + NCOEFF;
  struct timespec t_start, t_stop;
  double recon_avg = 0.0;
  int i;
  
  recon_result result;
  recon_result_part vmm;
  
  // Allocating plans and memory.
  
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
  log_info("Warm-up average %.0f ns.", nanoseconds(t_start, t_stop) / iters);
  
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
  char vma_desc[80];
  sprintf(vma_desc, "VMM slopes to actuators (%dx%d)", na, ns);
  vmm = rr_new_part(recon_avg, vma_desc);
  result = rr_new(1, vmm);
  // rr_print_result(result, "VMM");
  
error:
  cfree(s);
  cfree(a);
  cfree(m);
  return result;
}

recon_result dual_vmm_reconstructor(aperture ap, int navg, int iters)
{
  double * s, *p, * a, * m, * ma;
  int nn = ap->nx * ap->ny;
  int na = NACTUATORS + NCOEFF;
  int ns = 2 * ap->ni;
  int nm = (ap->nx / 2 + 1) * ap->ny + NCOEFF;
  struct timespec t_start, t_first, t_stop;
  double recon_avg = 0.0, apply_avg = 0.0, total = 0.0;
  int i;
  
  recon_result result;
  recon_result_part vmm, vma;
  
  // Allocating plans and memory.
  clock_gettime (CLOCK_REALTIME, &t_start); 
  
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
  log_info("Warm-up average %.0f ns.", nanoseconds(t_start, t_stop) / iters);
  
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
  
  char vmm_desc[80], vma_desc[80];
  sprintf(vmm_desc, "VMM slopes to modes (%dx%d)", nm, ns);
  vmm = rr_new_part(recon_avg, vmm_desc);
  sprintf(vma_desc, "VMM modes to actuators (%dx%d)", na, nm);
  vma = rr_new_part(apply_avg, vma_desc);
  result = rr_new(2, vmm, vma);
  // rr_print_result(result, "Dual VMM");
  
error:
  cfree(s);
  cfree(p);
  cfree(a);
  cfree(m);
  cfree(ma);
  return result;
}

recon_result hybrid_reconstructor(aperture ap, int navg, int iters)
{
  double * sx, * sy, * est;
  double * a, * m;
  struct timespec t_start, t_slopemanage, t_ftr, t_stop;
  double sm_avg = 0.0, ftr_avg = 0.0, mat_avg = 0.0;
  int nn = ap->nx * ap->ny;
  int i;
  int ns = 2*ap->ni;
  int na = NACTUATORS + NCOEFF;
  int nm = (ap->nx / 2 + 1) * ap->ny + NCOEFF;
  
  recon_result result;
  recon_result_part sm, ftr, vma;
  
  
  fftw_complex *gx, *gy;
  ftr_plan plan;
  sm_plan manage_plan;
  
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
  plan = ftr_plan_reconstructor(ap->nx, ap->ny, sx, sy, est);
  check(plan, "Allocating FTR Plan");
  
  ftr_set_filter(plan, gx, gy);
  manage_plan = slope_management_plan_from_aperture(ap);
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
  log_info("Warm-up average %.0f ns.", nanoseconds(t_start, t_stop) / iters);
  
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
  
  char vma_desc[80];
  sprintf(vma_desc, "VMM modes to actuators (%dx%d)", na, nm);
  vma = rr_new_part(mat_avg, vma_desc);
  sm = rr_new_part(sm_avg, "Slope Management");
  ftr = rr_new_part(ftr_avg, "FTR");
  result = rr_new(3, sm, ftr, vma);
  
error:
  if(plan) ftr_destroy(plan);
  if(manage_plan) slope_management_destroy(manage_plan);
  ffree(sx);
  ffree(sy);
  ffree(est);
  ffree(gx);
  ffree(gy);
  cfree(a);
  cfree(m);
  return result;
  
}

#define NTESTS 4
int main (int argc, char const *argv[])
{
  // DEFAULTS
  int navg = 1e3;
  int ng = 32, nacross = 30;
  int iters = 1e5;
  int nthread_max = 1;
  
  struct recon_test {
      reconstructor func;
      recon_result result;
      char * description;
  };
  double total, duration;
  size_t i, j;
  struct recon_test tests[NTESTS] = {
    { .func = pure_ftr_reconstructor, .result = NULL, .description = "FTR" },
    { .func = pure_vmm_reconstructor, .result = NULL, .description = "VMM" },
    { .func = dual_vmm_reconstructor, .result = NULL, .description = "Dual VMM" },
    { .func = hybrid_reconstructor,   .result = NULL, .description = "FTR + VMM" },
  };
  int order[NTESTS] = { 0, 1, 2, 3};
  shuffle(order, NTESTS);
  aperture ap;
  srand((unsigned) time(NULL));
  
  
  // Parse Arguments.
  if(argc > 1) iters = (int)atof(argv[1]);
  if(argc > 2) ng = (int)atoi(argv[2]);
  if(argc > 3) nacross = (int)atoi(argv[3]);
  if(argc > 4) nthread_max = (int)atoi(argv[4]);
  
  ap = ShaneAO_aperture(ng, nacross);
  printf("Simulating %d-across on a grid of %d x %d with %d subapertures.\n", nacross, ap->nx, ap->ny, ap->ni);
  
  if(nthread_max > 1){
    printf("Testing FTR threading, from 1 to %d threads.\n", nthread_max);
    ftr_init(1);
    struct recon_test * threadtests, test;
    char desc[MAX_STR_LEN];
    
    threadtests = malloc(sizeof(struct recon_test) * nthread_max);
    for(i = 0; i < nthread_max; ++i)
    {
      sprintf(desc, "FTR with %d threads", (int)i + 1);
      threadtests[i] = (struct recon_test){ .func = pure_ftr_reconstructor, .result = NULL, .description = NULL };
      threadtests[i].description = malloc(sizeof(char) * MAX_STR_LEN);
      strncpy(threadtests[i].description, desc, MAX_STR_LEN);
    }
    for(i = 0; i < nthread_max; ++i)
    {
      printf("--> %s\n", threadtests[i].description);
      fftw_plan_with_nthreads(i + 1);
      threadtests[i].result = threadtests[i].func(ap, navg, iters);
    }
    printf("--> Threading Summary:\n");
    printf("========================================== ========= ==========\n");
    for(i = 0; i < nthread_max; ++i)
    {
      rr_print_result(threadtests[i].result, threadtests[i].description);
    }
  }
  
  
  printf("Testing %d iterations.\n", iters);
  
  for(i = 0; i < NTESTS; ++i)
  {
      printf("--> %s\n", tests[order[i]].description);
      tests[order[i]].result = tests[order[i]].func(ap, navg, iters);
  }
  
  printf("--> Summary:\n");
  printf("========================================== ========= ==========\n");
  
  for(i = 0; i < NTESTS; ++i)
  {
      rr_print_result(tests[i].result, tests[i].description);
  }
  
  return 1;
}