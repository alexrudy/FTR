#include "ftr.h"
#include "slopemanage.h"
#include "ftr_example.h"
#include <cblas.h>


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

void init(double * pt, int nn, double value)
{
    int i;
    for(i = 0; i < nn; ++i)
    {   
        pt[i] = value;
    }
}

double pure_ftr_reconstructor(int nx, int ny, int navg, int iters)
{
  double * sx, * sy, *est;
  struct timespec t_start, t_slopemanage, t_ftr, t_stop;
  double sm_avg, ftr_avg;
  int nn = nx * ny;
  int i;
  int *ap;
  
  fftw_complex *gx, *gy;
  ftr_plan plan;
  sm_plan manage_plan;
  
  printf("Conducting %d test iterations of the pure FTR.\n", iters);
  
  // Allocating plans and memory.
  printf("Allocating...");
  clock_gettime (CLOCK_REALTIME, &t_start); 
  
  sx = malloc(sizeof(double) * nn);
  sy = malloc(sizeof(double) * nn);
  est = malloc(sizeof(double) * nn);
  ap = malloc(sizeof(int) * nn);
  
  gx = fftw_malloc(sizeof(fftw_complex) * nn);
  gy = fftw_malloc(sizeof(fftw_complex) * nn);
  
  clock_gettime(CLOCK_REALTIME, &t_stop);
  printf(" %.0f ns", nanoseconds(t_start, t_stop));
  printf("\n");
  
  make_aperture(ap, nx, ny);
  
  // Zeroing arrays. We don't time this becasue it isn't really relevant.
  for(i = 0; i < nn; ++i)
  {
      gx[i] = 1.0;
      gy[i] = 1.0;
      sx[i] = 1.0;
      sy[i] = 1.0;
      est[i] = 0.0;
  }
  
  printf("Planning...");
  clock_gettime (CLOCK_REALTIME, &t_start); 
  
  plan = ftr_plan_reconstructor(nx, ny, sx, sy, est);
  ftr_set_filter(plan, gx, gy);
  manage_plan = slope_management_plan(ny, nx, ap);
  clock_gettime(CLOCK_REALTIME, &t_stop);
  printf(" %.0f ns", nanoseconds(t_start, t_stop));
  printf("\n");
  
  // Re-zero arrays which might have been touched by the planners.
  for(i = 0; i < nn; ++i)
  {
      sx[i] = 1.0;
      sy[i] = 1.0;
      est[i] = 0.0;
  }
  
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
        printf("First iteration complete...\n");
      }
  }
  
  printf("Done timing.\n");
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
  
  
  printf("Conducting %d test iterations of the pure VMM.\n", iters);
  
  // Allocating plans and memory.
  printf("Allocating...");
  clock_gettime (CLOCK_REALTIME, &t_start); 
  
  s = malloc(sizeof(double) * ns);
  a = malloc(sizeof(double) * na);
  m = malloc(sizeof(double) * na * ns);
  
  clock_gettime(CLOCK_REALTIME, &t_stop);
  printf(" %.0f ns", nanoseconds(t_start, t_stop));
  printf("\n");
  
  init(s, ns, 1.0);
  init(a, na, 0.0);
  init(m, na * ns, 1.0);
  
  
  for(i = 0; i < iters; ++i)
  {
      // This part you might do many many times in a for-loop.
      if(i == 0)
      {
          printf("Reconstructing...\n");
      }
      clock_gettime (CLOCK_REALTIME, &t_start);
      cblas_dgemv(CblasRowMajor, CblasNoTrans, na, ns, 1.0, m, ns, s, 1, 0.0, a, 1);
      clock_gettime (CLOCK_REALTIME, &t_stop);
      if(i == 0)
      {
          recon_avg = nanoseconds(t_start, t_stop);
      }
      recon_avg = moving_average(recon_avg, nanoseconds(t_start, t_stop), navg);
      if(i == 0)
      {
        printf("First iteration complete...\n");
      }
  }
  printf("Done timing.\n");
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
  struct timespec t_start, t_stop;
  double recon_avg;
  int i;
  
  
  printf("Conducting %d test iterations of the dual VMM.\n", iters);
  
  // Allocating plans and memory.
  printf("Allocating...");
  clock_gettime (CLOCK_REALTIME, &t_start); 
  
  s = malloc(sizeof(double) * ns);
  p = malloc(sizeof(double) * na);
  a = malloc(sizeof(double) * na);
  m = malloc(sizeof(double) * na * ns);
  ma = malloc(sizeof(double) * na * na);
  
  clock_gettime(CLOCK_REALTIME, &t_stop);
  printf(" %.0f ns", nanoseconds(t_start, t_stop));
  printf("\n");
  
  init(s, ns, 1.0);
  init(a, na, 0.0);
  init(p, na, 0.0);
  init(m, na * ns, 1.0);
  init(ma, na * na, 1.0);
  
  for(i = 0; i < iters; ++i)
  {
      // This part you might do many many times in a for-loop.
      if(i == 0)
      {
          printf("Reconstructing...\n");
      }
      clock_gettime (CLOCK_REALTIME, &t_start);
      cblas_dgemv(CblasRowMajor, CblasNoTrans, na, ns, 1.0, m, ns, s, 1, 0.0, p, 1);
      cblas_dgemv(CblasRowMajor, CblasNoTrans, na, na, 1.0, ma, na, p, 1, 0.0, a, 1);
      clock_gettime (CLOCK_REALTIME, &t_stop);
      if(i == 0)
      {
          recon_avg = nanoseconds(t_start, t_stop);
      }
      recon_avg = moving_average(recon_avg, nanoseconds(t_start, t_stop), navg);
      if(i == 0)
      {
        printf("First iteration complete...\n");
      }
  }
  printf("Done timing.\n");
  printf("Averaged %.0f ns for matrix multiply.\n", recon_avg);
  printf("Total iteration time %.0f ns corresponds to a rate of %.0f kHz\n", recon_avg, 1e6 / recon_avg);
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
  
  fftw_complex *gx, *gy;
  ftr_plan plan;
  sm_plan manage_plan;
  
  printf("Conducting %d test iterations of the FTR + VMM.\n", iters);
  
  // Allocating plans and memory.
  printf("Allocating...");
  clock_gettime (CLOCK_REALTIME, &t_start); 
  
  sx = malloc(sizeof(double) * nn);
  sy = malloc(sizeof(double) * nn);
  est = malloc(sizeof(double) * nn);
  ap = malloc(sizeof(int) * nn);
  a = malloc(sizeof(double) * nn);
  m = malloc(sizeof(double) * nn * nn);
  
  gx = fftw_malloc(sizeof(fftw_complex) * nn);
  gy = fftw_malloc(sizeof(fftw_complex) * nn);
  
  clock_gettime(CLOCK_REALTIME, &t_stop);
  printf(" %.0f ns", nanoseconds(t_start, t_stop));
  printf("\n");
  
  make_aperture(ap, nx, ny);
  
  init(sx, nn, 1.0);
  init(sy, nn, 1.0);
  init(m, nn * nn, 1.0);
  init(a, nn, 0.0);
  init(est, nn, 0.0);
  
  // Zeroing arrays. We don't time this becasue it isn't really relevant.
  for(i = 0; i < nn; ++i)
  {
      gx[i] = 1.0;
      gy[i] = 1.0;
  }
  
  printf("Planning...");
  clock_gettime (CLOCK_REALTIME, &t_start); 
  
  plan = ftr_plan_reconstructor(nx, ny, sx, sy, est);
  ftr_set_filter(plan, gx, gy);
  manage_plan = slope_management_plan(ny, nx, ap);
  clock_gettime(CLOCK_REALTIME, &t_stop);
  printf(" %.0f ns", nanoseconds(t_start, t_stop));
  printf("\n");
  
  // Re-zero arrays which might have been touched by the planners.
  for(i = 0; i < nn; ++i)
  {
      sx[i] = 1.0;
      sy[i] = 1.0;
      est[i] = 0.0;
  }
  
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
        printf("Matrix Multiply...\n");
      }
      cblas_dgemv(CblasRowMajor, CblasNoTrans, nn, nn, 1.0, m, nn, est, 1, 0.0, a, 1);
      clock_gettime (CLOCK_REALTIME, &t_stop);
      
      if(i == 0)
      {
          sm_avg = nanoseconds(t_start, t_slopemanage);
          ftr_avg = nanoseconds(t_slopemanage, t_ftr);
          mat_avg = nanoseconds(t_ftr, t_stop);
          printf("Timing...\n");
      }
      sm_avg = moving_average(sm_avg, nanoseconds(t_start, t_slopemanage), navg);
      ftr_avg = moving_average(ftr_avg, nanoseconds(t_slopemanage, t_ftr), navg);
      mat_avg = moving_average(mat_avg, nanoseconds(t_ftr, t_stop), navg);
      if(i == 0)
      {
        printf("First iteration complete...\n");
      }
  }
  
  printf("Done timing.\n");
  printf("Averaged %.0f ns for slope management.\n", sm_avg);
  printf("Averaged %.0f ns for ftr.\n", ftr_avg);
  printf("Averaged %.0f ns for vmm.\n", mat_avg);
  printf("Total iteration time %.0f ns corresponds to a rate of %.0f kHz\n", (sm_avg + ftr_avg + mat_avg), 1e6 / (sm_avg + ftr_avg + mat_avg));
  return (sm_avg + ftr_avg + mat_avg);
}

int main (int argc, char const *argv[])
{
  int navg = 100;
  double ftr, vmm, hybrid, dual;
  int nx = 10;
  int ny = 10;
  int iters;
  printf("Parsing\n");
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
  printf("Running %d iterations\n", iters);
  
  ftr = pure_ftr_reconstructor(nx, ny, navg, iters);
  vmm = pure_vmm_reconstructor(nx, ny, navg, iters);
  hybrid = hybrid_reconstructor(nx, ny, navg, iters);
  dual = dual_vmm_reconstructor(nx, ny, navg, iters);
  
  printf("Summary:\n");
  printf("FTR:      %6.0f ns %.1f kHz\n", ftr, 1e6 / ftr);
  printf("VMM:      %6.0f ns %.1f kHz\n", vmm, 1e6 / vmm);
  printf("Hybrid:   %6.0f ns %.1f kHz\n", hybrid, 1e6 / hybrid);
  printf("Dual VMM: %6.0f ns %.1f kHz\n", dual, 1e6 / dual);
  
  return 1;
}