//
//  aperture.h
//  libFTR
//
//  Created by Alexander Rudy on 2015-10-20.
//  Copyright 2015 Alexander Rudy. All rights reserved.
//
#ifndef APERTURE_H
#define APERTURE_H

#include "dbg.h"
#include <stdlib.h>
#include <math.h>

/*
Tools for making apertures.
*/

// Structures for apertures.
typedef struct aperture_s * aperture;

typedef struct aperture_s {
  int nx, ny; // Dimensions of the aperture.
  int nm; // Number of modes.
  int * ap; // Pointer to the aperture array.
  int ni; // Number of illuminated subapertures.
};

/*
Function Prototypes
*/
aperture aperture_create_default(const int ny, const int nx);
aperture aperture_create_with_radii(const int ny, const int nx, const double outer_radius, const double inner_radius);
aperture aperture_create(const int ny, const int nx, const int *ap);
void aperture_print(const aperture ap);
void aperture_destroy(aperture ap);

/*
Functions to enable compatibility between slopemanagement
and aperture.h
*/
#ifdef SLOPEMANAGE_H
#define slope_management_plan_from_aperture(A) slope_management_plan(A->ny, A->nx, A->ap)
#endif

/*
Make a standard circular aperture with 1px padding around the edges
and no central obscuration.
Returns the aperture structure.
*/
aperture aperture_create_default(const int ny, const int nx)
{
    int n = nx > ny ? ny : nx;
    double outer_radius = ((double)n / 2.0) - 1.0;
    double inner_radius = outer_radius < 6 ? 0.0 : (outer_radius / 3.0);
    return aperture_create_with_radii(ny, nx, outer_radius, inner_radius);
}

/*
Make a standard circular aperture with an inner and outer radius provided.
Returns the number of illuminated subapertures.
*/
aperture aperture_create_with_radii(const int ny, const int nx, const double outer_radius, const double inner_radius)
{
  int i, j, *apid=NULL;
  double x, y;
  double radius;
  aperture ap=NULL;
  
  apid = malloc(sizeof(int) * nx * ny);
  check_mem(apid);
  
  for(i = 0; i < nx; ++i)
  {
      x = i - ((double)nx / 2.0) + 0.5;
      for(j = 0; j < ny; ++j)
      {
          y = j - ((double)ny / 2.0) + 0.5;
          radius = sqrt((x * x) + (y * y));
          if(radius < outer_radius && radius >= inner_radius) {
              apid[(i * nx) + j] = 1;
          }else{
              apid[(i * nx) + j] = 0;
          }
      }
  }
  ap = aperture_create(ny, nx, apid);
error:
  if(apid) free(apid);
  return ap;
}

aperture aperture_create(const int ny, const int nx, const int *ap)
{
  size_t i;
  aperture _ap = NULL;
  _ap = malloc(sizeof(struct aperture_s));
  check_mem(_ap);
  
  _ap->nx = nx;
  _ap->ny = ny;
  _ap->ni = 0;
  _ap->ap = malloc(sizeof(int) * nx * ny);
  check_mem(_ap->ap);
  memcpy(_ap->ap, ap, sizeof(int) * nx * ny);
  for(i = 0; i < (nx * ny); ++i)
  {
    check(_ap->ap[i] == 0 || _ap->ap[i] == 1, "Aperture arrays must be 1 or 0")
    _ap->ni += _ap->ap[i];
  }
  return _ap;
error:
  aperture_destroy(_ap);
  return NULL;
}

void aperture_destroy(aperture ap)
{
  if(ap){
    if(ap->ap) free(ap->ap);
    free(ap);
  }
  return;
}

void aperture_print(const aperture ap)
{
  int i, j, ns = 0;
  double x, y;
  double radius;
  double rguess = 0.0;
  for(i = 0; i < ap->nx; ++i)
  {
    x = i - ((double)ap->nx / 2.0) + 0.5;
    
      for(j = 0; j < ap->ny; ++j)
      {
        y = j - ((double)ap->ny / 2.0) + 0.5;
        radius = sqrt((x * x) + (y * y));
        printf("%d ", ap->ap[(i * ap->nx) + j]);
        ns += ap->ap[(i * ap->nx) + j];
        if(ap->ap[(i * ap->nx) + j] == 1 && radius > rguess) rguess = radius;
      }
      printf("\n");
  }
  printf("ns = %d radius? = %.1f\n", ns, rguess);
}

#endif /* end of include guard: APERTURE_H */