//
//  aperture.h
//  libFTR
//
//  Created by Alexander Rudy on 2015-10-20.
//  Copyright 2015 Alexander Rudy. All rights reserved.
//
#ifndef APERTURE_H_E17EC070
#define APERTURE_H_E17EC070

/*
Tools for making apertures.
*/
int make_aperture(const int ny, const int nx, int * ap);
int make_aperture_with_radii(const int ny, const int nx, int * ap, const double outer_radius, const double inner_radius);
void print_aperture(const int ny, const int nx, int * ap);

/*
Make a standard circular aperture with 1px padding around the edges
and no central obscuration.
Returns the number of illuminated subapertures.
*/
int make_aperture(const int ny, const int nx, int * ap)
{
    double outer_radius = ((double)nx / 2.0) - 1.0;
    double inner_radius = outer_radius / 6.0 + 0.5;
    return make_aperture_with_radii(ny, nx, ap, outer_radius, inner_radius);
}

/*
Make a standard circular aperture with an inner and outer radius provided.
Returns the number of illuminated subapertures.
*/
int make_aperture_with_radii(const int ny, const int nx, int * ap, const double outer_radius, const double inner_radius)
{
  int i, j, t = 0;
  double x, y;
  double radius;
  for(i = 0; i < nx; ++i)
  {
      x = i - ((double)nx / 2.0) + 0.5;
      for(j = 0; j < ny; ++j)
      {
          y = j - ((double)ny / 2.0) + 0.5;
          radius = sqrt((x * x) + (y * y));
          if(radius < outer_radius && radius >= inner_radius) {
              ap[(i * nx) + j] = 1;
          }else{
              ap[(i * nx) + j] = 0;
          }
          t += ap[(i * nx) + j];
      }
  }
  return t;
}

void print_aperture(const int ny, const int nx, int * ap)
{
  int i, j, ns = 0;
  double x, y;
  double radius;
  double rguess = 0.0;
  for(i = 0; i < nx; ++i)
  {
    x = i - ((double)nx / 2.0) + 0.5;
    
      for(j = 0; j < ny; ++j)
      {
        y = j - ((double)ny / 2.0) + 0.5;
        radius = sqrt((x * x) + (y * y));
        printf("%d ", ap[(i * nx) + j]);
        ns += ap[(i * nx) + j];
        if(ap[(i * nx) + j] == 1 && radius > rguess) rguess = radius;
      }
      printf("\n");
  }
  printf("ns = %d radius? = %.1f\n", ns, rguess);
}

#endif /* end of include guard: APERTURE_H_E17EC070 */