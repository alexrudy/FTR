#include "aperture.h"

int main (int argc, char const *argv[])
{
  aperture ap_16x;
  ap_16x = aperture_create_with_radii(32, 32, 15.0, 4.0);
  aperture_print(ap_16x);
  printf("\n\n\n");
  ap_16x = aperture_create_with_radii(20, 20, 7.5, 2.0);
  aperture_print(ap_16x);
  printf("\n\n\n");
  ap_16x = aperture_create_with_radii(10, 10, 3.6, 1.1);
  aperture_print(ap_16x);
  return 0;
}