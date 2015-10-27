//
//  ShaneAO_apertures.h
//  libFTR
//
//  Created by Alexander Rudy on 2015-10-27.
//  Copyright 2015 Alexander Rudy. All rights reserved.
//
#ifndef SHANEAO_APERTURES_H_437C7A2C
#define SHANEAO_APERTURES_H_437C7A2C

#include "aperture.h"

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
  aperture ap;
  switch (nacross) {
  case 8:
    outer = 3.6;
    break;
  default:
    outer = ((double)nacross / 2.0) - 1.0;
  }
  inner = 0.3 * outer;
  ap = aperture_create_with_radii(ng, ng, outer, inner);
  switch (nacross) {
  case 8:
    ap->nm = 60;
    break;
  case 16:
    ap->nm = 224;
    break;
  case 30:
    ap->nm = 840;
  default:
    ap->nm = (ap->nx - 1) * (ap->ny - 1);
  }
  
  return ap;
}


#endif /* end of include guard: SHANEAO_APERTURES_H_437C7A2C */
