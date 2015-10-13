//
//  slopemanage.h
//  FTR
//
//  Created by Alexander Rudy on 2015-10-10.
//  Copyright 2015 Alexander Rudy. All rights reserved.
//
#include <stdlib.h>

#ifndef SLOPEMANAGE_H_18A446D7
#define SLOPEMANAGE_H_18A446D7

typedef struct slope_management_plan_s * sm_plan;

sm_plan slope_management_plan(int nx, int ny, int *ap);
void slope_management_execute(sm_plan plan, double * sx, double * sy);
void slope_management(int nx, int ny, int *ap, double * sx, double * sy);
void slope_management_destroy(sm_plan);
#endif /* end of include guard: SLOPEMANAGE_H_18A446D7 */
