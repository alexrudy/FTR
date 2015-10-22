//
//  slopemanage.h
//  FTR
//
//  Created by Alexander Rudy on 2015-10-10.
//  Copyright 2015 Alexander Rudy. All rights reserved.
//
#include <stdlib.h>

#ifndef SLOPEMANAGE_H
#define SLOPEMANAGE_H

// Type definitions

// sm_plan is an opaque structure pointer used to memorize various
// aspects of plan information.
typedef struct slope_management_plan_s * sm_plan;

// Functions

// Create a slope managmenet plan for a given aperture.
sm_plan slope_management_plan(const int ny, const int nx, const int *ap);

// Execute a slope management plan on slope arrays.
void slope_management_execute(sm_plan plan, double * sy, double * sx);

// Create and exectue a plan in a single function, for when memory allocation
// isn't a limiting step.
void slope_management(const int ny, const int nx, const int *ap, double * sy, double * sx);

// Destroy the memory allocated by a slope management plan.
void slope_management_destroy(sm_plan);
#endif /* end of include guard: SLOPEMANAGE_H */
