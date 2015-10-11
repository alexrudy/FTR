//
//  slopemanage.h
//  FTR
//
//  Created by Alexander Rudy on 2015-10-10.
//  Copyright 2015 Alexander Rudy. All rights reserved.
//

#ifndef SLOPEMANAGE_H_18A446D7
#define SLOPEMANAGE_H_18A446D7

typedef struct slope_management_plan_s {
    int nx;
    int ny;
    int *row_any, *left, *right;
    int *col_any, *top, *bottom;
    double *y_row_sum, *x_col_sum;
    int *ap;
} * sm_plan;

sm_plan slope_management_plan(int nx, int ny, int *ap);
void slope_management_execute(sm_plan plan, double * sx, double * sy);
void slope_management(int nx, int ny, int *ap, double * sx, double * sy);
void slope_management_free(sm_plan);
#endif /* end of include guard: SLOPEMANAGE_H_18A446D7 */
