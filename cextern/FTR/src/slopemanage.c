//
//  slopemanage.c
//  FTR
//
//  Created by Alexander Rudy on 2015-10-10.
//  Copyright 2015 Alexander Rudy. All rights reserved.
//

#include "slopemanage.h"
#include <stdio.h>

struct slope_management_plan_s {
    int nx, ny, nn;
    int *row_any, *left, *right;
    int *col_any, *top, *bottom;
    double *y_row_sum, *x_col_sum;
    int *ap;
};

//TODO: This function should provide error checking for invalid apertures.
sm_plan slope_management_plan(int ny, int nx, int *ap)
{
    int r, c;
    int *row, *column;
    int cell;
    
    sm_plan plan;
    plan = malloc(sizeof(struct slope_management_plan_s));
    
    plan->ap = ap;
    plan->nx = nx;
    plan->ny = ny;
    plan->nn = plan->nx > plan->ny ? plan->nx : plan->ny;
    
    /* Figure out what is going on with the rows */
    plan->row_any = malloc(sizeof(int) * ny);
    plan->left = malloc(sizeof(int) * ny);
    plan->right = malloc(sizeof(int) * ny);
    
    for(r = 0; r < ny; ++r)
    {
        row = (ap + (r * nx));
        // Initialize these values so that we can see if we've already set them elsewhere.
        plan->left[r] = -1;
        plan->right[r] = -1;
        plan->row_any[r] = 0;
        for(c = 0; c < nx; ++c)
        {
            cell = *(row + c);
            if(cell != 0 && plan->left[r] == -1)
            {
                plan->left[r] = c;
            }
            if(plan->left[r] != -1 && cell != 0)
            {
                plan->right[r] = c;
            }
        }
        if(plan->left[r] != -1 && plan->right[r] != -1)
        {
            plan->row_any[r] = 1;
        }
    }
    
    plan->col_any = malloc(sizeof(int) * nx);
    plan->top = malloc(sizeof(int) * nx);
    plan->bottom = malloc(sizeof(int) * nx);
    
    for(c = 0; c < nx; ++c)
    {
        column = (ap + c);
        // Initialize these values so that we can see if we've already set them elsewhere.
        plan->top[c] = -1;
        plan->bottom[c] = -1;
        plan->col_any[c] = 0;
        
        for(r = 1; r < ny; ++r)
        {
            cell = *(column + (r * nx));
            if(cell != 0 && plan->top[c] == -1)
            {
                plan->top[c] = r;
            }
            if(plan->top[c] != -1 && cell != 0)
            {
                plan->bottom[c] = r;
            }
        }
        if(plan->top[c] != -1 && plan->bottom[c] != -1)
        {
            plan->col_any[c] = 1;
        }
    }
    
    return plan;
}

void slope_management_execute(sm_plan plan, double * sy, double * sx)
{
    int i, j;
    int n;
    double row_sum, col_sum;
    double *row, *col;
    
    for(i = 0; i < plan->nn; ++i)
    {
        if(i < plan->ny && plan->row_any[i] == 1)
        {
            // Sum the row, so we can use it to edge correct.
            row_sum = 0.0;
            for(j = 0; j < plan->nx; ++j)
            {
                n = (i * plan->nx) + j;
                if(plan->ap[n] == 1)
                {
                  row_sum += sy[n];
                }
            }
            row = (sy + (i * plan->nx));
            row[plan->left[i] - 1] = -0.5 * row_sum;
            row[plan->right[i] + 1] = -0.5 * row_sum;
        }
        
        if(i < plan->nx && plan->col_any[i] == 1)
        {
            // Sum the column, so we can use it to edge correct.
            col_sum = 0.0;
            for(j = 0; j < plan->ny; ++j)
            {
                n = (j * plan->nx) + i;
                
                if(plan->ap[n] == 1)
                {
                  col_sum += sx[n];
                }
            }
            j = (plan->top[i] - 1) * plan->nx;
            col = sx;
            col[i + j] = -0.5 * col_sum;
            j = (plan->bottom[i] + 1) * plan->nx;
            col[i + j] = -0.5 * col_sum;
        }
        
    }
    return;
}

// This method just combines the two methods above if you don't care
// about doing the memory allocation every time.
void slope_management(int ny, int nx, int *ap, double * sy, double * sx)
{
    sm_plan plan;
    plan = slope_management_plan(ny, nx, ap);
    slope_management_execute(plan, sy, sx);
    slope_management_destroy(plan);
}

void slope_management_destroy(sm_plan plan)
{
    // Free row variables.
    free(plan->row_any);
    free(plan->left);
    free(plan->right);
    
    // Free column variables.
    free(plan->col_any);
    free(plan->top);
    free(plan->bottom);
    
    // Free the plan itself.
    free(plan);
    
    return;
}
