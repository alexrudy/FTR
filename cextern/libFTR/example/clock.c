//
//  clock.c
//  FTR
//
//  Created by Alexander Rudy on 2015-10-20.
//  Copyright 2015 Alexander Rudy. All rights reserved.
//

#include "clock.h"

#ifdef __MACH__
#include <mach/mach_time.h>
#define CLOCK_REALTIME 0
#define CLOCK_MONOTONIC 0
int clock_gettime(int clk_id, struct timespec *t){
    mach_timebase_info_data_t timebase;
    mach_timebase_info(&timebase);
    uint64_t time;
    time = mach_absolute_time();
    double nseconds = ((double)time * (double)timebase.numer)/((double)timebase.denom);
    double seconds = ((double)time * (double)timebase.numer)/((double)timebase.denom * 1.0e9);
    t->tv_sec = seconds;
    t->tv_nsec = nseconds - (t->tv_sec * 1.0e9);
    return 0;
}

#else
#include <stdio.h>

int mywait(int ms) {
        double n,m,k;
        //n = 306000; // timed on rtc to take 1 ms
	m = (double)ms;
	while (m>0) {
	    m = m - 1.0;
	    n = 0.5e6;
	    while (n>0) {
		n = n - 1.0;
		k = n*m;
	  }
        }
	return (int)n;
}

#endif


double nanoseconds(const struct timespec start, const struct timespec end)
{
    double total;
    total = (end.tv_sec - start.tv_sec) * 1e9;
    total += (end.tv_nsec - start.tv_nsec);
    return total;
}

double moving_average(double average, const double duration, const int navg)
{
    return ((average * (double)(navg - 1)) + duration) / (double)navg;
}
