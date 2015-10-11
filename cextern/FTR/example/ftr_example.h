//
//  ftr_example.h
//  FTR
//
//  Created by Alexander Rudy on 2015-10-10.
//  Copyright 2015 Alexander Rudy. All rights reserved.
//

#ifndef FTR_EXAMPLE_H_9FF1C8EE
#define FTR_EXAMPLE_H_9FF1C8EE

#include <time.h>               /* struct timespec */

#ifdef __MACH__
int clock_gettime(int clk_id, struct timespec *t);
#define CLOCK_REALTIME 0
#define CLOCK_MONOTONIC 0
#else
#include <time.h>
#endif

#ifndef CLOCK_MONOTONIC_RAW     /* Not available on CentOS 5, kernel < 2.6.28 */
#define CLOCK_MONOTONIC_RAW CLOCK_MONOTONIC
#endif

#include <stdlib.h>
#include "clock.c"

#endif /* end of include guard: FTR_EXAMPLE_H_9FF1C8EE */
