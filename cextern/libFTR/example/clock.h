#ifndef CLOCK_C_5CCA3F39
#define CLOCK_C_5CCA3F39

#include <time.h>               /* struct timespec */

#ifdef __MACH__
int clock_gettime(int clk_id, struct timespec *t);
#define CLOCK_REALTIME 0
#define CLOCK_MONOTONIC 0
#endif

#ifndef CLOCK_MONOTONIC_RAW     /* Not available on CentOS 5, kernel < 2.6.28 */
#define CLOCK_MONOTONIC_RAW CLOCK_MONOTONIC
#endif

double moving_average(double average, double duration, int navg);
double nanoseconds(const struct timespec start, const struct timespec end);

#endif /* end of include guard: CLOCK_C_5CCA3F39 */