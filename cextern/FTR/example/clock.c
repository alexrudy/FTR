#ifdef __MACH__
#include <time.h>
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
