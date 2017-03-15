#ifndef GET_TIMESTAMP_H
#define GET_TIMESTAMP_H

#include <sys/time.h>

typedef unsigned long long ts_t;

ts_t t0__;
ts_t t1__;
double ms__;
ts_t   ns__;
#define START_TIMING()  t0__ = get_timestamp();
#define STOP_TIMING()   t1__ = get_timestamp();\
                        ms__ = (t1__ - t0__) / 1000.0L;\
                        ns__ = (t1__ - t0__);

/*
get_timestamp
    gets current timestamp in milliseconds.

usage:
    ts_t t0 = get_timestamp();
    ts_t t1 = get_timestamp();
    double secs = (t1 - t0) / 1000000.0L;
    double ms = (t1 - t0) / 1000.0L;
*/
static ts_t get_timestamp() {
    struct timeval now;
    gettimeofday (&now, NULL);
    return  now.tv_usec + (ts_t)now.tv_sec * 1000000;
}

#endif

