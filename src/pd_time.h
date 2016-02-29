#ifndef PD_TIME_H
#define PD_TIME_H


#include <time.h>


double
pd_time_diff_ms(struct timespec const *start, struct timespec const *end)
{
        return ((end->tv_sec*1000000000L + end->tv_nsec) -
                (start->tv_sec*1000000000L + start->tv_nsec))*0.000001;
}


#endif
