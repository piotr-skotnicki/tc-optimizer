#include "timer.h"

#include <stdlib.h>
#include <sys/time.h>

struct tc_timer
{
    struct timeval start_time;
    struct timeval end_time;
};

struct tc_timer* tc_timer_start()
{
    tc_timer* timer = (tc_timer*)malloc(sizeof(tc_timer));

    gettimeofday(&timer->start_time, NULL);

    return timer;
}

long tc_timer_stop(struct tc_timer* timer)
{
    gettimeofday(&timer->end_time, NULL);

    long elapsed = 1000L * (timer->end_time.tv_sec - timer->start_time.tv_sec)
                         + (timer->end_time.tv_usec - timer->start_time.tv_usec) / 1000L;

    free(timer);

    return elapsed;
}
