#ifndef TC_TIMER_H
#define TC_TIMER_H

struct tc_timer;

struct tc_timer* tc_timer_start();

long tc_timer_stop(struct tc_timer* timer);

#endif // TC_TIMER_H
