#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <crts.h>
#include "swarg.h"
#include <time.h>
#include <sys/time.h>

extern SLAVE_FUN(sw_site_propensity)(struct _sw_site_propensity_arg *);
extern SLAVE_FUN(sw_compute_vac)(struct _swarg *);

static inline unsigned long rpcc()
{
        unsigned long time;
        asm("rtc %0"
            : "=r"(time)
            :);
        return time;
}

#define CLOCKRATE 2.0E12

void kernel_site_propensity(struct _sw_site_propensity_arg *arg)
{
        CRTS_athread_spawn(sw_site_propensity, arg);
}

void kernel_func(struct _swarg *arg)
{
        CRTS_athread_spawn(sw_compute_vac, arg);
}

void c_call_athread_join()
{
        CRTS_athread_join();
}

void call_CRTS_sync_master_array()
{
        CRTS_sync_master_array();
}