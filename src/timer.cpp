
#include <stdio.h>
#include <sys/resource.h>
#include <sys/time.h>

#include "timer.h"


double cputime()
{
        struct rusage r;
        getrusage(RUSAGE_SELF, &r);
        return r.ru_utime.tv_sec + r.ru_stime.tv_sec + 1e-6 * (r.ru_utime.tv_usec + r.ru_stime.tv_usec);
}

double realtime()
{
        struct timeval tp;
        struct timezone tzp;
        gettimeofday(&tp, &tzp);
        return tp.tv_sec + tp.tv_usec * 1e-6;
}

void tic(timer_struct * t){
        t->cpu  = cputime();
        t->real = realtime();
}

void toc(timer_struct * t, timer_struct * stat_timer){
        stat_timer->cpu  += cputime() - t->cpu;
        stat_timer->real += realtime()  - t->real;
}
