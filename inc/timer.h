
#ifndef TIMER_H
#define TIMER_H

typedef struct {
        double real, cpu;
} timer_struct;

#ifdef __cplusplus
extern "C" {
#endif
void tic(timer_struct * t);
void toc(timer_struct * t, timer_struct * stat_timer);

#ifdef __cplusplus
}
#endif
#endif // end of definition of a timer
