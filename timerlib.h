#ifndef __TIMERLIB__
#define __TIMERLIB__

void init_timer();
void accum_timer(double val, int index);
void print_timers(double n, double nb);
void print_timers_mpi();
void get_cputime(double * laptime, double * splittime);
void print_current_time(char * message);
void init_current_time();

//#define rdtscl(low)						
//	__asm__ __volatile__("rdtsc" : "=a" (low) : : "edx")

inline static void rdtscl(unsigned long * count)
{
    unsigned int low,high;
    __asm__ __volatile__("rdtsc" : "=a" (low), "=d" (high));
    *count =  low + (((long) high)<<32);
}


#ifdef TIMETEST
#define BEGIN_TSC unsigned long rdtsc1;	rdtscl(&rdtsc1)
#define END_TSC(cycle,index) unsigned long rdtsc2;	rdtscl(&rdtsc2);	\
    double cycle =(double) rdtsc2-(double)rdtsc1;			\
    accum_timer(cycle,index)
#define BEGIN_TIMER(t)  unsigned long t;	rdtscl(&t)
#define END_TIMER(t,index,ops) {unsigned long rdtsc2; rdtscl(&rdtsc2);	\
	accum_timer_and_ops((double) rdtsc2-(double)(t),(index),(ops));}
#else
#define BEGIN_TSC 
#define END_TSC(cycle,index)  
#define BEGIN_TIMER
#define END_TIMER(cycle,index)  
#endif
#endif
