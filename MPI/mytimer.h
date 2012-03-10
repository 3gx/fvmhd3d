#include <sys/time.h>

#ifndef __TIMER__H__
#define __TIMER__H__

namespace mytimer {
	static double get_wtime(){
		struct timeval tv;
		gettimeofday(&tv, NULL);
		return tv.tv_sec + 1.e-6 * tv.tv_usec;
	}
};

struct Timer{
	const char *name;
	FILE       *fp;
	const char *format;
	double tstart;

	Timer(
			const char *_name,
			FILE       *_fp     = stdout,
			const char *_format = " %-10s : %f sec\n")
		: name(_name), fp(_fp), format(_format)
	{
		tstart = mytimer::get_wtime();
	}
	~Timer(){
		double tend = mytimer::get_wtime();
		fprintf(fp, format, name, tend - tstart);
		fflush(fp);
	}
};

#endif // __TIMER__H__
