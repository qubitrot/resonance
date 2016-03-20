#if defined(_WIN32) || defined(__WIN32__) || defined(__WINDOWS__)

#include "precise_timer.h"
#include "precise_timer_impl.h"
#include <windows.h>

struct TimerImpl_Win32 : TimerImpl {
	public:
		TimerImpl_Win32()
			: started(false)
			, stopped(true)
		{
			QueryPerformanceFrequency(&li);
			freq = double(li.QuadPart) / 1000000000.0;
		}

		virtual void start() {
			if (!started) {
				QueryPerformanceCounter(&li);
				start_time = li.QuadPart;
				started = true;
			}
			stopped = false;
		}

		virtual void stop() {
			QueryPerformanceCounter(&li);
			end_time = li.QuadPart;
			stopped = true;
		}

		virtual void reset() {
			started = false;
			stopped = true;
		}

		virtual uint64_t nanoseconds() {
			if (started) {
				if (stopped) {
					return calculateNanoseconds(start_time, end_time);
				} else {
					QueryPerformanceCounter(&li);
					return calculateNanoseconds(start_time, li.QuadPart);
				}
			} else {
				return 0;
			}
		}

	private:
		uint64_t calculateNanoseconds(__int64 start, __int64 end) {
			return uint64_t(double(end - start)/freq);
		}

	private:
		bool started;
		bool stopped;
		double freq;
		__int64 start_time;
		__int64 end_time;
		LARGE_INTEGER li;
};

Timer::Timer()
	: impl(new TimerImpl_Win32)
{
}

#endif
