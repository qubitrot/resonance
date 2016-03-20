#if defined(__GNUC__) || defined(__clang__)

#include "precise_timer.h"
#include "precise_timer_impl.h"
#include <time.h>

struct TimerImpl_Linux final : TimerImpl {
	public:
		TimerImpl_Linux()
			: started(false)
			, stopped(true)
		{
		}

		virtual void start() {
			if (!started) {
				clock_gettime(CLOCK_MONOTONIC, &start_time);
				started = true;
			}
			stopped = false;
		}

		virtual void stop() {
			clock_gettime(CLOCK_MONOTONIC, &end_time);
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
					struct timespec now_time;
					clock_gettime(CLOCK_MONOTONIC, &now_time);
					return calculateNanoseconds(start_time, now_time);
				}
			} else {
				return 0;
			}
		}

	private:
		uint64_t calculateNanoseconds(struct timespec& start, struct timespec& end) {
			uint64_t seconds = (end.tv_sec - start.tv_sec) * 1000000000;
			return seconds + (end.tv_nsec - start.tv_nsec);
		}

	private:
		bool started;
		bool stopped;
		struct timespec start_time;
		struct timespec end_time;
};

Timer::Timer()
	: impl(new TimerImpl_Linux)
{
}

#endif
