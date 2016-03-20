#pragma once
#ifndef PRECISE_TIMER_IMPL_H_
#define PRECISE_TIMER_IMPL_H_

#include <cstdint>

struct TimerImpl {
	virtual ~TimerImpl() {}

	virtual void start() = 0;
	virtual void stop() = 0;
	virtual void reset() = 0;

	virtual uint64_t nanoseconds() = 0;
};

#endif
