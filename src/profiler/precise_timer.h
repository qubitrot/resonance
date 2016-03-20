#pragma once
#ifndef PRECISE_TIMER_H_
#define PRECISE_TIMER_H_

#include <cstdint>

struct TimerImpl;

class Timer {
	public:
		Timer();
		~Timer();

		Timer& start();
		Timer& stop();
		Timer& restart();
		Timer& reset();

		uint64_t nanoseconds();

	private:
		TimerImpl* impl;
};

#endif
