#pragma once
#ifndef PROFILE_TIMER_H_
#define PROFILE_TIMER_H_

#include "precise_timer.h"

class ProfileTimer {
	public:
		ProfileTimer(const char* func_name);
		ProfileTimer(const char* func_name, const char* ident);
		~ProfileTimer();

	private:
		void start();

	private:
		const char* identifier;
		const char* function_name;
		Timer timer;
};

#endif
