#include "precise_timer.h"
#include "precise_timer_impl.h"

Timer::~Timer() {
	delete impl;
}

Timer& Timer::start() {
	impl->start();
	return *this;
}

Timer& Timer::stop() {
	impl->stop();
	return *this;
}

Timer& Timer::restart() {
	impl->reset();
	impl->start();
	return *this;
}

Timer& Timer::reset() {
	impl->reset();
	return *this;
}

uint64_t Timer::nanoseconds() {
	return impl->nanoseconds();
}
