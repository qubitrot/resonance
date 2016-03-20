#include "profiler.h"
#include "profile_timer.h"

ProfileTimer::ProfileTimer(const char* func_name)
	: identifier(nullptr)
	, function_name(func_name)
{
	start();
}

ProfileTimer::ProfileTimer(const char* func_name, const char* ident)
	: identifier(ident)
	, function_name(func_name)
{
	start();
}

void ProfileTimer::start() {
	Profiler::get().openRecord(function_name);
	timer.start();
}

ProfileTimer::~ProfileTimer() {
	timer.stop();

	Profiler::get().closeRecord(function_name, identifier, timer.nanoseconds());
}
