#pragma once
#ifndef PROFILER_H_
#define PROFILER_H_

#include <iostream>
#include <string>
#include <cstdint>
//#include "object_pool/objectpool.h"

class IWriter;

class Profiler {
	friend class ProfileTimer;

	public:
		struct Record;

		static void setLogFile(const std::string& filename);
		static void setLogFormat(IWriter* writer);

	private:
		Profiler();
		~Profiler();

		static Profiler& get();

		void openRecord(const char* func_name);
		void closeRecord(const char* func_name, const char* ident, uint64_t ns);

	public:
		struct Record {
			Record()
				: recursion(0)
				, max_recursion_depth(0)
				, total_count(0)
				, total_overflowed(false)
				, total_ns(0)
				, avg_ns(0)
				, func_name(nullptr)
				, identifier(nullptr)
				, parent(nullptr)
				, children(nullptr)
				, siblings(nullptr)
			{}
			uint32_t    recursion;
			uint32_t    max_recursion_depth;
			uint64_t    total_count      : 63;
			bool        total_overflowed : 1;
			uint64_t    total_ns;
			uint64_t    avg_ns;
			const char* func_name;
			const char* identifier;
			// linked lists
			Record*     parent;
			Record*     children;
			Record*     siblings;
		};

	private:
		       void    printRecord(Record* node, const std::string& indentation);
		       void    printRecordTree(Record* node, std::string& indentation);
		static Record* findChild(Record* node, const char* func_name);
		static void    addChild(Record* node, Record* child);

	private:
		Record* root_;
		Record* current_;

		IWriter*      writer_;
		std::ostream* out;
		bool          should_del_stream_;

		//ObjectPool<Record> objpool;
};

#include "currentfunc.h"
#include "profile_timer.h"

#ifndef NDEBUG
#define STRINGIZE(x) STRINGIZE_IMPL(x)
#define STRINGIZE_IMPL(x) #x

#define PROFILE() ProfileTimer profile_timer(CURRENT_FUNCTION) //, ##__VA_ARGS__)
#define PROFILE_EXTERNAL(func, ...) [&]() { \
	ProfileTimer profile_external_timer(#func "__" __FILE__ ":" STRINGIZE(__LINE__)); \
	return func(__VA_ARGS__); \
}()
#else
#define PROFILE()
#define PROFILE_EXTERNAL(func, ...) func(__VA_ARGS__)
#endif

#endif
