#pragma once
#ifndef IWRITER_H_
#define IWRITER_H_

#include <ostream>
#include <string>
#include <locale>
#include "profiler.h"

class IWriter {
	public:
		virtual ~IWriter() {};

		virtual void writeRecord(std::ostream& out, Profiler::Record* node, const std::string& indentation, float percent, float percent_total) =0;

		// Called before any records are written
		virtual void writeHeader(std::ostream& out) =0;
		// Called after all records have been written
		virtual void writeFooter(std::ostream& out) =0;

		struct comma_out : std::numpunct<char> {
			char do_thousands_sep()   const { return ',';  }
			std::string do_grouping() const { return "\3"; }
		};
};

#endif
