#pragma once
#ifndef PLAIN_WRITER_H_
#define PLAIN_WRITER_H_

#include "iwriter.h"

class PlainTextWriter : public IWriter {
		void writeRecord(std::ostream& out, Profiler::Record* node, const std::string& indentation, float percent, float percent_total) {
			out << indentation << '(' << percent << "%)  ";
			if (node->identifier)
				out << node->identifier << '[' << node->func_name << ']';
			else
				out << node->func_name;
			out << "  ran #" << node->total_count << " times";
			if (node->max_recursion_depth > 0) {
				out << "  recursion lvl: " << node->max_recursion_depth+1;
			}
			out << "  avg ns: " << node->avg_ns << std::endl;
		}

		void writeHeader(std::ostream&) {
		}

		void writeFooter(std::ostream&) {
		}
};

#endif
