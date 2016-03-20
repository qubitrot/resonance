#include "profiler.h"
#include <cassert>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "iwriter.h"
#include "plainwriter.h"

void Profiler::setLogFile(const std::string& filename) {
	std::ofstream* file = new std::ofstream;
	file->open(filename);
	if (file->is_open()) {
		Profiler::get().out = file;
		Profiler::get().should_del_stream_ = true;
	} else {
		delete file;
		std::cerr << "Failed to open file for write: " << filename << std::endl;
	}
}

void Profiler::setLogFormat(IWriter* writer) {
	delete Profiler::get().writer_;
	Profiler::get().writer_ = writer;
}

Profiler& Profiler::get() {
	static Profiler mgr;
	return mgr;
}

Profiler::Profiler()
	: root_(nullptr)
	, current_(nullptr)
	, writer_(new PlainTextWriter)
	, out(&std::cout)
	, should_del_stream_(false)
{
	root_ = new Record;
	root_->func_name = "";
	root_->identifier = "ROOT";
	current_ = root_;
}

Profiler::~Profiler() {
	if (root_) {
		// calculate root total time
		Record* child = root_->children;
		while(child) {
			root_->total_ns += child->total_ns;
			root_->avg_ns += child->avg_ns;
			child = child->siblings;
		}

		writer_->writeHeader(*out);

		*out << std::setprecision(3); // for the percentages
		std::string indent;
		printRecordTree(root_, indent);
		delete root_;

		writer_->writeFooter(*out);
	}

	delete writer_;
	if (should_del_stream_)
		delete out;
}

void Profiler::printRecord(Record* node, const std::string& indentation) {
	float percent = (node->parent) ? 100.f * (float)node->total_ns / node->parent->total_ns : 100.f;
	float percent_total = 100.f * (float)node->total_ns / root_->total_ns;
	writer_->writeRecord(*out, node, indentation, percent, percent_total);
}

void Profiler::printRecordTree(Record* node, std::string& indentation) {
	printRecord(node, indentation);

	Record* last_child;
	Record* child = node->children;
	uint64_t unprofiled_ns = node->total_ns;
	while (child) {
		indentation.push_back('\t');
		printRecordTree(child, indentation);
		indentation.pop_back();

		unprofiled_ns -= child->total_ns;

		last_child = child;
		child = child->siblings;
		delete last_child;
	}

	if (node->children) {
		Record tmp = *node;
		tmp.parent = node;
		indentation.push_back('\t');
		if (unprofiled_ns > 0) {
			tmp.func_name = "Unprofiled";
			tmp.identifier = nullptr;
			tmp.total_ns = unprofiled_ns;
			tmp.avg_ns = tmp.total_ns / tmp.total_count;
			tmp.children = nullptr;
			printRecord(&tmp, indentation);
		}

		tmp.func_name = "Total";
		tmp.total_ns = node->total_ns;
		//printRecord(&tmp, indentation);
		indentation.pop_back();
	}
}

static uint32_t max_uint32(uint32_t lhs, uint32_t rhs) {
	if (rhs > lhs) return rhs;
	else           return lhs;
}

void Profiler::openRecord(const char* func_name) {
	if (func_name == current_->func_name) {
		current_->max_recursion_depth = max_uint32(++current_->recursion, current_->max_recursion_depth);
		return;
	}

	Record* next = findChild(current_, func_name);
	if (next == nullptr) {
		// didn't find a child so create one
		next = new Record;
		next->parent = current_;
		next->func_name = func_name;
		addChild(current_, next);
	}

	current_ = next;
}

void Profiler::closeRecord(const char* func_name, const char* ident, uint64_t ns) {
	if (func_name != current_->func_name) {
		printf("%s != %s\n", func_name, current_->func_name);
		assert(func_name == current_->func_name);
	}
	if (current_->recursion > 0) {
		--current_->recursion;
		return;
	}

	current_->identifier = ident;
	uint64_t new_total = current_->total_ns + ns;
	if (new_total < current_->total_ns) {
		current_->total_overflowed = true;
		new_total = ns;
		current_->total_count = 0;
	}
	current_->total_ns = new_total;
	++current_->total_count;
	current_->avg_ns = current_->total_ns / current_->total_count;

	current_ = current_->parent;
}

Profiler::Record* Profiler::findChild(Record* node, const char* func_name) {
	if (node == nullptr)
		return nullptr;

	Record* child = node->children;
	while (child) {
		if (child->func_name == func_name)
			break;
		child = child->siblings;
	}
	return child;
}

void Profiler::addChild(Record* node, Record* child) {
	if (node->children == nullptr) {
		node->children = child;
		return;
	}

	Record* sibling = node->children;
	while (sibling->siblings) {
		sibling = sibling->siblings;
	}
	sibling->siblings = child;
}
