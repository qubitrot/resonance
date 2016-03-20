#pragma once
#ifndef HTML_WRITER_H_
#define HTML_WRITER_H_

#include "iwriter.h"
#include <cstring>

class HTMLWriter : public IWriter {
	public:
		HTMLWriter()
			: last_indent(1)
		{}

		void writeRecord(std::ostream& out, Profiler::Record* node,
				const std::string& indentation, float percent, float percent_total)
		{
			if (node->identifier && strcmp(node->identifier, "ROOT") == 0)
				return;

			if (last_indent < indentation.size()) {
				out << indentation << "<tr><td colspan=\"6\" style=\"padding-left:1cm;\">" << std::endl;
				out << indentation << "<table class=\""
					<< (indentation.size() <= 1 ? "func" : "sub_func")
					<< "\" width=\"100%\" cellpadding=\"0\" cellspacing=\"0\">" << std::endl;
				writeHeaderLine(out);
				last_indent = indentation.size();
			}

			if (last_indent > indentation.size()) {
				out << indentation << "</table>" << std::endl;
				out << indentation << "</td></tr>" << std::endl;
				last_indent = indentation.size();
			}

			out << indentation << "<tr";
			if (strcmp(node->func_name, "Unprofiled") == 0) {
				out << " bgcolor=\"#DA9BF0\""; // mmmmm.... tasty
			}
			out << "><td class=\"toggle\">";
			if (node->children)
				out << "<a href=\"#\">+</a>";
			out << "</td><td class=\"func_name\">";
			if (node->identifier) {
				out << "<span title=\"" << node->func_name << "\">" << node->identifier << "</span>";
			} else {
				out << node->func_name;
			}
			out.imbue(std::locale(std::cout.getloc(), new IWriter::comma_out));
			out << "</td><td class=\"total_count\">" << node->total_count
				<< "</td><td class=\"total_ns\">" << node->avg_ns
				<< "</td><td class=\"percent\">" << percent << '%'
				<< "</td><td class=\"total_percent\">" << percent_total << '%'
				<< "</td></tr>" << std::endl;
		}

		void writeHeader(std::ostream& out) {
			out << "<!DOCTYPE html>" << std::endl
				<< "<html>" << std::endl
				<< "<head>" << std::endl
				<< "<!-- Try loading both, because it could be either -->" << std::endl
				<< "<script src=\"jquery.js\"></script>" << std::endl
				<< "<script src=\"report.js\"></script>" << std::endl
				<< "<script src=\"Profiler/jquery.js\"></script>" << std::endl
				<< "<script src=\"Profiler/report.js\"></script>" << std::endl;
			out << "<style>" << std::endl
				<< "td             { padding: 0px 0px; border-left: solid black 1px; }" << std::endl
				<< ".toggle        { text-align: center; width: 1em; }" << std::endl
				<< ".toggle a      { color: black; text-decoration: none; }" << std::endl
				<< ".func_name     { text-align: left; padding-left: 5px; }" << std::endl
				<< ".total_count   { text-align: right; width: 2.4em; }" << std::endl
				<< ".total_ns      { text-align: right; width: 5em; }" << std::endl
				<< ".percent       { text-align: right; width: 4em; }" << std::endl
				<< ".total_percent { text-align: right; width: 4em; }" << std::endl
				<< ".header        { background-color: lightblue; text-align: right; font-weight: bold; }" << std::endl
				<< ".header td     { padding: 3px 5px; }" << std::endl
				<< ".sub_func      { border: solid black; border-width: 1px 0px; }" << std::endl
				<< ".func          { border: solid black; border-width: 1px 0px; }" << std::endl
				<< "</style>" << std::endl;
			out << "</head>" << std::endl
				<< "<body>" << std::endl;

			out << "<table style=\"border-right: solid black 1px;\" cellpadding=\"0\" cellspacing=\"0\">" << std::endl;
			writeHeaderLine(out);
		}

		void writeFooter(std::ostream& out) {
			out << "</table>" << std::endl
				<< "</body>" << std::endl
				<< "</html>" << std::endl;
		}

	private:
		void writeHeaderLine(std::ostream& out) {
				out << "<tr class=\"header\">"
					<< "<td class=\"toggle\"></td>"
					<< "<td class=\"func_name\">Function</td>"
					<< "<td class=\"total_count\"># Calls</td>"
					<< "<td class=\"total_ns\">Avg Nanosec</td>"
					<< "<td class=\"percent\">%</td>"
					<< "<td class=\"total_percent\">Total %</td>"
					<< "</tr>" << std::endl;
		}

	private:
		unsigned last_indent;
};

#endif
