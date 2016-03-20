#pragma once
#ifndef LATEX_WRITER_H_
#define LATEX_WRITER_H_

#include "iwriter.h"
#include <sstream>
#include <cstring>

class LaTeXWriter : public IWriter {
	public:
		void writeRecord(std::ostream& out, Profiler::Record* node, const std::string& indentation, float percent, float percent_total) {
			if (node->identifier && strcmp(node->identifier, "ROOT") == 0)
				return;

			bool total = false;
			bool unprofiled = false;

			if (strcmp(node->func_name, "Total") == 0)
				total = true;
			else if (strcmp(node->func_name, "Unprofiled") == 0)
				unprofiled = true;

			if (total)
				out << "\\rowcolor{red!15}" << std::endl;
			//else if (unprofiled)
			//	out << "\\rowcolor{purple!10}" << std::endl;
			int shade = 10000 - (25 * (indentation.size()-1));
			if (shade < 0)
				shade = 0;
			std::stringstream rowcolor;
			if ((shade / 100) % 2 == 0)
				rowcolor << "white!" << (shade%100) << "!darkgray";
			else
				rowcolor << "darkgray!" << (shade%100) << "!white";
			out << "\\rowcolor{" << rowcolor.str() << "!30}" << std::endl;

			out << "\\-\\hspace{" << (float)(indentation.size()-1)/2 << "cm}";
			out << "\\begin{lstlisting}" << std::endl;
			out << node->func_name << std::endl;
			out << "\\end{lstlisting}";
			out << " & " << node->total_count << " & " << node->avg_ns;
			out << " & " << percent << "\\%";
			out << "\\cellcolor{" << rowcolor.str() << "!" << int(100-percent) << "!halfred!30}";
			out << " & " << percent_total << "\\%";
			//out << "\\cellcolor{red!" << 35 * percent_total/100 << "}";
			out << "\\cellcolor{" << rowcolor.str() << "!" << int(100-percent_total) << "!red!30}";
			out << " \\\\" << std::endl;
		}

		void writeHeader(std::ostream& out) {
			out << "\\documentclass[11pt]{article}" << std::endl
				<< "\\usepackage{listings}" << std::endl
				//<< "\\usepackage{color}" << std::endl
				<< "\\usepackage[margin=1cm]{geometry}" << std::endl
				<< "\\usepackage[table]{xcolor}" << std::endl
				<< "\\begin{document}" << std::endl
				<< "\\lstset{language=C++}" << std::endl
				<< "\\colorlet{halfred}{red!65}" << std::endl;
			//out << "\\rowcolors{1}{gray}{white}" << std::endl;
			out << "\\begin{tabular*}{\\textwidth}{l | r | r | r | r\n}" << std::endl
				<< "\\textbf{Function} & \\textbf{Calls} & \\textbf{Avg Nanosec} & \\textbf{\\%} & \\textbf{\\%Total} \\\\" << std::endl
				<< "\\hline" << std::endl;
		}

		void writeFooter(std::ostream& out) {
			out << "\\end{tabular*}" << std::endl;
			out << "\\end{document}" << std::endl;
		}
};

#endif
