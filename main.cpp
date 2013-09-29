#include "BestFit.h"

#if !defined BEST_FIT_LIB

#include <boost/program_options.hpp>
#include <iostream>
#include <sstream>

namespace po = boost::program_options;

int main(int argc, char *argv[])
{
	// Declare the supported options.
	po::options_description desc("Allowed options");
	desc.add_options()
		("type,t", po::value<int>(), "type of equation (0=line, 1=circle, 2=ellipse)")
		("verbose,v", po::value<int>(), "0=simple, 1=default, 2=verbose, 3=more verbose, etc.")
		("help,h", "produce this message");

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);

	if (1 == argc || vm.count("help") || 0 == vm.count("type"))
		{
		std::cout << "Usage: best-fit [options]" << std::endl;
		std::cout << desc << std::endl;
		return 1;
		}

	int type = 0;
	if (vm.count("type"))
		type = vm["type"].as<int>();

	int verbosity = 1;
	if (vm.count("verbose"))
		verbosity = vm["verbose"].as<int>();

	// first read the number of observations to follow
	int numObs = 0;
	std::cin >> numObs;
	std::cin.precision(12);
	std::cout.precision(12);

	if (numObs > 0)
		{
		double *points = new double[numObs * 2];
		std::string str;
		char comma;

		// next read from stdin x,y coord pairs (comma-delimited)
		for (int i = 0; i < numObs; i++)
			{
			std::getline(std::cin, str);
			std::stringstream ss(str);
			ss >> points[i * 2 + 0] >> comma >> points[i * 2 + 1];
			}


		BestFitIO input, output;
		input.numPoints = numObs;
		input.points = &points[0];
		input.verbosity = verbosity;
		output.points = input.points;

		// instantiate computation object
		BestFit *b = BestFitFactory::Create(type, std::cout);

		// perform quasi-parametric least-squares adjustment
		if (b)
			b->Compute(input, output);
		else
			std::cout << "Invalid type.\n";

		delete [] points;

		return 0;
		}

	return 1;
}

#endif