#include <boost/program_options.hpp>
#include <iostream>
#include "Shapes.h"

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
		// instantiate computation object
		BestFit *b = BestFitFactory::Create(type, numObs);

		if (b)
			{
			// next read from stdin x,y coord pairs (comma-delimited)
			std::string str;
			for (int i = 0; i < numObs; i++)
				{
				std::getline(std::cin, str, ',');
				double x = std::strtod(str.c_str(), NULL);
				std::getline(std::cin, str);
				double y = std::strtod(str.c_str(), NULL);

				// add to the computation object
				b->AddObservation(x, y);
				}

			// perform quasi-parametric least-squares adjustment
			b->SetVerbosity(verbosity);
			b->Compute();
			delete b;
			}
		else
			std::cout << "Invalid type, or insufficient observations.\n";

		return 0;
		}

	return 1;
}
