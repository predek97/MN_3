#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <iterator>
#include <string>
#include <tuple>
#include "interpolation.h"

// loads data with ommiting n values
auto load_data(std::string const &fileName) {
	auto data = std::vector<interpolation::point>();
	auto file = std::ifstream(fileName);
	double a, b;
	char delimiter;

	while (file >> a >> delimiter >> b)
		data.push_back(interpolation::point(a, b));

	return data;
}


int main(int argc, char **argv) {
	if (argc != 3) {
		std::cout << "Invalid arguments. Usage:\n program [input_filename] [step]\n"
			"Input file should contain values separated by comma"
			"Output files are: lagrange_input_filename, cubic_spline_input_filename\n";
		return -1;
	}
	auto input_filename = std::string(argv[1]);
	auto step = std::stod(argv[2]);
	auto data = ::load_data("data1.txt");
	// process and save data
	std::cout << "Performing lagrange interpolation..." << std::flush;
	interpolation::lagrange(data, "lagrange_" + input_filename, 1);
	std::cout << "DONE\n" << std::flush;
	return 0;
}