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
	if (argc != 3  && argc != 4)
		return -1;
	auto input_filename = std::string(argv[1]);
	auto step = std::stod(argv[2]);
	auto interval = argc == 4 ? std::stoi(argv[3]) : 0;
	auto data = ::load_data(input_filename);
	// process and save data
	std::cout << "Interpolacja lagrange'a..." << std::flush;
	interpolation::lagrange(data, "lagrange_" + input_filename, 0.2, interval);
	std::cout << "OBLICZONO\nInterpolacja spline'ami..." << std::flush;
	interpolation::cubic_spline(data, "cubic_spline_" + input_filename);
	return 0;
}