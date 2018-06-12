#ifndef MN_INTERPOLATION_INTERPOLATION_TOOLS_H
#define MN_INTERPOLATION_INTERPOLATION_TOOLS_H

#include <vector>
#include <string>
#include <tuple>

namespace interpolation {
	using point = std::pair<double, double>;

	point lagrange_x(double x, std::vector<point> const &points, int interval);
	point lagrange_x(double x, std::vector<point> const &points);


	void lagrange(std::vector<point> const &points,
		std::string const &output_filename,
		double interpolationStep = 1,
		int interval = 0);

	void cubic_spline(std::vector<point> const &points,
		std::string const &output_filename,
		double interpolationStep = 1);

}
#endif //MN_INTERPOLATION_INTERPOLATION_TOOLS_H