#include "interpolation.h"
#include <future>
#include <fstream>
#include <iterator>
#include <tuple>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/lu.hpp>


interpolation::point interpolation::lagrange_x(double x, std::vector<point> const &points) {
	double sum = 0.0;
	for (int i = 0; i < points.size(); i++) {
		double product = points[i].second;
		for (int j = 0; j < points.size(); j++) {
			if (j == i) continue;
			product *= (x - points[j].first) / (points[i].first - points[j].first);
		}
		sum += product;
	}
	return point(x, sum);
}


void interpolation::lagrange(std::vector<point> const &points, std::string const &output_filename, double interpolationStep) {
	auto interpolatedPoints = std::vector<point>();
	interpolatedPoints.reserve(static_cast<unsigned long long>(points.back().first/interpolationStep + 1.0));
	double progress = 0.0;
	int steps = (points.back().first - points[0].first)/interpolationStep;
	for (auto x = points[0].first; x < points.back().first; x += interpolationStep) {
		interpolatedPoints.push_back(lagrange_x(x, points));
	}
	auto output = std::ofstream(output_filename);
	for (auto &ip : interpolatedPoints)
		output << ip.first << ',' << ip.second << '\n';
}

auto interpolation::build_equations_matrices(const std::vector<point> &points) {
	// matrix size
	auto N = 1u + 1u              // edge conditionals
		+ points.size() - 2  // inner points first derivative continuity
		+ points.size() - 2  // inner points second derivative continuity
		+ points.size() - 1  // f(x0) = y0
		+ points.size() - 1; // f(x1) = y1s
	auto A = boost::numeric::ublas::matrix<double>(N, N, 0);
	auto B = boost::numeric::ublas::vector<double>(N, 0);

	for (auto i = 0u; i < points.size() - 1; i++) {
		auto[x0, y0] = points[i];
		auto[x1, y1] = points[i + 1];
		auto h = x1 - x0;
		// generate X
		B[4 * i] = y0;
		B[4 * i + 1] = y1;
		// 1. x0 value
		A(4 * i + 0, 4 * i + 0) = 1;                // a
													// 2. x1 value
		A(4 * i + 1, 4 * i + 0) = 1;                // a
		A(4 * i + 1, 4 * i + 1) = h;                // b
		A(4 * i + 1, 4 * i + 2) = std::pow(h, 2);   // c
		A(4 * i + 1, 4 * i + 3) = std::pow(h, 3);   // d

		if (i >= points.size() - 2) { continue; } // check if not edge
												  // 3. x1 first derivative continuity
		A(4 * i + 2, 4 * i + 0) = 0;                    // a
		A(4 * i + 2, 4 * i + 1) = 1;                    // b
		A(4 * i + 2, 4 * i + 2) = 2 * h;                // c
		A(4 * i + 2, 4 * i + 3) = 3 * std::pow(h, 2);   // d
		A(4 * i + 2, 4 * i + 5) = -1; // == b1
									  // 4. x1 second derivative continuity
		A(4 * i + 3, 4 * i + 2) = 2;
		A(4 * i + 3, 4 * i + 3) = 6 * h;
		A(4 * i + 3, 4 * i + 6) = -2;
	}
	// second derivative = 0 at the beginning and at the end
	auto h = points[points.size() - 1].first - points[points.size() - 2].first;
	A(N - 2, 2) = 1;
	A(N - 1, N - 2) = 2;
	A(N - 1, N - 1) = 6 * h;
	return std::tuple(std::move(A), std::move(B));
}

void
interpolation::cubic_spline(std::vector<point> const &points, std::string const &output_filename,
	double interpolation_step) {
	if (interpolation_step <= 0) { throw std::runtime_error("interpolation step is negative"); }

	// compute coefficients
	auto coefficients = std::vector<std::tuple<double, double, double, double>>();
	{
		auto[A, B] = build_equations_matrices(points);
		auto X = boost::numeric::ublas::vector<double>(B.size());
		auto P = boost::numeric::ublas::permutation_matrix<double>(B.size());
		boost::numeric::ublas::lu_factorize(A, P);
		boost::numeric::ublas::lu_substitute(A, P, X);
		coefficients.reserve(X.size() / 4);
		for (auto i = 0u; i < X.size(); i += 4) {
			coefficients.emplace_back(X[i], X[i + 1], X[i + 2], X[i + 3]);
		}
	}

	// lambda which returns interpolated y value in x point
	auto f = [&](auto const x) {
		auto predicate = [&x](auto const &elem) { return elem.first >= x; };
		auto index = std::max(
			static_cast<int> (std::find_if(points.begin(), points.end(), predicate) - points.begin() - 1), 0);
		auto x0 = points[index].first;
		auto[a, b, c, d] = coefficients[index];
		auto result = a + b * (x - x0) + c * std::pow(x - x0, 2) + d * std::pow(x - x0, 3);
		return result;
	};

	// compute and save results
	auto inter_out = std::ofstream(output_filename);
	for (auto x = points.front().first; x <= points.back().first; x += interpolation_step) {
		inter_out << x << ',' << f(x) << '\n' << std::flush;
	}
}