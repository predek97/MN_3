#include "interpolation.h"
#include <future>
#include <fstream>
#include <iterator>
#include <tuple>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <Eigen/Dense>



interpolation::point interpolation::lagrange_x(double x, std::vector<point> const &points, int interval) {
	double sum = 0.0;
	if (interval == 0)
		return lagrange_x(x, points);
	if (x < interval || x >= points.size() - interval)
		return point(x, 0);
	for (int i = x-interval; i < x+interval; i++) {
		double product = points[i].second;
		for (int j = x-interval; j < x+interval; j++) {
			if (j == i) continue;
			product *= (x - points[j].first) / (points[i].first - points[j].first);
		}
		sum += product;
	}
	return point(x, sum);
}

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


void interpolation::lagrange(std::vector<point> const &points, std::string const &output_filename, double interpolationStep, int interval) {
	auto interpolatedPoints = std::vector<point>();
	interpolatedPoints.reserve(static_cast<unsigned long long>(points.back().first/interpolationStep + 1.0));
	double progress = 0.0;
	int steps = (points.back().first - points[0].first)/interpolationStep;
	for (auto x = points[0].first; x < points.back().first; x += interpolationStep) {
		interpolatedPoints.push_back(lagrange_x(x, points, interval));
	}
	auto output = std::ofstream(output_filename);
	for (auto &ip : interpolatedPoints)
		output << ip.first << ',' << ip.second << '\n';
}

void
interpolation::cubic_spline(std::vector<point> const &points, std::string const &output_filename, double interpolationStep) {
	auto size = 4 * (points.size() - 1);
	boost::numeric::ublas::matrix<double> A(size, size), B(size, 1);
	for(int i = 0;; i++){
		std::cout << "dupa" << std::endl;
		auto[x0, y0] = points[i];
		auto[x1, y1] = points[i + 1];
		auto h = x1 - x0;
		B(4 * i, 0) = y0;
		B(4 * i + 1, 0) = y1;
		B(4 * i + 2, 0) = 0;
		B(4 * i + 3, 0) = 0;
		A(4 * i, 4 * i) = 1;
		A(4 * i + 1, 4 * i) = 1;                // a
		A(4 * i + 1, 4 * i + 1) = h;                // b
		A(4 * i + 1, 4 * i + 2) = std::pow(h, 2);   // c
		A(4 * i + 1, 4 * i + 3) = std::pow(h, 3);   // d
		if (i >= points.size() - 2) { break; } // check if not edge
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
	auto const h = points[points.size() - 1].first - points[points.size() - 2].first;
	A(size - 2, 2) = 1;
	A(size - 1, size - 2) = 2;
	A(size - 1, size - 1) = 6 * h;

	auto X = boost::numeric::ublas::vector<double>(B.size1());
	auto P = boost::numeric::ublas::permutation_matrix<double>(B.size1());
	boost::numeric::ublas::lu_factorize(A, P);
	boost::numeric::ublas::lu_substitute(A, P, B);

	auto f = [&B, &points](auto const x) {
		auto index = 0u;
		while (points[index + 1].first < x) index++;
		auto x0 = points[index].first;
		return B(4 * index, 0) +                             // a
			B(4 * index + 1, 0) * (x - x0) +              // b
			B(4 * index + 2, 0) * std::pow(x - x0, 2) +   // c
			B(4 * index + 3, 0) * std::pow(x - x0, 3);    // d
	};

	// compute and save results
	auto inter_out = std::ofstream(output_filename);
	for (auto x = points.front().first; x <= points.back().first; x += interpolationStep) {
		std::cout << x << " " << f(x) << '\n';
		inter_out << x << ',' << f(x) << '\n' << std::flush;
	}
}