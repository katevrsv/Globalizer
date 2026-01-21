#ifndef LEAST_SQUARESOPTIMIZATION_H
#define LEAST_SQUARESOPTIMIZATION_H

#include "input_data.h"
#include "calc_data.h"
#include <vector>

class LeastSquaresOptimization {
public:
	static void least_squares_optimize(const input_data& input, const std::vector<calc_data>& c_data,
		std::vector<std::vector<double>>& k_calc, std::vector<double>& k0, std::vector<double>& Ea,
		const size_t t_ind_first, const size_t t_ind_last);
};
#endif