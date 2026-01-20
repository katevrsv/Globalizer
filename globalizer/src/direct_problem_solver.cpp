#include <iostream>
#include "input_data.h"
#include "direct_problem_solver.h"

double DirectProblemSolver::calc_y_and_get_error(const std::vector<double>& k,
	const std::vector<double>& y_exp, std::vector<double>& y_calc, std::vector<double>& W, const input_data& input,
	const std::vector<double>& degrees)
{
	std::array<std::vector<double>, n_KM_coeffs> KM_coeffs;  // Коэффициенты Кутты-Мерсона

	double dt = 0.001;
	double dt_prev = dt;
	const double eps = 1e-6;
	const double dt_min = 1e-9;  // Минимальный шаг
	const double dt_max = 0.1;   // Максимальный шаг
	double current_time = 0.;
	double total_time = input.t_exp;

	std::vector<double> y_prev = input.y0_conv;
	std::vector<double> y_prev_prev = input.y0_conv;  // Для модификации K3
	std::vector<double> y_current(input.n_dims);
	bool dt_changed = false;
	bool is_first_iteration = true;

	while (current_time < total_time)
	{
		if (!input.variable_parameter) // Подбираются константы скоростей
			get_KM_coeffs(k, y_prev, y_exp, dt, KM_coeffs, W, input, input.degrees, dt_changed, is_first_iteration);
		else // Подбираются порядки
			get_KM_coeffs(k, y_prev, y_exp, dt, KM_coeffs, W, input, degrees, dt_changed, is_first_iteration);

		std::vector<double> y_rk(input.n_dims); // y, вычисленное по 4 коэффициентам
		for (size_t i = 0; i < input.n_dims; ++i)
			y_rk[i] = y_prev[i] + (KM_coeffs[0][i] - 3 * KM_coeffs[2][i] + 4 * KM_coeffs[3][i]) / 2;

		if (!input.variable_parameter) // Подбираются константы скоростей
			KM_coeffs[4] = get_diff(k, y_rk, y_exp, W, input, input.degrees);
		else // Подбираются порядки
			KM_coeffs[4] = get_diff(k, y_rk, y_exp, W, input, degrees);

		for (size_t i = 0; i < input.n_dims; ++i)
			KM_coeffs[4][i] *= dt;

		for (size_t i = 0; i < input.n_dims; ++i)  // Рассчитываем концентрации на текущем шаге через концентрации на предыдущем шаге и коэффициенты Кутты-Мерсона
		{
			y_current[i] = y_prev[i] + (KM_coeffs[0][i] + 4 * KM_coeffs[3][i] + KM_coeffs[4][i]) / 6.;
		}

		/*double e_max = 0.0;
		for (size_t i = 0; i < input.n_dims; ++i)
		{
			double e_i = std::fabs(y_rk[i] - y_current[i]);
			if (e_i > e_max)
				e_max = e_i;
		}
		if (e_max > eps)
		{
			dt = std::max(dt / 2.0, dt_min);
			dt_changed = true;
			if (dt <= dt_min)
			{
				if (e_max < 2.0 * eps)  // Порог можно настроить
				{
					//std::cout << "Warning: Accepting step at dt_min with e_max=" << e_max << " (slightly above eps). Proceeding." << std::endl;
					dt_changed = false;  // Принимаем шаг
				}
				else
				{
					//std::cerr << "Error: dt reached dt_min, e_max=" << e_max << " >> eps. Cannot converge." << std::endl;
					break;  // Или throw
				}
			}
			else
			{
				continue;  // Продолжаем уменьшать dt
			}
		}
		else if (e_max < eps) {
			dt = std::min(dt * 2.0, dt_max);
			dt_changed = true;
		}
		else
		{
			dt_changed = false;
		}

		dt_changed = false;

		if (!dt_changed && !is_first_iteration)
		{
			KM_coeffs[2] = get_diff(k, y_prev_prev, y_exp, W, input, degrees);
			for (size_t i = 0; i < input.n_dims; ++i)
				KM_coeffs[2][i] = KM_coeffs[2][i] * dt + 8 * (KM_coeffs[0][i] + y_prev_prev[i] - y_prev[i]) / 3;
		}*/

		y_prev_prev = y_prev;
		y_prev = y_current;
		current_time += dt;

		is_first_iteration = false;

		if (std::fabs(current_time - input.t_exp) < dt / 2.0) {
			y_calc = y_current;
		}
	}

	double error = 0.;
	if (y_calc.size() > 0)
	{
		for (size_t i = 0; i < input.n_dims; ++i)
		{
			error += input.weights[i] * fabs(y_calc[i] - y_exp[i]);
		}
	}
	return error;
}

// без учета m и G
/*dy_i/dl = SUM(smm_ij * W_j) */
std::vector<double> DirectProblemSolver::get_diff(const std::vector<double>& k, const std::vector<double>& y, const std::vector<double>& y_exp,
	std::vector<double>& W, const input_data& input, const std::vector<double>& degrees)
{
	// Получаем значение производной для каждой из n_dims концентраций
	std::vector<double> diff(y.size());

	int d = 0;
	for (size_t i = 0; i < input.n_stages; ++i)
	{
		W[i] = k[i];
		for (size_t j = 0; j < input.n_dims; ++j)
		{
			if (input.smm[i][j] < 0)
			{
				W[i] *= pow(y[j], degrees[d]);
				d++;
			}
		}

		if (i == 0 && input.n_stages == 5)
		{
			double denominator = (y[1] + 1. * sqrt(y[3])) * (y[1] + 1. * sqrt(y[3]));
			W[i] = k[i] * y[0] * y[1] / denominator;
		}
		if (i == 2 && input.n_stages == 5)
		{
			double denominator = (1. * y[1] + 1. * sqrt(y[3])) * (1. * y[1] + 1. * sqrt(y[3]));
			W[i] = k[i] * y[3] * y[5] / denominator;
		}
		if (i == 1 && input.n_stages == 2)
		{
			double denominator = (1. * y[1] + 1. * sqrt(y[3])) * (1. * y[1] + 1. * sqrt(y[3]));
			W[i] = k[i] * y[0] * y[3] / denominator;
		}
	}

	// Корректировка для одностадийной реакции, чтобы порядок по H2 (y[1])
	// зависел от соотношения H2 и C2H6
	if (input.does_H2_degree_depend && input.n_stages == 1)
	{
		const std::vector<double>& y_for_degree = input.does_H2_degree_depend == 1 ? y : y_exp;
		if (y_for_degree[1] / y_for_degree[0] < 7.55) // H2 / C2H6 < 7.55
			W[input.stage_with_H2_input] = k[input.stage_with_H2_input] * y[0] * pow(y[1], 1);
		else // H2/C2H6 >= 7.55
			W[input.stage_with_H2_input] = k[input.stage_with_H2_input] * y[0] * pow(y[1], -1.);
	}

	// Корректировка для двухстадийной, трехстадийной и четырехстадийной реакций, чтобы порядок по H2 (y[3])
	// зависел от соотношения H2 и C2H6
	/*if (input.does_H2_degree_depend && (input.n_stages == 2 || input.n_stages == 3 || input.n_stages == 4))
	{
		const std::vector<double>& y_for_degree = input.does_H2_degree_depend == 1 ? y : y_exp;
		if (y_for_degree[3] / y_for_degree[0] < 7.55) // H2 / C2H6 < 7.55
			W[input.stage_with_H2_input] = k[input.stage_with_H2_input] * y[0] * pow(y[3], 1);
		else // H2/C2H6 >= 7.55
			W[input.stage_with_H2_input] = k[input.stage_with_H2_input] * y[0] * pow(y[3], -1);
	}*/

	double sum_smm_W = 0;
	for (size_t i = 0; i < input.n_dims; ++i)
	{
		sum_smm_W = 0;
		for (size_t j = 0; j < input.n_stages; ++j)
			sum_smm_W += input.smm[j][i] * W[j];
		diff[i] = sum_smm_W;// *mol_mass[i] / G;
	}

	return diff;
}

void DirectProblemSolver::get_RK_coeffs(const std::vector<double>& k, const std::vector<double>& y_prev, const std::vector<double>& y_exp,
	const double dt, std::array<std::vector<double>, n_RK_coeffs>& RK_coeffs, std::vector<double>& W, const input_data& input)
{
	// Каждый из 4 коэффициентов Рунге-Кутты получаем для каждой из n_dims концентраций
	RK_coeffs[0] = get_diff(k, y_prev, y_exp, W, input, input.degrees);
	std::vector<double> y_tmp(input.n_dims);
	for (size_t i = 0; i < input.n_dims; ++i)
		y_tmp[i] = y_prev[i] + dt * RK_coeffs[0][i] / 2.;
	RK_coeffs[1] = get_diff(k, y_tmp, y_exp, W, input, input.degrees);
	for (size_t i = 0; i < input.n_dims; ++i)
		y_tmp[i] = y_prev[i] + dt * RK_coeffs[1][i] / 2.;
	RK_coeffs[2] = get_diff(k, y_tmp, y_exp, W, input, input.degrees);
	for (size_t i = 0; i < input.n_dims; ++i)
		y_tmp[i] = y_prev[i] + dt * RK_coeffs[2][i];
	RK_coeffs[3] = get_diff(k, y_tmp, y_exp, W, input, input.degrees);
}

void DirectProblemSolver::get_KM_coeffs(const std::vector<double>& k, const std::vector<double>& y_prev, const std::vector<double>& y_exp,
	const double dt, std::array<std::vector<double>, n_KM_coeffs>& KM_coeffs, std::vector<double>& W, const input_data& input, const std::vector<double>& degrees,
	bool dt_changed, bool is_first_iteration)
{
	std::vector<double> y_tmp(input.n_dims);

	KM_coeffs[0] = get_diff(k, y_prev, y_exp, W, input, degrees);
	for (size_t i = 0; i < input.n_dims; ++i)
		KM_coeffs[0][i] *= dt;

	for (size_t i = 0; i < input.n_dims; ++i)
		y_tmp[i] = y_prev[i] + KM_coeffs[0][i] / 3;
	KM_coeffs[1] = get_diff(k, y_tmp, y_exp, W, input, degrees);
	for (size_t i = 0; i < input.n_dims; ++i)
		KM_coeffs[1][i] *= dt;

	if (dt_changed || is_first_iteration)
	{
		for (size_t i = 0; i < input.n_dims; ++i)
			y_tmp[i] = y_prev[i] + (KM_coeffs[0][i] + KM_coeffs[1][i]) / 6;
		KM_coeffs[2] = get_diff(k, y_tmp, y_exp, W, input, degrees);
		for (size_t i = 0; i < input.n_dims; ++i)
			KM_coeffs[2][i] *= dt;
	}

	for (size_t i = 0; i < input.n_dims; ++i)
		y_tmp[i] = y_prev[i] + (KM_coeffs[0][i] + 3 * KM_coeffs[2][i]) / 8;
	KM_coeffs[3] = get_diff(k, y_tmp, y_exp, W, input, degrees);
	for (size_t i = 0; i < input.n_dims; ++i)
		KM_coeffs[3][i] *= dt;
}

// Получить значения констант балансных соотношений
std::vector<double> DirectProblemSolver::get_balance_consts(const std::vector<double>& y, const std::vector<std::vector<double>>& amm)
{
	const size_t n_dims = amm.size();
	const size_t n_elems = amm[0].size();
	std::vector<double> balance_consts(n_elems);
	for (size_t i = 0; i < n_elems; ++i)
	{
		balance_consts[i] = 0;
		for (size_t j = 0; j < n_dims; ++j)
			balance_consts[i] += amm[j][i] * y[j];
	}

	return balance_consts;
}

// Сравнить два массива
bool DirectProblemSolver::equal(const std::vector<double>& array1, const std::vector<double>& array2)
{
	const size_t n_dims = array1.size() < array2.size() ? array1.size() : array2.size();
	for (size_t i = 0; i < n_dims; ++i)
		if (fabs(array1[i] - array2[i]) > EPS)
			return false;

	return true;
}
