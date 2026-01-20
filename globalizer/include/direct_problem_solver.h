#ifndef DIRECT_PROBLEM_SOLVER_H
#define DIRECT_PROBLEM_SOLVER_H

#include <vector>
#include <array>
#include "input_data.h"

const double EPS = 1.e-05;     // Точность при проверке балансных соотношений
const size_t n_RK_coeffs = 4;  // Количество коэффициентов в методе Рунге-Кутты 4 порядка
const size_t n_KM_coeffs = 5;  // Количество коэффициентов в методе Кутты-Мерсона

/***************Дифференцирование по высоте слоя катализатора********************/
/*const size_t n_steps = 1000;   // Количество шагов для высоты слоя катализатора, для которых будем рассчитывать концентрации
const double dl = 0.001;       // Шаг по длине реактора (высоте слоя катализатора)
const double l_exp = 0.01;     // Длина реактора (высота слоя катализатора), для которой заданы экспериментальные данные
const double l_start = 0.;     // Начальная длина реактора, с которой начинаем рассчитывать концентрации
const size_t ind_exp = size_t((l_exp - l_start) / dl); // Индекс, по которому в массиве l будет находиться экспериментальная длина реактора
*/

/**************Дифференцирование по времени *************************************/
const size_t n_steps = 1500;   // Количество шагов по времени, для которых будем рассчитывать концентрации
const double t_start = 0.;     // Начальный момент времени, с которого начинаем рассчитывать концентрации


class DirectProblemSolver {
public:
	static double calc_y_and_get_error(const std::vector<double>& k,
		const std::vector<double>& y_exp, std::vector<double>& y_calc, std::vector<double>& W, const input_data& input,
		const std::vector<double>& degrees);
	// Рассчитать производные (правую часть СДУ)
	static std::vector<double> get_diff(const std::vector<double>& k, const std::vector<double>& y, const std::vector<double>& y_exp,
		std::vector<double>& W, const input_data& input, const std::vector<double>& calc_degrees);
	// Получить коэффициенты Рунге-Кутты
	static void get_RK_coeffs(const std::vector<double>& k, const std::vector<double>& y_prev, const std::vector<double>& y_exp,
		const double dl, std::array<std::vector<double>, n_RK_coeffs>& RK_coeffs, std::vector<double>& W, const input_data& input);
	// Получить коэффициенты Кутты-Мерсона
	static void get_KM_coeffs(const std::vector<double>& k, const std::vector<double>& y_prev, const std::vector<double>& y_exp,
		const double dl, std::array<std::vector<double>, n_KM_coeffs>& KM_coeffs, std::vector<double>& W, const input_data& input, const std::vector<double>& degrees,
		bool dt_changed = false, bool is_first_iteration = false);
	// Получить значения констант балансных соотношений
	static std::vector<double> get_balance_consts(const std::vector<double>& y, const std::vector <std::vector<double>>& amm);
	// Сравнить два вектора
	static bool equal(const std::vector<double>& array1, const std::vector<double>& array2);
};
#endif
