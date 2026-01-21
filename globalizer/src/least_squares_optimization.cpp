#include "least_squares_optimization.h"
#include <fstream>
#include <string>
#include <algorithm>

// Метод наименьших квадратов для системы уравнений с двумя неизвестными
void
least_squares_method(const std::vector<double>& a, const std::vector<double>& b, const std::vector<double>& f,
	double& x1, double& x2)
{
	// k1 * x1 + k2 * x2 = b1
	// k3 * x1 + k4 * x2 = b2

	// x1 = (b1 - k2 * x2) / k1
	// x2 = (b2 - k3 * b1 / k1) / (k4 - k2 * k3 / k1)

	double k1 = 0, k2 = 0, k3 = 0, k4 = 0, b1 = 0, b2 = 0;
	for (size_t t = 0; t < a.size(); ++t)
	{
		k1 += a[t] * a[t];
		k2 += a[t] * b[t];
		k3 += a[t] * b[t];
		k4 += b[t] * b[t];
		b1 += a[t] * f[t];
		b2 += b[t] * f[t];
	}

	x2 = (b2 - k3 * b1 / k1) / (k4 - k2 * k3 / k1);
	x1 = (b1 - k2 * x2) / k1;
}

void LeastSquaresOptimization::least_squares_optimize(const input_data& input, const std::vector<calc_data>& c_data,
	std::vector<std::vector<double>>& k_calc, std::vector<double>& k0, std::vector<double>& Ea,
	const size_t t_ind_first, const size_t t_ind_last)
{
	// Уравнение Аррениуса типа 0:
	// ln(k0) - (E/R) * (1/T) = ln(k)
	// Уравнение Аррениуса типа 1:
	// ln(k0) - (E/R) * (1/T) = ln(k) - n * ln(T/298)
	// 1 * x1 + x2 * b = f
	std::vector<double> a, b, f, x1(input.n_stages), x2(input.n_stages);
	for (size_t t = t_ind_first; t < t_ind_last; ++t)
	{
		a.push_back(1.); // Коэффициент при x1
		b.push_back(1. / input.t_kelvin[t]); // Коэффициент при x2
	}

	size_t defined_ind = 0;
	for (size_t i = 0; i < input.n_stages; ++i)
	{
		// Для стадии не определена константа скорости, требуется расчёт
		if (input.is_stage_defined.find(i) == input.is_stage_defined.end())
		{
			f.clear();
			for (size_t t = t_ind_first; t < t_ind_last; ++t)
			{
				if (!input.arrenius_type) // arrenius_type = 0
					f.push_back(log(fabs(c_data[t].k[i])));
				else
					f.push_back(log(c_data[t].k[i]) - input.arrenius_degree * log(input.t_kelvin[t] / T_st));
			}
			least_squares_method(a, b, f, x1[i], x2[i]);

			std::ofstream fout("output_" + std::to_string(i) + ".csv");
			std::locale loc("ru_RU.UTF-8");

			fout << "T, C;" << "T, K;" << "1/T, 1/K;" << "ln(k);" << "ln(k0) - (E/R) * (1/T);" << std::endl;
			for (size_t t = t_ind_first; t < t_ind_last; ++t)
			{
				std::string res = "";
				res += std::to_string(input.T[t]) + ";";
				res += std::to_string(input.t_kelvin[t]) + ";";
				res += std::to_string(b[t - t_ind_first]) + ";";
				res += std::to_string(f[t - t_ind_first]) + ";";
				res += std::to_string(x1[i] + x2[i] * b[t - t_ind_first]) + ";";
				std::replace(res.begin(), res.end(), '.', ',');
				fout << res << std::endl;
			}
			fout << std::endl;
			fout.close();

			k0[i] = exp(x1[i]); // x1 = ln(k0)
			Ea[i] = -R * x2[i] * j_to_cal / 1000.; // x2 = -Ea/R => Ea = -R*x2

			// Рассчитываем k по уравнению Аррениуса через k0 и (-E/R), полученные методом наименьших квадратов
			// через k, найденные методом поиска по шаблону
			for (size_t t = t_ind_first; t < input.n_temps; ++t)
			{
				// В x2 используется Ea, выраженная в Дж/моль
				k_calc[t][i] = k0[i] * exp(x2[i] / input.t_kelvin[t]);
			}
		}
		// Для стадии определена константа скорости, копируем из начальных данных
		else
		{
			k0[i] = input.k0_defined[defined_ind];
			Ea[i] = input.Ea_defined[defined_ind] * j_to_cal / 1000.;
			for (size_t t = t_ind_first; t < input.n_temps; ++t)
			{
				k_calc[t][i] = input.k_defined[t][defined_ind];
			}
			defined_ind++;
		}
	}
}
