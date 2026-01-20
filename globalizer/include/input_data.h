#ifndef INPUT_DATA_H
#define INPUT_DATA_H

#include <vector>
#include <unordered_set>

const double j_to_cal = 0.239005736; // Калорий в одном Джоуле
const double R = 8.314;    // Универсальная газовая постоянная, Дж / (моль * K)
const double T_st = 298;   // Температура, при которой обычно приводятся в справочниках термохимические данные, K

enum need_convert
{
	NOT_NEED_CONVERT = 0,  // не нужно преобразовывать
	NEED_CONVERT_TO_MOL_L, // нужно преобразовать из процентов в мол/л
	NEED_CONVERT_TO_RATIOS  // нужно преобразовать из процентов в доли
};

struct input_data
{
	size_t n_dims = 0; // Размерность системы. Совпадает с количеством химических реагентов
	size_t n_stages = 0; // Количество стадий в реакции. Оно же количество констант скоростей реакций
	size_t n_temps = 0; // Количество температур, для которых есть экспериментальные данные
	size_t n_elems = 0; // Количество химических элементов
	size_t n_degrees = 0; // Количество порядков, то есть веществ, стоящих слева от знака равенства
	double t_exp = 0.; // Момент времени, для которого известны экспериментальные данные
	std::vector<double> molar_masses;
	std::vector<std::vector<double>> smm; // Стехиометрическая матрица
	std::vector<std::vector<double>> amm; // Атомно-молекулярная матрица
	std::vector<double> y0; // Концентрации в начальный момент времени
	std::vector<double> y0_conv; // Концентрации в начальный момент времени, сконвертированные в моль/л
	std::vector<double> T; // // Температура в C
	std::vector<double> t_kelvin; // Температура в K
	std::vector<std::vector<double>> y_exp; // Экспериментальные данные для разных температур
	std::vector<std::vector<double>> y_exp_conv; // Экспериментальные данные для разных температур, сконвертированные в моль/л
	std::vector<double> balance_consts; // Константы балансных соотношений
	std::vector<std::vector<double>> initial_k_base; // Начальное приближение для k
	size_t n_defined_stages = 0;     // Количество стадий, для которых k0 и Ea уже рассчитаны
	std::vector<size_t> defined_stages; // Индексы стадий, для которых k0 и Ea уже рассчитаны
	std::vector<double> k0_defined;  // k0 для уже рассчитанных стадий
	std::vector<double> Ea_defined;  // Ea для уже рассчитанных стадий
	std::vector<std::vector<double>> k_defined; // k, определённые для уже рассчитанных стадий, через k0 и Ea
	std::unordered_set<size_t> is_stage_defined; // Множество номеров стадий, для которых константа скорости определена
	int stage_with_H2_input; // Номер стадии, у которой на входе H2
	//double H2_degree_const; // Заданная константа для порядка H2
	int does_H2_degree_depend = 0; // Зависит ли порядок H2 от отношения H2 / C2H6
	int does_k_depend = 0; // Зависят ли k друг от друга
	size_t n_k_independent; // Стадия, от k которой зависит k другой стадии
	size_t n_k_dependent; // Стадия, k которой зависит от k другой стадии
	std::vector<double> mult_k; // Множители для отношений между k
	int arrenius_type = 0; // 0: k = k0 * exp(-Ea/(R * T)), 1: k = k0 * ((T/298)^n) * exp(-Ea/(R * T))
	double arrenius_degree = 0;
	double tol = 1e-10;
	std::vector<double> k0; // Предэкспоненциальные множители (для случая, когда k задаются, а подбираются порядки)
	std::vector<double> Ea; // Энергии активации в Дж/моль (для случая, когда k задаются, а подбираются порядки)
	std::vector<double> degrees; // Начальные порядки веществ (для случая, когда k задаются, а подбираются порядки)
	std::vector<double> weights; // Довески для веществ при оптимизации
	int variable_parameter = 0; // Параметр, который подбирается. 0 - константы скорости. 1 - порядки
};

#endif
