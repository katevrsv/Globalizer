#ifndef CALC_DATA_H
#define CALC_DATA_H

#include <string>

// Расчетные данные, соответствующие температуре
struct calc_data
{
	std::vector<double> k; // Константы скоростей стадий реакции
	std::vector<double> y; // Расчетные концентрации
	double H2_degree; // Порядок по водороду
	std::vector<double> W; // Скорости стадий реакции
	std::vector<double> degrees; // Порядки веществ (для случая, когда k задаются, а подбираются порядки)
	std::wstring degrees_string;
};

#endif
