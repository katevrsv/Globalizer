#include "File_Manager.h"
#include "direct_problem_solver.h"
#include <vector>
#include <locale>
#include <codecvt>

void FileManager::print_results(const std::string& filename, const input_data& input, const std::vector<calc_data>& c_data,
	const std::vector<std::vector<double>>& k_calc, const std::vector<double>& k0, const std::vector<double>& Ea)
{
	std::wofstream fout(filename);
	std::locale loc("ru_RU.UTF-8");
	fout.imbue(loc);
	if (input.n_stages == 1) {
		//C2H6, H2, CH4
		fout << L"Tem;Eksper C2H6;Raschet C2H6;Eksper H2;Raschet H2;Eksper CH4;Raschet CH4;H2 degree;"
			L"Сумма модулей отклонений (функция ошибки);k, полученные методом шаблона как минимизирующие функцию ошибки;"
			L"k, полученные по уравнению Аррениуса через k0 и - E / R, полученные методом наименьших квадратов; 1 / T;"
			L"ln k для k, полученных методом шаблона; ln k для k1, полученных по уравнению Аррениуса\n";
	}
	/*if (input.n_stages == 1) {
		// C2H6, H2, N2, CH4
		fout << L"Tem;Eksper C2H6;Raschet C2H6;Eksper H2;Raschet H2;Eksper N2;Raschet N2;Eksper CH4;Raschet CH4;"
 L"Сумма модулей отклонений (функция ошибки);k, полученные методом шаблона как минимизирующие функцию ошибки;"
 L"k, полученные по уравнению Аррениуса через k0 и - E / R, полученные методом наименьших квадратов;1/T;"
 L"ln k для k, полученных методом шаблона;ln k для k1, полученных по уравнению Аррениуса\n";
	}*/
	else if (input.n_stages == 3) {
		fout << L"Tem;Eksper C2H6;Raschet C2H6;Eksper H2O;Raschet H2O;Eksper CO2;Raschet CO2;Eksper H2;Raschet H2;Eksper CH4;Raschet CH4;"
			L"H2 degree; Сумма модулей отклонений(функция ошибки); k1, полученные методом шаблона как минимизирующие функцию ошибки;"
			L"k2, полученные методом шаблона как минимизирующие функцию ошибки;k3, полученные методом шаблона как минимизирующие функцию ошибки;"
			L"k1, полученные по уравнению Аррениуса через k0 и - E / R, полученные методом наименьших квадратов;"
			L"k2, полученные по уравнению Аррениуса через k0 и -E/R, полученные методом наименьших квадратов;"
			L"k3, полученные по уравнению Аррениуса через k0 и - E / R, полученные методом наименьших квадратов; 1 / T;"
			L"ln k для k1, полученных методом шаблона;ln k для k2, полученных методом шаблона;ln k для k3, полученных методом шаблона;"
			L"ln k для k1, полученных по уравнению Аррениуса; ln k для k2, полученных по уравнению Аррениуса; ln k для k3, полученных по уравнению Аррениуса\n";
	}
	else if (input.n_stages == 2) {
		fout << L"Tem;Eksper C2H6;Raschet C2H6;Eksper H2O;Raschet H2O;Eksper CO2;Raschet CO2;Eksper H2;Raschet H2;Eksper CH4;Raschet CH4;"
			L"H2_degree;Сумма модулей отклонений (функция ошибки);k1, полученные методом шаблона как минимизирующие функцию ошибки;"
			L"k2, полученные методом шаблона как минимизирующие функцию ошибки;"
			L"k1, полученные по уравнению Аррениуса через k0 и -E/R, полученные методом наименьших квадратов;"
			L"k2, полученные по уравнению Аррениуса через k0 и - E / R, полученные методом наименьших квадратов;"
			L"1/T;"
			L"ln k для k1, полученных методом шаблона; ln k для k2, полученных методом шаблона;"
			L"ln k для k1, полученных по уравнению Аррениуса;ln k для k2, полученных по уравнению Аррениуса;"
			L"balans1_exp;"
			L"balans1_calc;"
			L"balans2_exp;"
			L"balans2_calc;"
			L"balans3_exp;"
			L"balans3_calc;"
			L"w1;"
			L"w2\n";
	}
	else if (input.n_stages == 4) {
		fout << L"Tem;Eksper C2H6;Raschet C2H6;Eksper H2O;Raschet H2O;Eksper CO2;Raschet CO2;Eksper H2;Raschet H2;"
			L"Eksper CH4;Raschet CH4;Eksper C2H6;Raschet C2H6;H2 degree;"
			L"Сумма модулей отклонений(функция ошибки); k1, полученные методом шаблона как минимизирующие функцию ошибки;"
			L"k2, полученные методом шаблона как минимизирующие функцию ошибки;k3, полученные методом шаблона как минимизирующие функцию ошибки;"
			L"k4, полученные методом шаблона как минимизирующие функцию ошибки;"
			L"k1, полученные по уравнению Аррениуса через k0 и -E/R, полученные методом наименьших квадратов;"
			L"k2, полученные по уравнению Аррениуса через k0 и - E / R, полученные методом наименьших квадратов;"
			L"k3, полученные по уравнению Аррениуса через k0 и - E/R, полученные методом наименьших квадратов;"
			L"k4, полученные по уравнению Аррениуса через k0 и - E / R, полученные методом наименьших квадратов; 1 / T;"
			L"ln k для k1, полученных методом шаблона;ln k для k2, полученных методом шаблона;"
			L"ln k для k3, полученных методом шаблона; ln k для k4, полученных методом шаблона;"
			L"ln k для k1, полученных по уравнению Аррениуса;ln k для k2, полученных по уравнению Аррениуса;"
			L"ln k для k3, полученных по уравнению Аррениуса; ln k для k3, полученных по уравнению Аррениуса\n";
	}
	else if (input.n_stages == 5) {
		fout << L"Tem;Eksper C3H8;Raschet C3H8;Eksper H2O;Raschet H2O;Eksper CO2;Raschet CO2;Eksper H2;Raschet H2;"
			L"Eksper CH4;Raschet CH4;Eksper C2H6;Raschet C2H6;H2 degree;"
			L"Сумма модулей отклонений(функция ошибки); k1, полученные методом шаблона как минимизирующие функцию ошибки;"
			L"k2, полученные методом шаблона как минимизирующие функцию ошибки;k3, полученные методом шаблона как минимизирующие функцию ошибки;"
			L"k4, полученные методом шаблона как минимизирующие функцию ошибки;k5, полученные методом шаблона как минимизирующие функцию ошибки;"
			L"k1, полученные по уравнению Аррениуса через k0 и -E/R, полученные методом наименьших квадратов;"
			L"k2, полученные по уравнению Аррениуса через k0 и - E / R, полученные методом наименьших квадратов;"
			L"k3, полученные по уравнению Аррениуса через k0 и - E/R, полученные методом наименьших квадратов;"
			L"k4, полученные по уравнению Аррениуса через k0 и - E / R, полученные методом наименьших квадратов;"
			L"k5, полученные по уравнению Аррениуса через k0 и - E / R, полученные методом наименьших квадратов;"
			L"1 / T;"
			L"ln k для k1, полученных методом шаблона;ln k для k2, полученных методом шаблона;"
			L"ln k для k3, полученных методом шаблона; ln k для k4, полученных методом шаблона;"
			L"ln k для k5, полученных по уравнению Аррениуса;ln k для k1, полученных по уравнению Аррениуса;"
			L"ln k для k2, полученных по уравнению Аррениуса; ln k для k3, полученных по уравнению Аррениуса;"
			L"ln k для k4, полученных по уравнению Аррениуса; ln k для k5, полученных по уравнению Аррениуса;"
			L"balans1_exp;"
			L"balans1_calc;"
			L"balans2_exp;"
			L"balans2_calc;"
			L"balans3_exp;"
			L"balans3_calc;"
			L"w1;"
			L"w2;"
			L"w3;"
			L"w4;"
			L"w5;"
			L"degrees\n";
	}
	for (size_t t = 0; t < input.n_temps; ++t)
	{
		std::wstring res = L"";
		res += std::to_wstring(input.T[t]) + L";";
		std::vector<double> balance_consts_exp(input.n_elems);
		std::vector<double> balance_consts_calc(input.n_elems);
		double error = 0;
		for (size_t i = 0; i < input.n_dims; ++i)
		{
			res += std::to_wstring(input.y_exp_conv[t][i]) + L";";
			res += std::to_wstring(c_data[t].y[i]) + L";";
			error += fabs(input.y_exp_conv[t][i] - c_data[t].y[i]);
		}
		res += std::to_wstring(c_data[t].H2_degree) + L";";
		res += std::to_wstring(error) + L";";
		for (size_t i = 0; i < input.n_stages; ++i)
			res += std::to_wstring(c_data[t].k[i]) + L";";
		for (size_t i = 0; i < input.n_stages; ++i)
			res += std::to_wstring(k_calc[t][i]) + L";";
		res += std::to_wstring(1 / input.T[t]) + L";";
		for (size_t i = 0; i < input.n_stages; ++i)
			res += std::to_wstring(log(c_data[t].k[i])) + L";";
		for (size_t i = 0; i < input.n_stages; ++i)
			res += std::to_wstring(log(k_calc[t][i])) + L";";
		balance_consts_exp = DirectProblemSolver::get_balance_consts(input.y_exp_conv[t], input.amm);
		balance_consts_calc = DirectProblemSolver::get_balance_consts(c_data[t].y, input.amm);
		for (size_t i = 0; i < input.n_elems; ++i)
		{
			res += std::to_wstring(balance_consts_exp[i]) + L";";
			res += std::to_wstring(balance_consts_calc[i]) + L";";
		}
		for (size_t i = 0; i < input.n_stages; ++i)
		{
			res += std::to_wstring(c_data[t].W[i]) + L";";
		}
		if (input.variable_parameter == 1) // Подбираются порядки
			res += c_data[t].degrees_string;
		std::replace(res.begin(), res.end(), '.', ',');
		fout << res << std::endl;
	}

	fout << "\nk0; Ea, kkal/mol\n";
	std::wstring out_data = L"";
	for (size_t i = 0; i < input.n_stages; ++i)
	{
		out_data = std::to_wstring(k0[i]) + L";" + std::to_wstring(Ea[i]);
		std::replace(out_data.begin(), out_data.end(), '.', ',');
		fout << out_data << std::endl;
	}
	fout.close();
}

template <typename T>
void read_vector(std::ifstream& fin, const size_t n, std::vector<T>& v)
{
	std::string line;
	std::getline(fin, line);
	if (!line.size())
		return;
	std::stringstream ss(line);
	T x;
	for (size_t i = 0; i < n; ++i) {
		ss >> x;
		v.push_back(x);
	}
}

template <typename T>
void read_matrix(std::ifstream& fin, const size_t n_rows, const size_t n_cols, std::vector<std::vector<T>>& m)
{
	std::string line;
	m.resize(n_rows);
	for (size_t i = 0; i < n_rows; ++i)
		m[i].resize(n_cols);
	for (size_t i = 0; i < n_rows; ++i) {
		std::getline(fin, line);
		std::stringstream ss(line);
		for (size_t j = 0; j < n_cols; ++j)
			ss >> m[i][j];
	}
}

void process_keyword(const std::string& keyword, std::string& data, std::ifstream& fin, std::stringstream& ss, input_data& input)
{
	if (keyword == "NDIMS")
	{
		std::getline(fin, data);
		ss.clear();
		ss.str(data);
		ss >> input.n_dims;
	}
	else if (keyword == "NSTAGES")
	{
		std::getline(fin, data);
		ss.clear();
		ss.str(data);
		ss >> input.n_stages;
	}
	else if (keyword == "NTEMPS")
	{
		std::getline(fin, data);
		ss.clear();
		ss.str(data);
		ss >> input.n_temps;
	}
	else if (keyword == "NELEMS")
	{
		std::getline(fin, data);
		ss.clear();
		ss.str(data);
		ss >> input.n_elems;
	}
	else if (keyword == "MOLMASSES")
	{
		// Molar masses
		read_vector(fin, input.n_elems, input.molar_masses);
	}
	else if (keyword == "STOICHIOM")
	{
		// Stoichiometric matrix
		read_matrix(fin, input.n_stages, input.n_dims, input.smm);
	}
	else if (keyword == "ATOMMOLECUL")
	{
		// Atomic-molecular matrix
		read_matrix(fin, input.n_dims, input.n_elems, input.amm);
	}
	else if (keyword == "INITCONCEN")
	{
		// Concentrations in initial moment
		read_vector(fin, input.n_dims, input.y0);
	}
	else if (keyword == "TEMP")
	{
		// Temperature
		read_vector(fin, input.n_temps, input.T);
	}
	else if (keyword == "EXPTIME")
	{
		std::getline(fin, data);
		ss.clear();
		ss.str(data);
		// Experimental time moment
		ss >> input.t_exp;
	}
	else if (keyword == "EXPCONCEN")
	{
		// Experimental concentrations
		read_matrix(fin, input.n_temps, input.n_dims, input.y_exp);
	}
	else if (keyword == "H2STAGE")
	{
		std::getline(fin, data);
		ss.clear();
		ss.str(data);
		// Stage with H2 input
		ss >> input.stage_with_H2_input;
	}
	/*else if (keyword == "H2DEGREE")
	{
		std::getline(fin, data);
		ss.clear();
		ss.str(data);
		ss >> input.H2_degree_const;
	}*/
	else if (keyword == "H2OPTION")
	{
		std::getline(fin, data);
		ss.clear();
		ss.str(data);
		ss >> input.does_H2_degree_depend;
	}
	else if (keyword == "DEFSTAGES")
	{
		std::getline(fin, data);
		ss.clear();
		ss.str(data);
		// Count of defined stages
		ss >> input.n_defined_stages;
		// Stage numbers with defined k0 and Ea
		read_vector(fin, input.n_defined_stages, input.defined_stages);
		// k0 for defined stages
		read_vector(fin, input.n_defined_stages, input.k0_defined);
		// Ea for defined stages
		read_vector(fin, input.n_defined_stages, input.Ea_defined);
		// Переводим Ea из ккал/моль в Дж/моль
		for (size_t i = 0; i < input.n_defined_stages; ++i)
			input.Ea_defined[i] = input.Ea_defined[i] * 1000 / j_to_cal;
	}
	else if (keyword == "KRELATIONS")
	{
		input.does_k_depend = 1;
		std::getline(fin, data);
		ss.clear();
		ss.str(data);
		// Стадия, от k которой зависит k другой стадии
		ss >> input.n_k_independent;
		std::getline(fin, data);
		ss.clear();
		ss.str(data);
		// Стадия, k которой зависит от k другой стадии
		ss >> input.n_k_dependent;
		// Множители для отношений между k
		read_vector(fin, input.n_temps, input.mult_k);
	}
	else if (keyword == "ARRENTYPE")
	{
		std::getline(fin, data);
		ss.clear();
		ss.str(data);
		// Тип для уравнения Аррениуса
		ss >> input.arrenius_type;
	}
	else if (keyword == "ARRENDEGREE")
	{
		std::getline(fin, data);
		ss.clear();
		ss.str(data);
		// Степень для уравнения Аррениуса типа 1
		ss >> input.arrenius_degree;
	}
	else if (keyword == "INITKBASE")
	{
		// Initial k_base
		read_matrix(fin, input.n_temps, input.n_stages, input.initial_k_base);
	}
	else if (keyword == "DEGREES")
	{
		// Начальные порядки веществ
		input.n_degrees = 0;
		for (size_t i = 0; i < input.n_stages; ++i)
			for (size_t j = 0; j < input.n_dims; ++j)
				if (input.smm[i][j] < 0)
					input.n_degrees++;
		read_vector(fin, input.n_degrees, input.degrees);
	}
	else if (keyword == "WEIGHTS")
	{
		// Довески для веществ при оптимизации
		read_vector(fin, input.n_dims, input.weights);
	}
}

void FileManager::read_input_data(const std::string& file_name, input_data& input)
{
	std::ifstream fin(file_name);
	std::string keyword, data;
	std::stringstream ss;

	while (std::getline(fin, keyword))
	{
		if (!keyword.empty())
			process_keyword(keyword, data, fin, ss, input);
	}
}

void FileManager::print_results_for_direct_problem(const std::vector<double>& l, const std::vector<std::vector<double>>& y, const int print_step)
{
	std::ofstream fout("output.csv");

	for (size_t s = 0; s < y.size(); s += print_step)
	{
		fout << l[s] << ";";
		for (size_t i = 0; i < y[s].size(); ++i)
			fout << y[s][i] << ";";
		fout << std::endl;
	}
	fout.close();
}
