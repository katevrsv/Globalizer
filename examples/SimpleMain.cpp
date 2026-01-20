/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//             LOBACHEVSKY STATE UNIVERSITY OF NIZHNY NOVGOROD             //
//                                                                         //
//                       Copyright (c) 2025 by UNN.                        //
//                          All Rights Reserved.                           //
//                                                                         //
//  File:      SimpleMain.cpp                                              //
//                                                                         //
//  Purpose:   Console version of Globalizer system                        //
//                                                                         //
//  Author(s): Lebedev I., Barkalov K.,                                    //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////


#include "Globalizer.h"
#include "calc_data.h"
#include "input_data.h"
#include "FileManager.h"
#include "direct_problem_solver.h"

const double T0 = 273.15;

void convert_concentrations(input_data& input, need_convert convert_param);
void calc_k_defined(input_data& input);
void get_defined_stages_map(input_data& input);

// ------------------------------------------------------------------------------------------------
int main(int argc, char* argv[])
{

  GlobalizerInitialization(argc, argv);
  
  //Инициализация - чтение из файла, запись в input и т.д.
  //взять из "химического" main
  input_data input;

  parameters.MaxNumOfPoints[0] = 10;
  parameters.Dimension = 5; // Размерность задачи
  IProblem* problem = nullptr;

	//"C:\Users\Ms.evil_Cat\Desktop\Downloads\input_data_new_5_1_temp.txt"
  std::string input_path = "C:\\Users\\Ms.evil_Cat\\Desktop\\Downloads\\input_data_new_5_1_temp.txt";

  FileManager::read_input_data(input_path, input);

  need_convert convert_param = NEED_CONVERT_TO_RATIOS;
  if (argc > 3)
    convert_param = need_convert(atoi(argv[3]));

  convert_concentrations(input, convert_param);

	input.t_kelvin.resize(input.n_temps); // Температура в K
	for (size_t i = 0; i < input.n_temps; ++i)
		input.t_kelvin[i] = input.T[i] + T0;

	calc_k_defined(input); // Получаем k, рассчитанные через заданные k0 и Ea
	get_defined_stages_map(input);

	// Расчетные данные для разных температур
	std::vector<calc_data> c_data(input.n_temps);
	for (size_t t = 0; t < input.n_temps; t++)
		c_data[t].W.resize(input.n_stages);


	if (!input.n_degrees) // Если порядки не заданы пользователем, берутся порядки из стехимометрической матрицы
	{
		input.degrees.clear();
		for (size_t i = 0; i < input.n_stages; ++i)
			for (size_t j = 0; j < input.n_dims; ++j)
				if (input.smm[i][j] < 0)
					input.degrees.push_back(-input.smm[i][j]);
		input.n_degrees = input.degrees.size();
	}

	if (!input.weights.size())
	{
		input.weights.resize(input.n_dims, 1.);
	}

  problem = new ProblemFromFunctionPointers(parameters.Dimension, // размерность задачи
		{ 11, 0.2, 0.26, 1.9, 1.2 }, // нижняя граница
		{ 13, 8, 3.2, 6.7, 13.8 }, //  верхняя граница
    std::vector<std::function<double(const double*)>>(1, [&input, &c_data](const double* y)
      {
				std::vector<double> result(5);

				// k1: [11, 13] -> [1e11, 5e13]
				result[0] = 1e11 + (y[0] - 11.0) * (5e13 - 1e11) / (13.0 - 11.0);

				// k2: [0.2, 8] -> [20, 800]
				result[1] = 20.0 + (y[1] - 0.2) * 100.0;

				// k3: [0.26, 3.2] -> [26000, 320000]
				result[2] = 26000.0 + (y[2] - 0.26) * 100000.0;

				// k4: [1.9, 6.7] -> [1900, 6700]
				result[3] = 1900.0 + (y[3] - 1.9) * 1000.0;

				// k5: [1.2, 13.8] -> [12000, 138000]
				result[4] = 12000.0 + (y[4] - 1.2) * 10000.0;

				input.initial_k_base[0] =  { 0.1, 0.1, 0.1, 0.1, 0.1 };//result;

				double res = DirectProblemSolver::calc_y_and_get_error(input.initial_k_base[0], input.y_exp_conv[0], c_data[0].y, c_data[0].W,
				input, input.degrees);
					return res;
      }) // критерий 
  );

  problem->Initialize();

  // Решатель
  Solver solver(problem);

  // Решаем задачу
  if (solver.Solve() != SYSTEM_OK)
    throw EXCEPTION("Error: solver.Solve crash!!!");

	//Значения параметров для лучшей найденной точки (в нашем диапазоне)
	const double* resY = solver.GetSolutionResult()->BestTrial->y;

	std::vector<double> result(5);

	// k1: [11, 13] -> [1e11, 5e13]
	result[0] = 1e11 + (resY[0] - 11.0) * (5e13 - 1e11) / (13.0 - 11.0);

	// k2: [0.2, 8] -> [20, 800]
	result[1] = 20.0 + (resY[1] - 0.2) * 100.0;

	// k3: [0.26, 3.2] -> [26000, 320000]
	result[2] = 26000.0 + (resY[2] - 0.26) * 100000.0;

	// k4: [1.9, 6.7] -> [1900, 6700]
	result[3] = 1900.0 + (resY[3] - 1.9) * 1000.0;

	// k5: [1.2, 13.8] -> [12000, 138000]
	result[4] = 12000.0 + (resY[4] - 1.2) * 10000.0;

	input.initial_k_base[0] =  { 0.1, 0.1, 0.1, 0.1, 0.1 };//result;

	DirectProblemSolver::calc_y_and_get_error(input.initial_k_base[0], input.y_exp_conv[0],
		c_data[0].y, c_data[0].W, input, input.degrees);

	//result_y[t].resize(input.n_dims);
	std::cout << "t = " << input.T[0] << std::endl;
	for (int i = 0; i < input.n_dims; ++i) {
		std::cout << "y[" << i << "] = " << c_data[0].y[i] << std::endl;
		//result_y[t][i] = c_data[t].y[i];
	}

  return 0;
}

void convert_concentrations(input_data& input, need_convert convert_param)
{
	input.y0_conv.resize(input.n_dims);
	input.y_exp_conv.resize(input.n_temps);
	for (size_t i = 0; i < input.n_temps; ++i)
		input.y_exp_conv[i].resize(input.n_dims);

	if (convert_param == NOT_NEED_CONVERT)
	{
		for (size_t i = 0; i < input.n_dims; ++i)
			input.y0_conv[i] = input.y0[i];
		for (size_t t = 0; t < input.n_temps; ++t)
			for (size_t i = 0; i < input.n_dims; ++i)
				input.y_exp_conv[t][i] = input.y_exp[t][i];
	}
	else if (convert_param == NEED_CONVERT_TO_RATIOS)
	{
		for (size_t i = 0; i < input.n_dims; ++i)
			input.y0_conv[i] = input.y0[i] / 100.;
		for (size_t t = 0; t < input.n_temps; ++t)
			for (size_t i = 0; i < input.n_dims; ++i)
				input.y_exp_conv[t][i] = input.y_exp[t][i] / 100.;
	}
	else if (convert_param == NEED_CONVERT_TO_MOL_L)
	{
		const double P0 = 1.05;
		const double P = P0 * 1.e5; // Давление в Па
		const double R = 8314; // R выражена в (Па * л)/(моль * K), чтобы концентрации получались в моль/л
		double Sx = 0, Sy0 = 0, Sy = 0;
		std::vector<double>mol_masses(input.n_dims);
		std::vector<double>xi0(input.n_dims);
		std::vector<double>yi0(input.n_dims);
		std::vector<double>yi(input.n_dims);

		for (size_t i = 0; i < input.n_dims; ++i)
		{
			for (size_t j = 0; j < input.n_elems; ++j)
				mol_masses[i] += input.amm[i][j] * input.molar_masses[j];

			mol_masses[i] /= 1000.;
		}

		for (size_t i = 0; i < input.n_dims; ++i)
		{
			xi0[i] = input.y0[i] / 100.;
			Sx += xi0[i] * mol_masses[i];
		}

		for (size_t i = 0; i < input.n_dims; ++i)
		{
			yi0[i] = xi0[i] * mol_masses[i] / Sx;
			Sy0 += yi0[i] / mol_masses[i];
		}

		for (size_t i = 0; i < input.n_dims; ++i)
		{
			yi[i] = (yi0[i] / mol_masses[i]) * P0 / Sy0;
			input.y0_conv[i] = yi[i] * P / (R * (input.T[0] + T0));
		}

		for (size_t t = 0; t < input.n_temps; ++t)
		{
			Sx = 0;
			Sy0 = 0;
			Sy = 0;
			for (size_t i = 0; i < input.n_dims; ++i)
			{
				xi0[i] = input.y_exp[t][i] / 100.;
				Sx += xi0[i] * mol_masses[i];
			}

			for (size_t i = 0; i < input.n_dims; ++i)
			{
				yi0[i] = xi0[i] * mol_masses[i] / Sx;
				Sy0 += yi0[i] / mol_masses[i];
			}

			for (size_t i = 0; i < input.n_dims; ++i)
			{
				yi[i] = (yi0[i] / mol_masses[i]) * P0 / Sy0;
				input.y_exp_conv[t][i] = yi[i] * P / (R * (input.T[t] + T0));
			}
		}
	}
	else
		return;
}

void calc_k_defined(input_data& input)
{
	input.k_defined.resize(input.n_temps);
	for (size_t t = 0; t < input.n_temps; ++t)
	{
		input.k_defined[t].resize(input.n_defined_stages);
		for (size_t s = 0; s < input.n_defined_stages; ++s)
			input.k_defined[t][s] = input.k0_defined[s] * exp((-input.Ea_defined[s] / R) / input.t_kelvin[t]);
	}
}

void get_defined_stages_map(input_data& input)
{
	for (size_t i = 0; i < input.defined_stages.size(); ++i)
		input.is_stage_defined.insert(input.defined_stages[i]);
}

// - end of file ----------------------------------------------------------------------------------
