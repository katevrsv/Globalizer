#ifndef FILEMANAGER_H
#define FILEMANAGER_H

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include "input_data.h"
#include "calc_data.h"

class FileManager {
public:
  // Считать данные из файла
  static void read_input_data(const std::string& file_name, input_data& input);
  // Вывести в файл результаты решения прямой задачи
  static void print_results_for_direct_problem(const std::vector<double>& l, const std::vector<std::vector<double>>& y, const int print_step);
  // Вывести результаты в файл
  static void print_results(const std::string& filename, const input_data& input, const std::vector<calc_data>& c_data,
    const std::vector<std::vector<double>>& k_calc, const std::vector<double>& k0, const std::vector<double>& Ea);
};
#endif
