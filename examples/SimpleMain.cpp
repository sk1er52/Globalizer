/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//             LOBACHEVSKY STATE UNIVERSITY OF NIZHNY NOVGOROD             //
//                                                                         //
//                       Copyright (c) 2015 by UNN.                        //
//                          All Rights Reserved.                           //
//                                                                         //
//  File:      Main.cpp                                                    //
//                                                                         //
//  Purpose:   Console version of Globalizer system                        //
//                                                                         //
//  Author(s): Sysoyev A., Barkalov K.                                     //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif

#include "Solver.h"
#include "GlobalizerProblem.h"

#include <algorithm>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <iostream>

// ------------------------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  parameters.Init(argc, argv, true);

  // Инициализация системы вывода и печати ошибок
  OutputMessage::Init(true, parameters.logFileNamePrefix, parameters.GetProcNum(),
    parameters.GetProcRank());


  parameters.Dimension = 2;

  auto problem = new ProblemFromFunctionPointers(parameters.Dimension, // размерность задачи
    std::vector<double>(parameters.Dimension, -2.2), // верхняя граница
    std::vector<double>(parameters.Dimension, 1.8), // нижняя граница
    std::vector<std::function<double(const double*)>>(1, [](const double* y)
      {
        double pi_ = 3.14159265358979323846;
        double sum = 0.;
        for (int j = 0; j < parameters.Dimension; j++)
          sum += y[j] * y[j] - 10. * cos(2.0 * pi_ * y[j]) + 10.0;
        return sum;
      }), // критерий
    true, // определен ли оптимум
    0, // значение глобального оптимума
    std::vector<double>(parameters.Dimension, 0).data() // координаты глобального минимума
  );

  problem->Initialize();

  // Решатель
  Solver solver(problem);

  // Решаем задачу
  if (solver.Solve() != SYSTEM_OK)
    throw EXCEPTION("Error: solver.Solve crash!!!");

  return 0;
}
// - end of file ----------------------------------------------------------------------------------
