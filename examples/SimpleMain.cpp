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

enum ProblemName { RASTRIGIN, STRONGINC3_LAMBDA_EXPRESSION, STRONGINC3_FUNCTION_POINTER};

double StronginC3Functionals(const double* y, int fNumber)
{
  double res = 0.0;
  double x1 = y[0], x2 = y[1];
  switch (fNumber)
  {
  case 0: // constraint 1
    res = 0.01 * ((x1 - 2.2) * (x1 - 2.2) + (x2 - 1.2) * (x2 - 1.2) - 2.25);
    break;
  case 1: // constraint 2
    res = 100.0 * (1.0 - ((x1 - 2.0) / 1.2) * ((x1 - 2.0) / 1.2) -
      (x2 / 2.0) * (x2 / 2.0));
    break;
  case 2: // constraint 3
    res = 10.0 * (x2 - 1.5 - 1.5 * sin(6.283 * (x1 - 1.75)));
    break;
  case 3: // criterion
  {
    double t1 = pow(0.5 * x1 - 0.5, 4.0);
    double t2 = pow(x2 - 1.0, 4.0);
    res = 1.5 * x1 * x1 * exp(1.0 - x1 * x1 - 20.25 * (x1 - x2) * (x1 - x2));
    res = res + t1 * t2 * exp(2.0 - t1 - t2);
    res = -res;
  }
  break;
  }

  return res;
}

// ------------------------------------------------------------------------------------------------
int main(int argc, char* argv[])
{

  GlobalizerInitialization(argc, argv);

  parameters.Dimension = 2;
  ProblemName problemName = STRONGINC3_FUNCTION_POINTER;
  IProblem* problem = nullptr;

  if (problemName == RASTRIGIN)
  {
    problem = new ProblemFromFunctionPointers(parameters.Dimension, // размерность задачи
      std::vector<double>(parameters.Dimension, -2.2), // нижняя граница
      std::vector<double>(parameters.Dimension, 1.8), //  верхняя граница
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
      std::vector<double>(parameters.Dimension, 0) // координаты глобального минимума
    );
  }
  else if (problemName == STRONGINC3_LAMBDA_EXPRESSION)
  {
    problem = new ProblemFromFunctionPointers(parameters.Dimension, // размерность задачи
      {0.0, -1.0}, // нижняя граница
      {4.0, 3.0}, // верхняя граница
      std::vector<std::function<double(const double*)>>({ 
        [](const double* y) { return 0.01 * ((y[0] - 2.2) * (y[0] - 2.2) + (y[1] - 1.2) * (y[1] - 1.2) - 2.25); }, // ограничение 0
        [](const double* y) { return 100.0 * (1.0 - ((y[0] - 2.0) / 1.2) * ((y[0] - 2.0) / 1.2) - (y[1] / 2.0) * (y[1] / 2.0)); }, // ограничение 1
        [](const double* y) { return 10.0 * (y[1] - 1.5 - 1.5 * sin(6.283 * (y[0] - 1.75))); }, // ограничение 2
        [](const double* y) 
        { 
          double t1 = pow(0.5 * y[0] - 0.5, 4.0);
          double t2 = pow(y[1] - 1.0, 4.0);
          return -((1.5 * y[0] * y[0] * exp(1.0 - y[0] * y[0] - 20.25 * (y[0] - y[1]) * (y[0] - y[1]))) + t1 * t2 * exp(2.0 - t1 - t2));
        } // ограничение 2
        }), 
      true, // определен ли оптимум
      -1.489444, // значение глобального оптимума
      { 0.941176,  0.941176 } // координаты глобального минимума
    );
  }
  else if (problemName == STRONGINC3_FUNCTION_POINTER)
  {
    problem = new ProblemFromFunctionPointers(parameters.Dimension, // размерность задачи
      { 0.0, -1.0 }, // нижняя граница
      { 4.0, 3.0 }, // верхняя граница
      StronginC3Functionals, // задача
      4 // количество функций (3 ограничения + 1 критерий)      
    );
  }

  problem->Initialize();

  // Решатель
  Solver solver(problem);

  // Решаем задачу
  if (solver.Solve() != SYSTEM_OK)
    throw EXCEPTION("Error: solver.Solve crash!!!");

  return 0;
}
// - end of file ----------------------------------------------------------------------------------
