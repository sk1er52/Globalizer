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


#include <algorithm>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <iostream>

#ifndef WIN32
#include <unistd.h>
#endif


// ------------------------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  std::cout << "\n\n" << std::endl;
  for (int i = 1; i < argc; i++)
  {
    std::cout << argv[i] << " ";
  }
  std::cout << "\n\n" << std::endl;
  

  MPI_Init(&argc, &argv);

#ifdef ONE_MPI_PROCESS_PER_NODE
  int	mProcNum = -1;
  int mProcRank = -1;

  if (MPI_Comm_size(MPI_COMM_WORLD, &mProcNum) != MPI_SUCCESS)
  {
    throw ("Error in MPI_Comm_size call");
  }
  std::cout << "mProcNum = " << mProcNum << std::endl;
  if (MPI_Comm_rank(MPI_COMM_WORLD, &mProcRank) != MPI_SUCCESS)
  {
    throw ("Error in MPI_Comm_rank call");
  }
  std::cout << "mProcRank = " << mProcRank << std::endl;

  if (mProcRank != 0)
  {
    std::cout << "!!! mProcRank = " << mProcRank << "Exit !!!" << std::endl;
    return 0;
  }
#endif

  parameters.Init(argc, argv, true);
  if (!parameters.IsStart())
  {
    print << "Need command-line arguments!";
    return 0;
  }

  // Инициализация системы вывода и печати ошибок
  OutputMessage::Init(true, parameters.logFileNamePrefix, parameters.GetProcNum(),
    parameters.GetProcRank());

  ProblemManager manager;
  IProblem* problem = 0;

  //parameters.Dimension = 2;

  //problem = new ProblemFromFunctionPointers(parameters.Dimension, // размерность задачи
  //  std::vector<double>(parameters.Dimension, -2.2), // верхняя граница
  //  std::vector<double>(parameters.Dimension, 1.8), // нижняя граница
  //  std::vector<std::function<double(const double*)>>(1, [](const double* y)
  //    {
  //      double pi_ = 3.14159265358979323846;
  //      double sum = 0.;
  //      for (int j = 0; j < parameters.Dimension; j++)
  //        sum += y[j] * y[j] - 10. * cos(2.0 * pi_ * y[j]) + 10.0;
  //      return sum;
  //    }), // критерий
  //  true, // определен ли оптимум
  //  0, // значение глобального оптимума
  //  std::vector<double>(parameters.Dimension, 0).data() // координаты глобального минимума
  //  
  //  );
  //problem->Initialize();
  if (InitProblem(manager, problem, argc, argv, 1))
  {
    print << "Error during problem initialization\n";
    return 0;
  }

  if (parameters.GetProcRank() == 0 && !parameters.disablePrintParameters)
  {
    parameters.PrintParameters();

    unsigned long size = 256;
    char* CompName = 0;
#ifdef WIN32
    LPWSTR buffer = new wchar_t[size];
    for (unsigned long i = 0; i < size; i++)
      buffer[i] = 0;
    GetComputerNameW(buffer, &size);
    CompName = new char[size + 1];
    for (unsigned long i = 0; i < size; i++)
      CompName[i] = (char)buffer[i];
    CompName[size] = 0;
    size++;
#else
    char* hostname = new char[size];
    for (unsigned long i = 0; i < size; i++)
      hostname[i] = 0;
    gethostname(hostname, 256);
    size = (unsigned long)strlen(hostname);
    CompName = new char[size + 1];
    for (unsigned long i = 0; i < size; i++)
      CompName[i] = (char)hostname[i];
    CompName[size] = 0;
    size++;
#endif

    printf("%s\tProcRank=%d\tProcNum=%d\n", CompName, parameters.GetProcRank(), parameters.GetProcNum());
  }
  
//#ifdef WIN32
//  unexpected_function prev, cur;
//  cur = &Unexpected;
//  prev = set_unexpected(cur);
//  set_terminate(Terminate);
//#endif

  // Решатель
  Solver solver(problem);
  // Решаем задачу
  if (solver.Solve() != SYSTEM_OK)
    throw EXCEPTION("Error: solver.Solve crash!!!");


  MPI_Finalize();
  return 0;
}
// - end of file ----------------------------------------------------------------------------------
