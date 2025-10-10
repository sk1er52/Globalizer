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
#include "SeparableOptimizationSolver.h"
#include "GlobalizerProblem.h"

#include <algorithm>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <iostream>

#ifdef _GLOBALIZER_BENCHMARKS
#include "IGlobalOptimizationProblem.h"
#include "GlobalOptimizationProblemManager.h"
#endif // _GLOBALIZER_BENCHMARKS



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

  parameters.Init(argc, argv, true);
  if (!parameters.IsStart())
  {
    print << "Need command-line arguments!";
    return 0;
  }

  // Инициализация системы вывода и печати ошибок
  OutputMessage::Init(true, parameters.logFileNamePrefix, parameters.GetProcNum(),
    parameters.GetProcRank());


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


#ifdef _GLOBALIZER_BENCHMARKS
  GlobalOptimizationProblemManager manager;
  IGlobalOptimizationProblem* problem = 0;
  if (InitProblemGlobalizerBenchmarks(manager, problem))
  {
    print << "Error during problem initialization\n";
    return 0;
  }

  // Решатель
  Solver solver(problem);
  // Решаем задачу
  if (solver.Solve() != SYSTEM_OK)
    throw EXCEPTION("Error: solver.Solve crash!!!");


#endif


  MPI_Finalize();
  return 0;
}
// - end of file ----------------------------------------------------------------------------------
