/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//             LOBACHEVSKY STATE UNIVERSITY OF NIZHNY NOVGOROD             //
//                                                                         //
//                       Copyright (c) 2016 by UNN.                        //
//                          All Rights Reserved.                           //
//                                                                         //
//  File:      init_problem.cpp                                            //
//                                                                         //
//  Purpose:   Source file for program                                     //
//                                                                         //
//  Author(s): Sysoyev A., Lebedev I., Sovrasov V.                         //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

#include "InitProblem.h"

int InitProblem(ProblemManager& problemManager, IProblem*& problem,
  int argc, char* argv[], bool isMPIInit)
{
  if (problemManager.LoadProblemLibrary(parameters.libPath) != ProblemManager::OK_)
  {
    //сообщение об ошибке печатает manager
    return 1;
  }

  IProblem* baseProblem = problemManager.GetProblem();
  CombinableBaseParameters* newProblem = dynamic_cast<CombinableBaseParameters*>(baseProblem);
  problem = baseProblem;
  if (newProblem == 0)
  {
    problem->SetConfigPath(parameters.libConfigPath);
    if (problem->Initialize() != ProblemManager::OK_)
    {
      printf("Error during problem initialization\n");
      return 1;
    }

    if (parameters.Dimension.GetIsChange())
    {
      //вообще, вызов SetDimension лучше убрать и получать размерность из конфигурационного файла для всех задач, где она не фиксирована
      if (problem->SetDimension(parameters.Dimension) != ProblemManager::OK_)
      {
        printf("Unsupported problem dimension!\n");
        return 1;
      }
    }
    //размерность задачи из конфигурационного файла имеет приоритет над значением из командной строки
    parameters.Dimension = problem->GetDimension();
  }
  else
  {
    problem->SetConfigPath(parameters.libConfigPath);
    newProblem->SetInitParam(argc, argv, isMPIInit);
    if (problem->Initialize() != ProblemManager::OK_)
    {
      printf("Error during problem initialization\n");
      return 1;
    }

    parameters.CombineOptions(newProblem->GetOptions(), newProblem->GetOptionsCount());
    newProblem->CombineOptions(parameters.GetOptions(), parameters.GetOptionsCount());

    newProblem->InitDataByParameters();
  }
  return 0;
}

#ifdef _GLOBALIZER_BENCHMARKS
// ------------------------------------------------------------------------------------------------
int InitProblemGlobalizerBenchmarks(GlobalOptimizationProblemManager& problemManager, IGlobalOptimizationProblem*& problem)
{
  std::string libPath = parameters.libPath;
  if (problemManager.LoadProblemLibrary(libPath) != GlobalOptimizationProblemManager::OK_)
  {
    //сообщение об ошибке печатает manager
    return 1;
  }

  IGlobalOptimizationProblem* baseProblem = problemManager.GetProblem();

  if (parameters.Dimension.GetIsChange())
    baseProblem->SetDimension(parameters.Dimension);
  else
    parameters.Dimension = baseProblem->GetDimension();

  baseProblem->Initialize();

  problem = baseProblem;
  return 0;
}
#endif // _GLOBALIZER_BENCHMARKS

// - end of file ----------------------------------------------------------------------------------
