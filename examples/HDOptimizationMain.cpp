/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//             LOBACHEVSKY STATE UNIVERSITY OF NIZHNY NOVGOROD             //
//                                                                         //
//                       Copyright (c) 2025 by UNN.                        //
//                          All Rights Reserved.                           //
//                                                                         //
//  File:      HDOptimizationMain.cpp                                      //
//                                                                         //
//  Purpose:   Console version of Globalizer system                        //
//                                                                         //
//  Author(s): Lebedev I., Barkalov K.                                     //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////


#include "Globalizer.h"

// ------------------------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  GlobalizerInitialization(argc, argv);

#ifdef _GLOBALIZER_BENCHMARKS
  GlobalOptimizationProblemManager manager;
  IGlobalOptimizationProblem* problem = 0;
  if (InitGlobalOptimizationProblem(manager, problem, parameters.libPath))
  {
    print << "Error during problem initialization\n";
    return 0;
  }

  if (parameters.Dimension.GetIsChange())
    problem->SetDimension(parameters.Dimension);
  else
    parameters.Dimension = problem->GetDimension();

  std::vector<double> y(problem->GetDimension());
  std::vector<std::string> u;
  std::vector<double> values(problem->GetNumberOfFunctions());

  int err = problem->GetStartTrial(y, u, values);
  if (err == IGlobalOptimizationProblem::PROBLEM_OK)
  {
    parameters.startPoint.SetSize(problem->GetDimension());
    for (int i = 0; i < problem->GetDimension(); i++)
    {
      parameters.startPoint[i] = y[i];
    }
  }

  // Решатель
  HDSolver solver(problem);
  // Решаем задачу
  if (solver.Solve() != SYSTEM_OK)
    throw EXCEPTION("Error: solver.Solve crash!!!");

#endif


  MPI_Finalize();
  return 0;
}
// - end of file ----------------------------------------------------------------------------------
