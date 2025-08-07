/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//             LOBACHEVSKY STATE UNIVERSITY OF NIZHNY NOVGOROD             //
//                                                                         //
//                       Copyright (c) 2021 by UNN.                        //
//                          All Rights Reserved.                           //
//                                                                         //
//  File:      mpi_calculation.h                                           //
//                                                                         //
//  Purpose:   Header file for Async MPI calculation class                 //
//                                                                         //
//  Author(s): Lebedev I.                                                  //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

#ifndef __MPI_CALCULATION_ASYNC_H__
#define __MPI_CALCULATION_ASYNC_H__

#include "Calculation.h"

class MPICalculationAsync : public Calculation
{
protected:
  bool isFirst = true;

  std::vector<Trial*> vecTrials;

  /// MPI-номер потомка, закончившего решение выделенной задачи
  int ChildNumRecv; 

  /// "Внутренний" номер потомка среди всех потомков данного процесса
  int ChildNum;     

  void FirstStartCalculate(InformationForCalculation& inputSet, TResultForCalculation& outputSet);
  void StartCalculate(InformationForCalculation& inputSet, TResultForCalculation& outputSet);
  void RecieveCalculatedFunctional();

public:
  static void AsyncFinilize();

  MPICalculationAsync(Task& _pTask) : Calculation(_pTask)
  {
  }

  // Вычисляет функции (ограничения и критерий) и индекс невыполненного ограничения по координатам
  virtual void Calculate(InformationForCalculation& inputSet, TResultForCalculation& outputSet);
};

#endif
// - end of file ----------------------------------------------------------------------------------
