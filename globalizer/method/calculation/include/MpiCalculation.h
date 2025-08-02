/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//             LOBACHEVSKY STATE UNIVERSITY OF NIZHNY NOVGOROD             //
//                                                                         //
//                       Copyright (c) 2021 by UNN.                        //
//                          All Rights Reserved.                           //
//                                                                         //
//  File:      mpi_calculation.h                                           //
//                                                                         //
//  Purpose:   Header file for MPI calculation class                       //
//                                                                         //
//  Author(s): Lebedev I.                                                  //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

#ifndef __MPI_CALCULATION_H__
#define __MPI_CALCULATION_H__

#include "Calculation.h"

class MPICalculation : public Calculation
{
protected:

  bool isFirst = true;

  //Производим вычисления
  void StartCalculate(InformationForCalculation& inputSet, TResultForCalculation& outputSet);

  //Производим вычисления на концах изначального отрезка при нужных параметрах запуска
  //Метод не будет работать при n < 3
  void StartCalculateInBorder(InformationForCalculation& inputSet, TResultForCalculation& outputSet);

public:
  MPICalculation(Task& _pTask) : Calculation(_pTask)
  {
  }

  // Вычисляет функции (ограничения и критерий) и индекс невыполненного ограничения по координатам
  virtual void Calculate(InformationForCalculation& inputSet, TResultForCalculation& outputSet);
};

#endif
// - end of file ----------------------------------------------------------------------------------
