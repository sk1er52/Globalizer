/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//             LOBACHEVSKY STATE UNIVERSITY OF NIZHNY NOVGOROD             //
//                                                                         //
//                       Copyright (c) 2015 by UNN.                        //
//                          All Rights Reserved.                           //
//                                                                         //
//  File:      omp_calculation.h                                           //
//                                                                         //
//  Purpose:   Header file for OpenMP calculation class                    //
//                                                                         //
//  Author(s): Lebedev I.                                                  //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

#ifndef __OMP_CALCULATION_H__
#define __OMP_CALCULATION_H__

#include "Calculation.h"

class OMPCalculation : public Calculation
{
protected:

  /// Проверяем что индекс в допустимом диапозоне, иначе перераспределяем память
  //void CheckSize(int _indexPoint);

  /// 1 - если значение функции было вычислено
  //int* calc;

  void StartCalculate(InformationForCalculation& inputSet, TResultForCalculation& outputSet);

public:
  OMPCalculation(Task& _pTask) : Calculation(_pTask)
  {
    //calc = 0;
  }

  // 0 - Инициализирует массивы
  //virtual void Init(int _numPoint, int _N, int _numFunc);

  // 1 - Задает координаты точки испытания
  //virtual void SetCoordinate(int _indexPoint, double* _coordinates, bool isCalculate = true);

  // Вычисляет функции (ограничения и критерий) и индекс невыполненного ограничения по координатам
  virtual void Calculate(InformationForCalculation& inputSet, TResultForCalculation& outputSet);

  // 3 - Возвращает вычисленные значения функций
  //virtual void GetFuncValue(int _indexPoint, double* _funcVal, int& index);

  //virtual void Clean();
};

#endif
// - end of file ----------------------------------------------------------------------------------
