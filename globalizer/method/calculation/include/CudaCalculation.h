/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//             LOBACHEVSKY STATE UNIVERSITY OF NIZHNY NOVGOROD             //
//                                                                         //
//                       Copyright (c) 2015 by UNN.                        //
//                          All Rights Reserved.                           //
//                                                                         //
//  File:      cuda_calculation.h                                          //
//                                                                         //
//  Purpose:   Header file for CUDA calculation class                      //
//                                                                         //
//  Author(s): Lebedev I.                                                  //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

#ifndef __CUDA_CALCULATION_H__
#define __CUDA_CALCULATION_H__

#include "Calculation.h"

class CUDACalculation : public Calculation
{
protected:
  ///Кол-во устройств
  int deviceCount;
  ///Кол-во данных на каждом устройстве
  int* dataSize;
  ///Начало данных для устройств
  int* dataStart;
  //индексы используемых устройств
  int* devicesIndex;

  double* coordinates;
  int coordinatesSize;
  double* FuncValues;
  int FuncValuesSize;

  bool mIsInitialized;

  
  void StartCalculate(InformationForCalculation& inputSet, TResultForCalculation& outputSet);

public:
  CUDACalculation(Task& _pTask) : Calculation(_pTask)
  {
    mIsInitialized = false;
	coordinates = 0;
    FuncValues = 0;
	coordinatesSize = 0;
	FuncValuesSize = 0;
  }

  //~CUDACalculation();

  /// 0 - Инициализирует массивы
  //virtual void Init(int _numPoint, int _N, int _numFunc);

  /// 1 - Задает координаты точки испытания
  //virtual void SetCoordinate(int _indexPoint, double* _coordinates, bool isCalculate = true);

  /// 2 - Вычисляет функции (ограничения и критерий) и индекс невыполненного ограничения по координатам
  virtual void Calculate(InformationForCalculation& inputSet, TResultForCalculation& outputSet);

  /// 3 - Возвращает вычисленные значения функций
  //virtual void GetFuncValue(int _indexPoint, double* _funcVal, int& index);
};

#endif
// - end of file ----------------------------------------------------------------------------------
