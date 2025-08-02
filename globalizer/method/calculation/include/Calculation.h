/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//             LOBACHEVSKY STATE UNIVERSITY OF NIZHNY NOVGOROD             //
//                                                                         //
//                       Copyright (c) 2015 by UNN.                        //
//                          All Rights Reserved.                           //
//                                                                         //
//  File:      calculation.h                                               //
//                                                                         //
//  Purpose:   Header file for calculation base class                      //
//                                                                         //
//  Author(s): Lebedev I.                                                  //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

#ifndef __CALCULATION_H__
#define __CALCULATION_H__

#include "Common.h"
#include "Parameters.h"
#include "Task.h"
#include "InformationForCalculation.h"
#include "Trial.h"
#include "SearchData.h"


/**
Базовый класс для проведения испытаний
*/
class Calculation
{
protected:

  /// Указатель на решаемую задачу
  Task* pTask;

  /// Данные по решению
  SearchData* pData;


  /// Количество вычислений за итерацию
  static int countCalculation;

  static bool isStartComputingAway;

  static InformationForCalculation inputCalculation;

  static TResultForCalculation resultCalculation;

public:

  /// Вычислитель который вызывался изначально
  static Calculation* firstCalculation;
  /// Вычислитель который вызывался на листьях
  static Calculation* leafCalculation;


  Calculation(Task& _pTask);

  void SetCountCalculation(int c);

  virtual ~Calculation() {}

  virtual void ContinueComputing() 
  {};

  /// Вычисляет функции fNumber и индекс невыполненного ограничения по координатам
  virtual void Calculate(InformationForCalculation& inputSet,
    TResultForCalculation& outputSet) = 0;
  /// Задасть используемую задачу
  void SetTask(Task* _pTask);

  void SetSearchData(SearchData* _pData);

  //Сброс данных
  virtual void Reset();
};

#endif
// - end of file ----------------------------------------------------------------------------------
