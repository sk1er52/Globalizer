/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//             LOBACHEVSKY STATE UNIVERSITY OF NIZHNY NOVGOROD             //
//                                                                         //
//                       Copyright (c) 2015 by UNN.                        //
//                          All Rights Reserved.                           //
//                                                                         //
//  File:      MixedIntegerMethod.h                                        //
//                                                                         //
//  Purpose:   Header file for method class                                //
//                                                                         //
//  Author(s):  Lebedev.i                                                  //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////


/**
\file MixedIntegerMethod.h

\authors  Лебедев И.
\date 2025-2026
\copyright ННГУ им. Н.И. Лобачевского

\brief Объявление класса #MixedIntegerMethod

\details Объявление класса #MixedIntegerMethod и сопутствующих типов данных
*/


#ifndef __MIXED_INTEGER_METHOD_H__
#define __MIXED_INTEGER_METHOD_H__

#include "Method.h"


/**
Базовый класс, реализующий алгоритм глобального поиска.

В классе #Method реализованы основные функции, определяющие работу алгоритма глобального поиска.
*/
class MixedIntegerMethod : public Method
{
protected:



  // ----------------------------------------------------------------------------
  // Внутренние данные метода
  // ----------------------------------------------------------------------------



  ///// Входные данные для вычислителя, формирубтся в CalculateFunctionals()
  InformationForCalculation inputSet;
  ///// Выходные данные вычислителя, обрабатывается в CalculateFunctionals()
  TResultForCalculation outputSet;
  /// информация о данных текущей итерации
  SearchIteration iteration;
  /// Была получена точка в окрестности глобального оптимума
  bool isFoundOptimalPoint;


  /// достигнутая точность
  double            AchievedAccuracy;
  /** Коэффициент локальной адаптации

  Диапазон значений параметра alfa от 1 (глобальный) до 20 (локальный) поиск
  Рекомендуемое значение alfa = 15.
  */
  double alfa;

  /// Число вычисленных значений каждой функции
  std::vector<int> functionCalculationCount;

  /// нужно ли искать интервал
  bool isFindInterval;

  std::vector< double > globalM;

  ///Новая точка устанавливается в интервал принадлежащий окрестности локального минимума
  bool isSetInLocalMinimumInterval;

  /// количество точек вычисленных локальным методом
  int localPointCount;
  /// число запусков локально метода
  int numberLocalMethodtStart;
  /// Нужно останавливаться
  bool isStop;

  /**
  Количество дискретных значений
  Произведение числа значений всех дискретных переменных.
  Равно числу интервалов.
  */
  int mDiscreteValuesCount;
  /// Значения дискретных параметров
  std::vector< std::vector< double > > mDiscreteValues;
  /// Индекс первого дискретного параметра
  int startDiscreteVariable;


  /// найденные локальные минимумы
  std::vector<Trial*> localMinimumPoints;

  //=====================================================================================================================================================
  //Для методов локального уточнения нужны миксимумы
  double* Xmax;
  double* mu;
  SearchInterval* intervalXMax;
  bool isSearchXMax;
  //=====================================================================================================================================================

  /// Массив для сохранения точек для последующей печати и рисования
  std::vector<Trial*> printPoints;


  /// Вычисление координат точек испытания для основной\единственной развертки
  virtual void CalculateCurrentPoint(Trial& pCurTrialsj, SearchInterval* BestIntervalsj);


  /// Задать значения дискретного параметра
  virtual void SetDiscreteValue(int u, std::vector< std::vector <double> > dvs);

  /// Получаем поисковую информацию, важно для адаптивного метода
  virtual SearchData* GetSearchData(Trial* trial);


public:

  MixedIntegerMethod(Task& _pTask, SearchData& _pData,
    Calculation& _Calculation, Evolvent& _Evolvent);
  virtual ~MixedIntegerMethod();

  /** Функция выполняет первую итерацию метода
  */
  virtual void FirstIteration();

  /** Вычисления точек очередной итерации

Вычисленные точки очередной итерации записываются в массив #iteration.pCurTrials
*/
  virtual void CalculateIterationPoints();



};

#endif