/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//             LOBACHEVSKY STATE UNIVERSITY OF NIZHNY NOVGOROD             //
//                                                                         //
//                       Copyright (c) 2015 by UNN.                        //
//                          All Rights Reserved.                           //
//                                                                         //
//  File:      parameters.h                                                //
//                                                                         //
//  Purpose:   Header file for parameters class                            //
//                                                                         //
//  Author(s): Lebedev I.                                                  //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

/**
\file parameters.h

\authors Лебедев И.
\date 2015-2016
\copyright ННГУ им. Н.И. Лобачевского

\brief Объявление общих параметров

\details Объявление общих параметров, реализованны в виде класса #Parameters
*/


#ifndef __PARAMETERS_H__
#define __PARAMETERS_H__

#include "Messages.h"
#include "Common.h"

#include "BaseParameters.h"

#include <string>
#include <stdio.h>





class Parameters : public BaseParameters<Parameters>
{
#undef OWNER_NAME
#define OWNER_NAME Parameters

protected:
  /// Определить параметры MPI
  void DetermineProc();
  /// Расставить номера на устройсва
  void SetDeviceIndex();

  /**
  Задание значений по умолчанию для всех параметров
  Пример:
  InitOption(имя параметра, значение по умолчанию, "короткая команда", "справка по параметру", кол-во элементов);
  *кол-во элементов для не массивов всегда равно 1.
  InitOption(Separator,_, "-Separator", "eparator", 1);
  */
  virtual void SetDefaultParameters();

  /// Индекс процесса
  int mProcRank;
  /// Кол-во процессов
  int mProcNum;
  /// Необходимое число процессов
  int mNeedMPIProcessorCount;

public:
  /// Уровень текущего процесса
  int MyLevel;
  /// Номер
  int MyMap;
  /// последний сохраненный номер итерации
  int iterationNumber;

  //Параметры командной строки
  TInt<Parameters> NumPoints; // число точек, порождаемых методом на 1 итерации
  TInt<Parameters> StepPrintMessages; //
  /// Через какое количество итераций сохранять точки
  TInt<Parameters> StepSavePoint; 
  TETypeMethod<Parameters> TypeMethod;
  TETypeCalculation<Parameters> TypeCalculation;
  TETypeProcess<Parameters> TypeProcess; //
  TInt<Parameters> NumThread;
  TInt<Parameters> SizeInBlock; //размер CUDA блока
  /// Печатать ли отчетет в файл
  TBool<Parameters> IsPrintFile;
  /// Файл для печати результата, если "000" то не печатаем
  TString<Parameters> ResulLog; 
  TInt<Parameters> Dimension; //размерность исходной задачи
  TDouble<Parameters> r; // надежность метода (> 1)
  ///Добавка при динамичеки изменяемом r, r = r + rDynamic / (Iteration ^ (1/N))
  TDouble<Parameters> rDynamic;
  TDouble<Parameters> rEps; //параметр eps-резервирования
  TDouble<Parameters> Epsilon; //единая точность
  TString<Parameters> Comment; //Коментарий к эксперименту

  TDoubles<Parameters> M_constant;
  TInt<Parameters> m; // плотность построения развертки (точность 1/2^m по к-те)
  TInt<Parameters> deviceCount; //кол-во используемых ускорителей
  TEMapType<Parameters> MapType; // тип развертки (сдвиговая, вращаемая)
  TInt<Parameters> DebugAsyncCalculation; // Флаг для проверки работы асинхронной схемы, если не 0, то вычисления проводятся в строго заданном порядке

   /// Печатать ли информацию о сечении в многошаговой схеме
  TBool<Parameters> IsPrintSectionPoint;

  TInts<Parameters> MaxNumOfPoints; // максимальное число итераций для процессов на каждом уровне //  размер - NumOfProcLevels{100, 100, 100, 100};// // параметры метода
  TFlag<Parameters> HELP;
  TFlag<Parameters> IsPlot;
  TInt<Parameters> PlotGridSize;
  /// Число испытаний за итерацию будет вычисляться на каждой итерации в методе CalculateNumPoint()
  TFlag<Parameters> IsCalculateNumPoint;
  TBool<Parameters> IsSetDevice; //Назначать каждому процессу свое устройство (ускоритель)
  ///Индекс используемого устройства (ускорителей), если -1 используется первые deviceCount устройств
  TInt<Parameters> deviceIndex;

  TInt<Parameters> ProcRank;

  TELocalMethodScheme<Parameters> localVerificationType; //cпособ использования локального метода(только для синхронного типа процесса)

  /// Количество итераций локального метода
  TInt<Parameters> localVerificationIteration;
  /// Точность локального метода
  TDouble<Parameters> localVerificationEpsilon; 
  /// Количество точек точек параллельно вычисляемых локальным методом
  TInt<Parameters> localVerificationNumPoint;

  TInt<Parameters> localMix; //параметр смешивания в локально-глобальном алгоритме
  TDouble<Parameters> localAlpha; //степень локальной адаптации в локально-глобальном алгоритме
  TInts<Parameters> calculationsArray; //Распределение типов вычислений по
  TESeparableMethodType<Parameters> sepS;  //флаг сепарабельного поиска на первой итерации
  TBool<Parameters> rndS;  //флаг случайного поиска на первой итерации
  TString<Parameters> libPath;  //путь к библиотеке с задачей
  TString<Parameters> libConfigPath; //путь к библиотеке с задачей
  /// тип критерия остановки
  TEStopCondition<Parameters> stopCondition; 
  /// Критерий применим только к верхнему уровню или к любому (для адаптивной схемы)
  TBool<Parameters> isStopByAnyLevel;
  TString<Parameters> iterPointsSavePath; //путь, по которому будут сохранены многомерные точки, поставленные методом корневого процесса
  TFlag<Parameters> printAdvancedInfo; //флаг, включающий печать дополнительной статистики: оценки констант Гёльдера и значения функций в точке оптимума
  TFlag<Parameters> disablePrintParameters; //флаг, выключающий печать параметров при запуске системы
  TString<Parameters> logFileNamePrefix; //префикс в имени лог-файла

  TETypeSolver<Parameters> TypeSolver;
  TInts<Parameters> DimInTask; // размерности каждой из подзадач в режиме сепарабильного или сикуенсального поиска

  /// Размер блока, отправляемого в MPI другим процессам
  TInt<Parameters> mpiBlockSize;
  /// Использовать ли специальное вычисление R как характеристики задачи
  TBool<Parameters> isUseTaskR;
  /// Использовать глобальный пересчет характеристик при изменение M или Z
  TBool<Parameters> isUseFullRecount;
  /// Использовать ли специальное вычисление R как характеристики интервалов
  TBool<Parameters> isUseIntervalR;
  /// Использовать ли глобальное Z
  TBool<Parameters> isUseGlobalZ;
  /// Не использовать Z
  TBool<Parameters> isNotUseZ;

  /// Тип локального метода (0 - Хука-Дживас)
  TETypeLocalMethod<Parameters> TypeLocalMethod;

  /// Тип добавления точек локального уточнения (0 - как обычные точки, 1 - точки локального метода не учитываются в критерии остановки по точности)
  TETypeAddLocalPoint<Parameters> TypeAddLocalPoint;
  /// Максимальное Кол-во точек устанавлиемых локальным методом
  TInt<Parameters> maxCountLocalPoint;
  /// Вычислять ли значения функции в крайних точках интервала
  TBool<Parameters> isCalculationInBorderPoint;
  /// Тип локального уточнения: 0 - без него; 1 - минимаксное; 2 - адаптивное; 3 - адаптивно-минимаксное
  TELocalTuningType<Parameters> LocalTuningType;
  /// Параметр кси, используемый в локальном уточении
  TDouble<Parameters> ltXi;  

  /// Загружать начальныеточки из файла или распределять их равномерно
  TBool<Parameters> isLoadFirstPointFromFile;
  /// Путь откуда будут считаны начальные точки испытания
  TString<Parameters> FirstPointFilePath;
  /// Тип распределения начальных точек
  TETypeDistributionStartingPoints<Parameters> TypeDistributionStartingPoints;

  /// Множитель перед функцие определяющий минимизируем или максимизируем функцию
  TDoubles<Parameters> functionSignMultiplier;

  /// Начальная точка для решения задачи оптимизации
  TDoubles<Parameters> startPoint;

  /// Проверка правильности при изменение параметров
  virtual int CheckValueParameters(int index = 0);
  /// Возвращает номер текущего процесса
  int GetProcRank() const;
  /// Возвращает общее число процессов
  int GetProcNum();
  /// Возвращает имя файла для сохранения картинки построенных линий уровней
  std::string GetPlotFileName();
  /// Печать текущих значений параметров
  void PrintParameters();

  ///Печать текущих значений параметров в файл
  void PrintParametersToFile(FILE* pf);

  /// Инициализация параметров
  virtual void Init(int argc, char* argv[], bool isMPIInit = false);
  Parameters();
  Parameters(Parameters& _parameters);
  virtual ~Parameters();

  /// Является ли класс задачей
  virtual bool IsProblem();

};

extern Parameters parameters;

#endif
// - end of file ----------------------------------------------------------------------------------
