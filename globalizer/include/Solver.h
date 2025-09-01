#ifndef __SOLVER_H__
#define __SOLVER_H__

#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif

#include <mpi.h>
#include <exception>

#include "Common.h"
#include "GlobProcess.h"
#include "Exception.h"
#include "InitProblem.h"
#include "OutputSystem.h"
#include "Messages.h"
#include "SolutionResult.h"
#include "SolverInterface.h"

#ifdef _GLOBALIZER_BENCHMARKS
#include "IGlobalOptimizationProblem.h"
#include "GlobalOptimizationProblemManager.h"
#endif // _GLOBALIZER_BENCHMARKS



/**
* Базовые классы для решения задач глобальной оптимизации
**/
// ------------------------------------------------------------------------------------------------



// ------------------------------------------------------------------------------------------------
/**
Базового класса
*/
class Solver : public ISolver
{
protected:
  ///Процесс решающий задачу
  Process* mProcess;
  /// Задача
  IProblem* mProblem;

  /// Общее описание задачи
  Task* pTask;
  /// База данных(поисковая информация)
  SearchData* pData;

  /// Дерева распараллеливания
  //TParallelTree* parTree;
  /// Результат работы системы
  SolutionResult* result;

  /// Очистить данные
  virtual void ClearData();
  /** Инициализация чисел с расширенной точностью

  Функция автоматически выбирает используемый тип данных в зависимости от размерности решаемой
  задачи и плотности развертки
  */
  virtual void InitAutoPrecision();
  /// Создание процесса и всего остального
  virtual int CreateProcess();
  /// Проверяет что параметры солвера соответствуют решаемой задаче
  int CheckParameters();

  /// Решатель производящий вычисление при параллельном АГП (распараллеливание по точкам) на MPI
  void MpiCalculation();
  /// Решатель используемый при ассинхронной схеме
  void AsyncCalculation();

public:
  Solver(IProblem* problem);

  Solver(IGlobalOptimizationProblem* problem);

  /// Решение задачи по умолчанию
  virtual int Solve();
  /// Решение подзадачи с указанными параметрами
  virtual int Solve(Task* task);
  virtual ~Solver();

  void SetProblem(IProblem* problem);
  IProblem* GetProblem();
  SolutionResult* GetSolutionResult();
};

#endif //solver.h
