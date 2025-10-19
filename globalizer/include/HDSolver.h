#pragma once

#include "Solver.h"
#include "HDTask.h"


/**
 Базовый класс для решателя задач большой размерности
**/
class HDSolver : public ISolver
{
protected:

  /// Решатели для оптимизации по группам параметров
  std::vector< Solver*> solvers;
  /// Решатель для объединения всех остальных Решателей
  Solver* finalSolver;
  /// Размерости групп парамеров, по умолчанию по 1
  std::vector<int> dimensions;
  /// Решение задачи оптимизации
  SolutionResult* solutionResult;
  /// Размерность исходной задачи
  int originalDimension;
  /// Задачи для оптимизации по группам параметров
  std::vector <HDTask*> tasks;

  /// задача оптимизации
  IProblem* problem;

  /// Задачть размерности
  void SetDimentions(std::vector<int> _dimentions);

  /// Создать начальную точку решения задач
  void CreateStartPoint();

  /// заполняет основные поля класса
  void Construct();

public:
  HDSolver(IProblem* problem, std::vector<int> _dimentions = {});

#ifdef _GLOBALIZER_BENCHMARKS

  HDSolver(IGlobalOptimizationProblem* problem, std::vector<int> _dimentions = {});

#endif

  virtual ~HDSolver();

  /// Решение задачи по умолчанию
  virtual int Solve();



  SolutionResult* GetSolutionResult();

  /// Добавляет точки испытаний
  virtual void SetPoint(std::vector<Trial*>& points);
  /// Возврящает все имеющиеся точки испытаний
  virtual std::vector<Trial*>& GetAllPoint();
};