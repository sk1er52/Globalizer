#pragma once

#include "Task.h"

class SeparableOptimizationTask : public Task
{
protected:

  /// Начало блока переменных
  int startParameterNumber;


public:
  SeparableOptimizationTask(IProblem* _problem, int _ProcLevel);
  SeparableOptimizationTask();

  /// Создает копию класса
  virtual Task* Clone();

  /// Возвращает левую границу области поиска
  virtual const double* GetA() const;
  /// Возвращает правую границу области поиска
  virtual const double* GetB() const;

  /**
 Возвращает априори известные координаты точки глобального минимума
 Перед первым вызовом нужно вызвать resetOptimumPoint()
 */
  virtual const double* GetOptimumPoint() const;
  /// Вычисляет значение функции с номером fNumber в точке y
  virtual double CalculateFuncs(const double* y, int fNumber);
  /**
  Вычисляет numPoints значений функции с номером fNumber, в координатах y в массив values
  Работает только если задача является наследником IGPUProblem
  */
  virtual void CalculateFuncsInManyPoints(double* y, int fNumber, int numPoints, double* values);

  /// Задает начало блока переменных
  void SetStartParameterNumber(int _startParameterNumber);

  /**
  * \brief Копирует координаты точки из массива, согласно имеющимся правилам
  * \param[in] y имеющисся координаты.
  * \param[out] point точка назначения.
  * \return true, если значение допустимо, иначе false.
  */
  virtual void CopyPoint(double* y, Trial* point);

};