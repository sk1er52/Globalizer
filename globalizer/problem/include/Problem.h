/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//             LOBACHEVSKY STATE UNIVERSITY OF NIZHNY NOVGOROD             //
//                                                                         //
//                       Copyright (c) 2016 by UNN.                        //
//                          All Rights Reserved.                           //
//                                                                         //
//  File:      problem_interface.h                                         //
//                                                                         //
//  Purpose:   Header file for Globalizer problem interface                //
//                                                                         //
//                                                                         //
//  Author(s): Lebedev I.                                                  //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

/**
\file problem_interface.h
\authors Соврасов В.
\date 2016
\copyright ННГУ им. Н.И. Лобачевского
\brief Объявление абстрактного класса #TIProblem
\details Объявление абстрактного класса #TIProblem и сопутствующих типов данных
*/

#ifndef __PROBLEM_H__
#define __PROBLEM_H__

#include <vector>
#include <string>
#include <cmath>
#include <functional>
#include <iostream>

#include "Common.h"
#include "ProblemPar.hpp"
#include "ProblemInterface.h"
#include "BaseParameters.h"

#include "Exception.h"


/**
Базовый класс-интерфейс, от которого наследуются классы, описывающие задачи оптимизации.
В классе #TIProblem описаны прототипы методов, которые должны быть реализованы в подключамых модулях с задачами.
*/
template <class Owner>
class BaseProblem : virtual public IIntegerProgrammingProblem, public BaseParameters<Owner>, virtual public ProblemPar
{
#undef OWNER_NAME
#define OWNER_NAME BaseProblem
protected:

  /// Максимальная допустимая размерность
  int mMaxDimension;
  /// минимальная допустимая размерность
  int mMinDimension;


  /** Левая граница области поиска,
  используется если левая граница одна и таже для всех размерностей
  */
  double mLeftBorder;
  /** Правая граница области поиска,
  используется если правая граница одна и таже для всех размерностей
  */
  double mRightBorder;

  /// Число значений дискретного параметра
  int* mNumberOfValues;
  /// Значение по умолчанию для числа значений дискретного параметра
  int mDefNumberOfValues;

  /// Очищает номер текущего значения для дискретного параметра
  virtual void ClearCurrentDiscreteValueIndex(int** mCurrentDiscreteValueIndex)
  {
    if (NumberOfDiscreteVariable > 0)
    {
      if (*mCurrentDiscreteValueIndex != 0)
        delete[] * mCurrentDiscreteValueIndex;

      *mCurrentDiscreteValueIndex = new int[NumberOfDiscreteVariable];
      for (int i = 0; i < NumberOfDiscreteVariable; i++)
        (*mCurrentDiscreteValueIndex)[i] = 0;
    }
  }

  /// Проверка правильности после окончания чтения параметров
  virtual int CheckValue(int index = -1)
  {

    BaseParameters<Owner>::CheckValue(index);
    if ((Dimension < mMinDimension) || (Dimension > mMaxDimension))
    {
      Dimension = mMinDimension;
    }

    if (NumberOfDiscreteVariable > 0)
    {
      if (mNumberOfValues != nullptr)
      {
        delete[] mNumberOfValues;
        mNumberOfValues = nullptr;
      }

      mNumberOfValues = new int[NumberOfDiscreteVariable];
      for (int i = 0; i < NumberOfDiscreteVariable; i++)
      {
        mNumberOfValues[i] = mDefNumberOfValues;
      }
    }

    return 0;
  }

  /// Задание значений по умолчанию базовых параметров
  virtual void SetBaseDefaultParameters()
  {
    BaseParameters<Owner>::SetBaseDefaultParameters();

    int sizeVal = 1;
    int bd = mDim;
    Owner* aqqq = (Owner*)this;
    Dimension.InitializationParameterProperty(aqqq, &BaseProblem::GetDim, &BaseProblem::SetDim,
      &BaseProblem::CheckValue,
      this->mOptionsCount, this->Separator, sizeVal, "Dimension", "Problem dimension", "-N", "0");
    if (bd != 0)
      Dimension = bd;
    else
      Dimension = mMinDimension;
    mDim = Dimension;
    this->AddOption((BaseProperty<Owner>*)(&Dimension));

    IsPrintParameters.InitializationParameterProperty(aqqq, &BaseProblem::CheckValue,
      this->mOptionsCount, this->Separator, sizeVal,
      "IsPrintParameters", "Is Print Parameters", "-IsPrintParameters", "false");
    this->AddOption((BaseProperty<Owner>*)(&IsPrintParameters));

    NumberOfDiscreteVariable.InitializationParameterProperty(aqqq, &BaseProblem::CheckValue, this->mOptionsCount,
      this->Separator, sizeVal,
      "NumberOfDiscreteVariable", "Number Of Discrete Variable", "-NDV", "0");
    this->AddOption((BaseProperty<Owner>*)(&NumberOfDiscreteVariable));

    IsMultInt.InitializationParameterProperty(aqqq, &BaseProblem::CheckValue, this->mOptionsCount,
      this->Separator, sizeVal,
      "IsMultInt", "IsMultInt", "-IMI", "true");
    this->AddOption((BaseProperty<Owner>*)(&IsMultInt));

    Dimension.mIsEdit = true;
    IsPrintParameters.mIsEdit = true;
    NumberOfDiscreteVariable.mIsEdit = true;
    IsMultInt.mIsEdit = true;

  }
  /**
  Задание значений по умолчанию для всех параметров
  Пример:
  InitOption(имя параметра, значение по умолчанию, "короткая команда", "справка по параметру", кол-во элементов);
  *кол-во элементов для не массивов всегда равно 1.
  InitOption(Separator,_, "-Separator", "eparator", 1);
  */
  virtual  void SetDefaultParameters()
  {
  }
public:

  /// Используемая размерность
  TInt<Owner> Dimension;
  /// Печатать ли параметры задачи на консоль
  TBool<Owner> IsPrintParameters;
  /// Число дискретных параметров
  TInt<Owner> NumberOfDiscreteVariable;
  /// Умножать ли на пароболоидцелевую функцию
  TBool<Owner> IsMultInt;

  int GetDim() const
  {
    return mDim;
  }

  void SetDim(int dim)
  {
    mDim = dim;
  }

  ///// Код ошибки, возвращаемый, если операция завершена успешно
  // static const int OK = 0;
  // /** Код ошибки, возвращаемый методами #GetOptimumValue и #GetOptimumPoint,
  // если соответствующие параметры задачи не определены,
  // */
  // static const int UNDEFINED = -1;
  // /// Код ошибки, возвращаемый, если операция не выполнена
  static const int PROBLEM_ERROR = -2;

  BaseProblem() : BaseParameters<Owner>()
  {
    this->mIsInit = false;
    mMaxDimension = MaxDim;
    mMinDimension = 2;
    mNumberOfCriterions = 1;
    mNumberOfConstraints = 0;
    std::string configPath = "";
    mLeftBorder = -1.8;
    mRightBorder = 2.2;
    //mCurrentDiscreteValueIndex = 0;
    mNumberOfValues = 0;
    mDim = 0;
    mDefNumberOfValues = -1;
  }
  /** Задание пути до конфигурационного файла
  Данный метод должн вызываться перед #Initialize
  \param[in] configPath строка, содержащая путь к конфигурационному файлу задачи
  \return Код ошибки
  */
  virtual int SetConfigPath(const std::string& configPath)
  {
    if (this->mIsInit)
    {
      this->ConfigPath = configPath;
    }
    else
    {
      this->mConfigPath = configPath;
    }
    return IProblem::OK;
  }

  /** Метод задаёт размерность задачи
  Данный метод должен вызываться перед #Initialize. Размерность должна быть в
  списке поддерживаемых.
  \param[in] dimension размерность задачи
  \return Код ошибки
  */
  virtual int SetDimension(int dimension)
  {
    if (!this->mIsInit)
    {
      mDim = dimension;
      Dimension.SetIsPreChange(true);
    }
    //return BaseProblem::PROBLEM_ERROR;
    if (dimension >= mMinDimension && dimension <= mMaxDimension)
    {
      Dimension = dimension;
      return IProblem::OK;
    }
    else
      return BaseProblem::PROBLEM_ERROR;
  }
  ///Возвращает размерность задачи, можно вызывать после #Initialize
  virtual int GetDimension() const
  {
    return mDim;
  }

  /// Инициализация параметров
  virtual void Init(int argc, char* argv[], bool isMPIInit = false)
  {
    BaseParameters<Owner>::Init(argc, argv, false);
    if (IsPrintParameters)
    {
      this->PrintParameters();
    }
  }

  ///Инициализация задачи
  virtual int Initialize(int argc, char* argv[], bool isMPIInit = false)
  {
    Init(argc, argv, false);
    return IProblem::OK;
  }

  ///Инициализация задачи
  virtual int Initialize()
  {
    Init(0, 0, false);
    return IProblem::OK;
  }

  /** Метод возвращает границы области поиска
  */
  virtual void GetBounds(double* lower, double* upper)
  {
    if (this->mIsInit)
      for (int i = 0; i < Dimension; i++)
      {
        lower[i] = mLeftBorder;
        upper[i] = mRightBorder;
      }
  }
  /** Метод возвращает значение целевой функции в точке глобального минимума
  \param[out] value оптимальное значение
  \return Код ошибки (#OK или #UNDEFINED)
  */
  virtual int GetOptimumValue(double& value) const
  {
    return BaseProblem<Owner>::UNDEFINED;
  }
  /** Метод возвращает значение функции с номером index в точке глобального минимума
  \param[out] value оптимальное значение
  \return Код ошибки (#OK или #UNDEFINED)
  */
  virtual int GetOptimumValue(double& value, int index) const
  {
    if (index == GetNumberOfConstraints())
      return GetOptimumValue(value);
    return IProblem::UNDEFINED;
  }


  /** Метод возвращает число общее функций в задаче (оно равно число ограничений + число критериев)
  \return Число функций
  */
  virtual int GetNumberOfFunctions() const
  {
    return mNumberOfCriterions + mNumberOfConstraints;
  }

  /** Метод возвращает число ограничений в задаче
  \return Число ограничений
  */
  virtual int GetNumberOfConstraints() const
  {
    return mNumberOfConstraints;
  }
  /** Метод возвращает число критериев в задаче
  \return Число критериев
  */
  virtual int GetNumberOfCriterions() const
  {
    return mNumberOfCriterions;
  }

  ///Деструктор
  virtual ~BaseProblem()
  {
    delete[] mNumberOfValues;
  }

  virtual std::string ProblemName()
  {
    return std::string("This Problem\n");
  }

  /// Возвращает число дискретных параметров, дискретные параметры всегда последние в векторе y
  virtual int GetNumberOfDiscreteVariable()
  {
    return NumberOfDiscreteVariable;
  }
  /**
  Возвращает число значений дискретного параметра discreteVariable.
  GetDimension() возвращает общее число параметров.
  (GetDimension() - GetNumberOfDiscreteVariable()) - номер начальной дискретной переменной
  Для не дискретных переменных == -1
  */
  virtual int GetNumberOfValues(int discreteVariable)
  {
    if ((discreteVariable > GetDimension()) ||
      (discreteVariable < (GetDimension() - GetNumberOfDiscreteVariable())))
      return -1;
    if (mNumberOfValues == 0)
      return -1;
    return mNumberOfValues[discreteVariable - (GetDimension() - GetNumberOfDiscreteVariable())];
  }
  /**
  Определяет значения дискретного параметра с номером discreteVariable
  Возвращает код ошибки.
  \param[out] values массив, в который будут сохранены значения дискретного параметра
  */
  virtual int GetAllDiscreteValues(int discreteVariable, double* values)
  {
    if ((discreteVariable > GetDimension()) ||
      (discreteVariable < (GetDimension() - GetNumberOfDiscreteVariable())))
      return IIntegerProgrammingProblem::ERROR_DISCRETE_VALUE;
    int* mCurrentDiscreteValueIndex = 0;
    ClearCurrentDiscreteValueIndex(&mCurrentDiscreteValueIndex);

    // сбрасываем значение индекса текущего значения и задаем левую границу
    GetNextDiscreteValues(mCurrentDiscreteValueIndex, values[0], discreteVariable, -1);
    int numVal = GetNumberOfValues(discreteVariable);
    // определяем все остальные значения
    for (int i = 1; i < numVal; i++)
    {
      GetNextDiscreteValues(mCurrentDiscreteValueIndex, values[i], discreteVariable);
    }
    return IProblem::OK;
  }
  /**
  Определяет значения дискретного параметра с номером discreteVariable после номера previousNumber
  Возвращает код ошибки.
  \param[in] previousNumber - номер значения после которого возвращается значение
  -2 - значение по умолчанию, возвращает следующее значение
  -1 - возвращает после -1, т.е. левую границу области
  \param[out] value переменная в которую сохраняется значение дискретного параметра
  */
  virtual int GetNextDiscreteValues(int* mCurrentDiscreteValueIndex, double& value, int discreteVariable, int previousNumber = -2)
  {
    if ((discreteVariable > GetDimension()) ||
      (discreteVariable < (GetDimension() - GetNumberOfDiscreteVariable())) ||
      (mCurrentDiscreteValueIndex == 0) ||
      (mNumberOfValues == 0))
      return IIntegerProgrammingProblem::ERROR_DISCRETE_VALUE;
    // если -1 то сбрасываем значение текущего номера
    if (previousNumber == -1)
    {
      mCurrentDiscreteValueIndex[discreteVariable - GetNumberOfDiscreteVariable()] = 0;
      value = mLeftBorder;
      return IProblem::OK;
    }
    else if (previousNumber == -2)
    {
      double d = (mRightBorder - mLeftBorder) /
        (mNumberOfValues[discreteVariable - (GetDimension() - GetNumberOfDiscreteVariable())] - 1);
      mCurrentDiscreteValueIndex[discreteVariable - GetNumberOfDiscreteVariable()]++;
      value = mLeftBorder + d *
        mCurrentDiscreteValueIndex[discreteVariable - GetNumberOfDiscreteVariable()];
      return IProblem::OK;
    }
    else
    {
      double d = (mRightBorder - mLeftBorder) /
        (mNumberOfValues[discreteVariable - (GetDimension() - GetNumberOfDiscreteVariable())] - 1);
      mCurrentDiscreteValueIndex[discreteVariable - GetNumberOfDiscreteVariable()] =
        previousNumber;
      mCurrentDiscreteValueIndex[discreteVariable - GetNumberOfDiscreteVariable()]++;
      value = mLeftBorder + d * mCurrentDiscreteValueIndex[discreteVariable -
        GetNumberOfDiscreteVariable()];
      return IProblem::OK;
    }
  }
  /// Проверяет является ли value допустимым значением для параметра с номером discreteVariable
  virtual bool IsPermissibleValue(double value, int discreteVariable)
  {
    if ((discreteVariable > GetDimension()) ||
      (discreteVariable < (GetDimension() - GetNumberOfDiscreteVariable())) ||
      (mNumberOfValues == 0))
      return false;
    double d = (mRightBorder - mLeftBorder) /
      (mNumberOfValues[discreteVariable - (GetDimension() - GetNumberOfDiscreteVariable())] - 1);
    double v = 0;
    for (int i = 0; i < mNumberOfValues[discreteVariable - (GetDimension() - GetNumberOfDiscreteVariable())]; i++)
    {
      v = mLeftBorder + d * i;
      if (fabs(v - value) < AccuracyDouble)
      {
        return true;
      }
    }
    return false;
  }
};

template <class Owner>
class Problem : public BaseProblem<Owner>
{
public:
  /** Метод возвращает координаты точки глобального минимума целевой функции
\param[out] y точка, в которой достигается оптимальное значение
\return Код ошибки (#OK или #UNDEFINED)
*/
  virtual int GetOptimumPoint(double* y) const
  {
    return IProblem::UNDEFINED;
  }

  /// Возвращает точку глобального оптимума для функции fNumber
  virtual int GetConstraintOptimumPoint(double* point, int fNumber)
  {
    return IProblem::UNDEFINED;
  }

  template <class ClassType, class Type>
  void ProblemParameterCheckSize2(ClassType& par, Type defVal, Type left, Type right, int size, bool isChackVal = false)
  {
    if (par.GetSize() < size)
    {
      if (par.GetSize() == 1)
      {
        Type d = defVal;
        if (par.GetSize() > 0)
          d = par[0];

        if (isChackVal)
          if ((d > right) || (d < left))
            d = defVal;
        par.SetSize(size);
        for (int i = 0; i < size; i++)
          par[i] = d;
      }
      else
      {
        int oldSize = par.GetSize();
        Type d = defVal;
        if (par.GetSize() > 0)
          d = par[0];

        if (isChackVal)
          if ((d > right) || (d < left))
            d = defVal;
        par.SetSize(size);
        for (int i = oldSize; i < size; i++)
          (par.GetData())[i] = d;
      }
    }
  }
};



class ProblemFromFunctionPointers : public Problem<ProblemFromFunctionPointers>
{
#undef OWNER_NAME
#define OWNER_NAME ProblemFromFunctionPointers
protected:
  std::vector<std::function<double(const double*)>> function;
  double optimumValue = 0;
  std::vector<double> optimumCoordinate;
  std::vector<double> lowerBounds;
  std::vector<double> upperBounds;
  bool isSetOptimum;
public:

  // ------------------------------------------------------------------------------------------------
  ProblemFromFunctionPointers(int dimention, std::vector<double> lower_, std::vector<double> upper_,
    std::vector<std::function<double(const double*)>> function_,
    bool isSetOptimum_ = false, double optimumValue_ = 0, double* optimumCoordinate_ = nullptr)
  {
    this->mDim = dimention;
    this->mOwner = this;
    this->mMinDimension = 1;
    this->mMaxDimension = 50;
    this->mNumberOfConstraints = 0;
    this->mLeftBorder = -1.0;
    this->mRightBorder = 1.0;
    this->mNumberOfCriterions = 1;
    function = function_;
    lowerBounds = lower_;
    upperBounds = upper_;
    isSetOptimum = isSetOptimum_;
    this->mNumberOfConstraints = function_.size() - this->mNumberOfCriterions;

    if (isSetOptimum)
    {
      optimumValue = optimumValue_;
      if (optimumCoordinate_ != nullptr)
      {
        optimumCoordinate.resize(dimention);
        for (int i = 0; i < dimention; i++)
          optimumCoordinate[i] = optimumCoordinate_[i];
      }
    }

  }

  // ------------------------------------------------------------------------------------------------
  /// Инициализация параметров
  void Init(int argc, char* argv[], bool isMPIInit)
  {
    mIsInit = false;
    BaseProblem<ProblemFromFunctionPointers>::Init(argc, argv, false);

    mIsInit = true;
  }

  // ------------------------------------------------------------------------------------------------
  int GetOptimumValue(double& value) const
  {
    if (!mIsInit || !isSetOptimum)
      return IProblem::UNDEFINED;

    value = optimumValue;

    return IProblem::OK;
  }

  // ------------------------------------------------------------------------------------------------
  int GetOptimumPoint(double* point) const
  {
    if (!mIsInit || !isSetOptimum)
      return IProblem::UNDEFINED;

    for (int i = 0; i < mDim; i++)
      point[i] = optimumCoordinate[i];

    return IProblem::OK;
  }
  /** Метод возвращает границы области поиска
*/
  virtual void GetBounds(double* lower, double* upper)
  {
    if (this->mIsInit)
      for (int i = 0; i < Dimension; i++)
      {
        lower[i] = lowerBounds[i];
        upper[i] = upperBounds[i];
      }
  }

  virtual double CalculateFunctionals(const double* x, int fNumber)
  {
    if (fNumber >= function.size())
      throw EXCEPTION("Error function number");
    return function[fNumber](x);
  }

};


#ifdef _GLOBALIZER_BENCHMARKS
#include "IGlobalOptimizationProblem.h"

class GlobalizerBenchmarksProblem : public Problem<GlobalizerBenchmarksProblem>
{
#undef OWNER_NAME
#define OWNER_NAME GlobalizerBenchmarksProblem

protected:
  IGlobalOptimizationProblem* problem;

  std::vector< std::vector<std::string>> DiscreteVariableValues;

  std::vector<double> A;
  std::vector<double> B;
  std::vector<double> optPoint;
  //std::vector<double> discreteValues;
  std::vector <std::vector<double>> discreteValues;


  void XtoYU(const double* x, std::vector<double>& y, std::vector<std::string>& u) const
  {
    y.resize(problem->GetNumberOfContinuousVariable());
    int i = 0;
    for (; i < problem->GetNumberOfContinuousVariable(); i++)
    {
      y[i] = x[i];
    }
    for (int j = 0; j < problem->GetNumberOfDiscreteVariable(); j++, i++)
    {
      if (x[i] < 0 || x[i] >  DiscreteVariableValues[j].size())
        throw - 1;
      else
      {
        u[j] = DiscreteVariableValues[j][int(x[i])];
      }
    }
  }

  void YUtoX(std::vector<double>& y, std::vector<std::string>& u, double* x) const
  {
    int i = 0;
    for (; i < problem->GetNumberOfContinuousVariable(); i++)
    {
      x[i] = y[i];
    }
    for (int j = 0; j < problem->GetNumberOfDiscreteVariable(); j++, i++)
    {
      bool f = false;
      for (int k = 0; k < DiscreteVariableValues[j].size(); k++)
      {
        if (DiscreteVariableValues[j][k] == u[j])
        {
          f = true;
          x[i] = k;
          break;
        }
      }
      if (!f)
        throw - 1;
    }
  }


public:
  GlobalizerBenchmarksProblem(IGlobalOptimizationProblem* _problem) : problem(_problem)
  {
    mDim = problem->GetDimension();
    mMaxDimension = MaxDim;
    mMinDimension = 1;
    mNumberOfCriterions = problem->GetNumberOfCriterions();
    mNumberOfConstraints = problem->GetNumberOfConstraints();

    NumberOfDiscreteVariable = problem->GetNumberOfDiscreteVariable();

    std::vector<double> x(this->mDim);
    A.resize(this->mDim);
    B.resize(this->mDim);
    problem->GetBounds(A, B);

    std::vector< std::vector<std::string>> values;
    problem->GetDiscreteVariableValues(DiscreteVariableValues);
    if (NumberOfDiscreteVariable > 0)
      mDefNumberOfValues = DiscreteVariableValues[0].size();

    mNumberOfValues = new int[GetNumberOfDiscreteVariable()];
    discreteValues.resize(GetNumberOfDiscreteVariable());
    for (int i = 0; i < GetNumberOfDiscreteVariable(); i++)
    {
      mNumberOfValues[i] = DiscreteVariableValues[i].size();
      discreteValues[i].resize(mNumberOfValues[i]);
      for (int j = 0; j < mNumberOfValues[i]; j++)
      {
        discreteValues[i][j] = j;
      }
    }

    mLeftBorder = A[0];
    mRightBorder = B[0];

    //discreteValues.resize(GetNumberOfDiscreteVariable());
    ////double d = (mRightBorder - mLeftBorder) / (mDefNumberOfValues - 1);
    //for (int i = 0; i < mDefNumberOfValues; i++)
    //{
    //  discreteValues[i] = mLeftBorder + d * i;
    //}

  }

  virtual int GetOptimumValue(double& value) const
  {
    return problem->GetOptimumValue(value);
  }

  virtual int GetOptimumPoint(double* x) const
  {
    std::vector<double> y(problem->GetNumberOfContinuousVariable());
    std::vector<std::string> u(problem->GetNumberOfContinuousVariable());

    int err = problem->GetOptimumPoint(y, u);

    YUtoX(y, u, x);

    return err;
  }

  virtual double CalculateFunctionals(const double* x, int fNumber)
  {
    std::vector<double> y(problem->GetNumberOfContinuousVariable());
    std::vector<std::string> u(problem->GetNumberOfContinuousVariable());

    XtoYU(x, y, u);

    return problem->CalculateFunctionals(y, u, fNumber);
  }


  virtual int SetConfigPath(const std::string& configPath)
  {
    return problem->SetConfigPath(configPath);
  }

  virtual int SetDimension(int dimension)
  {
    int err = problem->SetDimension(dimension);
    Dimension = problem->GetDimension();
    return err;
  }
  virtual int GetDimension() const
  {
    return problem->GetDimension();
  }


  virtual void GetBounds(double* lower, double* upper)
  {
    std::vector<double> lower_(mDim);
    std::vector<double> upper_(mDim);
    problem->GetBounds(lower_, upper_);
    for (int i = 0; i < mDim; i++)
    {
      lower[i] = lower_[i];
      upper[i] = upper_[i];

    }
  }

  /**
Возвращает число значений дискретного параметра discreteVariable.
GetDimension() возвращает общее число параметров.
(GetDimension() - GetNumberOfDiscreteVariable()) - номер начальной дискретной переменной
Для не дискретных переменных == -1
*/
  virtual int GetNumberOfValues(int discreteVariable)
  {
    if ((discreteVariable > GetDimension()) ||
      (discreteVariable < (GetDimension() - GetNumberOfDiscreteVariable())))
      return -1;
    if (mNumberOfValues == 0)
      return -1;
    return mNumberOfValues[discreteVariable - (GetDimension() - GetNumberOfDiscreteVariable())];
  }

  /// Очищает номер текущего значения для дискретного параметра
  virtual void ClearCurrentDiscreteValueIndex(int** mCurrentDiscreteValueIndex)
  {
    if (NumberOfDiscreteVariable > 0)
    {
      if (*mCurrentDiscreteValueIndex != 0)
        delete[] * mCurrentDiscreteValueIndex;

      *mCurrentDiscreteValueIndex = new int[NumberOfDiscreteVariable];
      for (int i = 0; i < NumberOfDiscreteVariable; i++)
        (*mCurrentDiscreteValueIndex)[i] = 0;
    }
  }

  /**
  Определяет значения дискретного параметра с номером discreteVariable
  Возвращает код ошибки.
  \param[out] values массив, в который будут сохранены значения дискретного параметра
  */
  virtual int GetAllDiscreteValues(int discreteVariable, double* values)
  {
    if ((discreteVariable > GetDimension()) ||
      (discreteVariable < (GetDimension() - GetNumberOfDiscreteVariable())))
      return IIntegerProgrammingProblem::ERROR_DISCRETE_VALUE;
    int* mCurrentDiscreteValueIndex = 0;
    ClearCurrentDiscreteValueIndex(&mCurrentDiscreteValueIndex);

    // сбрасываем значение индекса текущего значения и задаем левую границу
    GetNextDiscreteValues(mCurrentDiscreteValueIndex, values[0], discreteVariable, -1);
    int numVal = GetNumberOfValues(discreteVariable);
    // определяем все остальные значения
    for (int i = 1; i < numVal; i++)
    {
      GetNextDiscreteValues(mCurrentDiscreteValueIndex, values[i], discreteVariable);
    }
    return IProblem::OK;
  }
  /**
  Определяет значения дискретного параметра с номером discreteVariable после номера previousNumber
  Возвращает код ошибки.
  \param[in] previousNumber - номер значения после которого возвращается значение
  -2 - значение по умолчанию, возвращает следующее значение
  -1 - возвращает после -1, т.е. левую границу области
  \param[out] value переменная в которую сохраняется значение дискретного параметра
  */
  virtual int GetNextDiscreteValues(int* mCurrentDiscreteValueIndex, double& value, int discreteVariable, int previousNumber = -2)
  {
    if ((discreteVariable > GetDimension()) ||
      (discreteVariable < (GetDimension() - GetNumberOfDiscreteVariable())) ||
      (mCurrentDiscreteValueIndex == 0) ||
      (mNumberOfValues == 0))
      return IIntegerProgrammingProblem::ERROR_DISCRETE_VALUE;
    // если -1 то сбрасываем значение текущего номера
    if (previousNumber == -1)
    {
      mCurrentDiscreteValueIndex[discreteVariable - GetNumberOfDiscreteVariable()] = 0;
      value = 0;
      return IProblem::OK;
    }
    else if (previousNumber == -2)
    {      
      mCurrentDiscreteValueIndex[discreteVariable - GetNumberOfDiscreteVariable()]++;
      value = mCurrentDiscreteValueIndex[discreteVariable - GetNumberOfDiscreteVariable()];
      return IProblem::OK;
    }
    else
    {     
      mCurrentDiscreteValueIndex[discreteVariable - GetNumberOfDiscreteVariable()] =
        previousNumber;
      mCurrentDiscreteValueIndex[discreteVariable - GetNumberOfDiscreteVariable()]++;
      value = mCurrentDiscreteValueIndex[discreteVariable -
        GetNumberOfDiscreteVariable()];
      return IProblem::OK;
    }
  }

};
#endif


#endif
// - end of file ----------------------------------------------------------------------------------
