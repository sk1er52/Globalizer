#pragma once

#include "Problem.h"


/** Класс задач принимающий функции как параметры
*/
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
/** Класс задачи обертки для интерфейса IGlobalOptimizationProblem
*/
class GlobalizerBenchmarksProblem : public Problem<GlobalizerBenchmarksProblem>
{
#undef OWNER_NAME
#define OWNER_NAME GlobalizerBenchmarksProblem

protected:
  /// Новая задача
  IGlobalOptimizationProblem* problem;

  /// Все комбинации дискретных параметров
  std::vector< std::vector<std::string>> DiscreteVariableValues;
  /// Левая граница области определения
  std::vector<double> A;
  /// Правая граница области определения
  std::vector<double> B;
  /// Известная точка оптимума
  std::vector<double> optPoint;
  /// Номера значений дискретных параметров
  std::vector <std::vector<double>> discreteValues;

  /** Метод, получает по координатам непрерывных параметров и номерам дискретных два вектора с непрерывными и дискретными параметрами соответственно

  \param[in] x координаты точки, сначала непрерывные координаты, затем номера дискретных параметров
  \param[out] y непрерывные координаты точки, в которой необходимо вычислить значение
  \param[out] u целочисленые координаты точки, в которой необходимо вычислить значение
  */
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

  /** Метод, преобразует вектора с непрерывными и дискретными параметрами объединенный массив с координатами

  \param[in] y непрерывные координаты точки, в которой необходимо вычислить значение
  \param[in] u целочисленые координаты точки, в которой необходимо вычислить значение
  \param[out] x координаты точки, сначала непрерывные координаты, затем номера дискретных параметров
  */
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

  }
  /** Метод возвращает значение целевой функции в точке глобального минимума
  \param[out] value оптимальное значение
  \return Код ошибки (#OK или #UNDEFINED)
  */
  virtual int GetOptimumValue(double& value) const
  {
    return problem->GetOptimumValue(value);
  }
  /** Метод возвращает координаты точки глобального минимума целевой функции
  \param[out] x точка, в которой достигается оптимальное значение
  \return Код ошибки (#OK или #UNDEFINED)
  */
  virtual int GetOptimumPoint(double* x) const
  {
    std::vector<double> y(problem->GetNumberOfContinuousVariable());
    std::vector<std::string> u(problem->GetNumberOfContinuousVariable());

    int err = problem->GetOptimumPoint(y, u);

    YUtoX(y, u, x);

    return err;
  }
  /** Метод, вычисляющий функции задачи

  \param[in] x Точка, в которой необходимо вычислить значение
  \param[in] fNumber Номер вычисляемой функции. 0 соответствует первому ограничению,
  #GetNumberOfFunctions() - 1 -- последнему критерию
  \return Значение функции с указанным номером
  */
  virtual double CalculateFunctionals(const double* x, int fNumber)
  {
    std::vector<double> y(problem->GetNumberOfContinuousVariable());
    std::vector<std::string> u(problem->GetNumberOfContinuousVariable());

    XtoYU(x, y, u);

    return problem->CalculateFunctionals(y, u, fNumber);
  }

  /** Задание пути до конфигурационного файла

  Данный метод должн вызываться перед #Initialize
  \param[in] configPath строка, содержащая путь к конфигурационному файлу задачи
  \return Код ошибки
  */
  virtual int SetConfigPath(const std::string& configPath)
  {
    return problem->SetConfigPath(configPath);
  }

  /** Метод задаёт размерность задачи

  Данный метод должен вызываться перед #Initialize. Размерность должна быть в
  списке поддерживаемых.
  \param[in] dimension размерность задачи
  \return Код ошибки
  */
  virtual int SetDimension(int dimension)
  {
    int err = problem->SetDimension(dimension);
    Dimension = problem->GetDimension();
    return err;
  }

  ///Возвращает размерность задачи, можно вызывать после #Initialize
  virtual int GetDimension() const
  {
    return problem->GetDimension();
  }

  /** Метод возвращает границы области поиска
  */
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
