#ifdef _GLOBALIZER_BENCHMARKS
#include "IOptProblem.hpp"
#include "Problem.h"

/** Класс задачи обертки для интерфейса IOptProblem 
*/
class GlobalizerOptProblem : public Problem<GlobalizerOptProblem>
{
#undef OWNER_NAME
#define OWNER_NAME GlobalizerOptProblem

protected:
  /// Внутренняя задача оптимизации
  IOptProblem* problem;

  /// Левая граница области определения
  std::vector<double> A;
  /// Правая граница области определения
  std::vector<double> B;
  /// Известная точка глобального минимума
  std::vector<double> optPoint;

public:

  /** Конструктор обертки
  \param[in] _problem указатель на задачу типа IOptProblem
  */
  GlobalizerOptProblem(IOptProblem* _problem) : problem(_problem)
  {
    int dim = problem->GetDimension();

    mDim = dim;
    mMaxDimension = MaxDim;
    mMinDimension = 1;
    A.resize(this->mDim);
    B.resize(this->mDim);
    problem->GetBounds(A, B);

    mNumberOfCriterions = 1;
    mNumberOfConstraints = 0;

    optPoint = problem->GetOptimumPoint();

    mLeftBorder = A[0];
    mRightBorder = B[0];
  }

  /** Метод возвращает значение целевой функции в точке глобального минимума
  \param[out] value оптимальное значение
  \return Код ошибки (#OK или #UNDEFINED)
  */
  virtual int GetOptimumValue(double& value) const override
  {
    value = problem->GetOptimumValue();
    return IProblem::OK;
  }

  /** Метод возвращает координаты точки глобального минимума целевой функции
  \param[out] x точка, в которой достигается оптимальное значение
  \return Код ошибки (#OK или #UNDEFINED)
  */
  virtual int GetOptimumPoint(double* x) const override
  {
  std::vector<double> y = problem->GetOptimumPoint();
    if (static_cast<int>(y.size()) != mDim)
        return IProblem::UNDEFINED;

    for (int i = 0; i < mDim; ++i)
        x[i] = y[i];

        return IProblem::OK;
  }

  /** Метод, вычисляющий функции задачи
  \param[in] x Точка, в которой необходимо вычислить значение
  \param[in] fNumber Номер вычисляемой функции. 0 соответствует целевой функции, другие значения не поддерживаются.
  \return Значение функции с указанным номером
  */
  virtual double CalculateFunctionals(const double* x, int fNumber) override
  {
    if (fNumber != 0)
        return std::numeric_limits<double>::quiet_NaN();

    std::vector<double> y(x, x + mDim);
    return problem->ComputeFunction(y);
  }

  /** Задание пути до конфигурационного файла

  Данный метод должен вызываться перед #Initialize.
  \param[in] configPath строка, содержащая путь к конфигурационному файлу задачи
  \return Код ошибки
  */
  virtual int SetConfigPath(const std::string& configPath) override
  {
    (void)configPath;
    return IProblem::UNDEFINED;
  }

  /** Метод задаёт размерность задачи

  Данный метод должен вызываться перед #Initialize. Размерность должна быть в списке поддерживаемых.
  \param[in] dimension размерность задачи
  \return Код ошибки
  */
  virtual int SetDimension(int dimension) override
  {
    if (dimension != problem->GetDimension())
        return IProblem::ERROR;
    Dimension = dimension;
    return IProblem::OK;
  }

  /// Возвращает размерность задачи, можно вызывать после #Initialize
  virtual int GetDimension() const override
  {
      return problem->GetDimension();
  }

  /** Метод возвращает границы области поиска
  */
  virtual void GetBounds(double* lower, double* upper) override
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
};

#endif 