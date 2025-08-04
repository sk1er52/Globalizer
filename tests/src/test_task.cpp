#include "Task.h"
#include "Problem.h"
#include "test_config.h"

#include <gtest/gtest.h>
#include <string>
#include <cstdlib>

/**
  Вспомогательный класс, помогающий задать начальную конфигурацию объекта класса #Task,
  которая будет использоваться в тестах
 */
class TaskTest : public ::testing::Test
{
protected:
  /// Размерность задачи
  static const int n = 5;
  /// Размерность подзадачи
  static const int freeN = 2;
  /// Число функционалов
  static const int numOfFunc = 1;
  /// Левая граница области поиска
  double A[MaxDim];
  /// Правая граница области поиска
  double B[MaxDim];
  /// Указатель на задачу
  Task* task;
  //Функционалы
  IProblem* problem;

  void SetUp()
  {

    
    problem = new ProblemFromFunctionPointers(n, // размерность задачи
      std::vector<double>(parameters.Dimension, -2.2), // верхняя граница
      std::vector<double>(parameters.Dimension, 1.8), // нижняя граница
      std::vector<std::function<double(const double*)>>(1, [](const double* y)
        {
          double M_PI = 3.14159265358979323846;
          double sum = 0.;
          for (int j = 0; j < parameters.Dimension; j++)
            sum += y[j] * y[j] - 10. * cos(2.0 * M_PI * y[j]) + 10.0;
          return sum;
        }), // критерий
      true, // определен ли оптимум
      0, // значение глобального оптимума
      std::vector<double>(parameters.Dimension, 0).data() // координаты глобального минимума

    );
    parameters.Dimension = 5;
    task = new Task(problem, 0);

  }

  void TearDown()
  {
    delete task;
  }
};

/**
 * Проверка параметра размерности задачи N
 * 1<= N <= MaxDim
 */
TEST_F(TaskTest, throws_when_create_with_negative_N)
{
  int oldN = parameters.Dimension;
  parameters.Dimension = -1;
  ASSERT_ANY_THROW(Task testTask( problem, 0));
  parameters.Dimension = oldN;
}

TEST_F(TaskTest, throws_when_create_with_null_N)
{
  int oldN = parameters.Dimension;
  parameters.Dimension = 0;
  ASSERT_ANY_THROW(Task testTask( problem, 0));
  parameters.Dimension = oldN;
}

TEST_F(TaskTest, throws_when_create_with_too_large_N)
{
  int oldN = parameters.Dimension;
  parameters.Dimension = 500;
  ASSERT_ANY_THROW(Task testTask(problem, 0));
  parameters.Dimension = oldN;
}

/**
 * Создание задачи с корректными входными параметрами
 */
TEST_F(TaskTest, can_create_with_correct_values)
{
  ASSERT_NO_THROW(Task testTask( problem, 0));
}

